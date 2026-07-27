"""
HTTP server that provides REST API endpoints for data reconciliation problems.

# Features

- RESTful API for solving material flow analysis and substance flow analysis problems.
- Health endpoint for monitoring server status and performance.
- Cross-origin resource sharing (CORS) support.
- Request tracking and statistics.

# Endpoints

- `POST /solve`: Solves a data reconciliation problem:
  - Accepts JSON payload with `content` field.
  - Returns solution with objective value, variable values, and residual.

- `GET /health`: Returns server health information:
  - Uptime in seconds.
  - Request statistics (total, average time, max time).

# Usage

```julia
using OpenSEFA

# Start the server on default host (0.0.0.0) and port (8080)
start_server()

# The server can be configured via environment variables:
# - SEFA_SOLVER_HTTP_HOST: Host to bind to
# - SEFA_SOLVER_HTTP_PORT: Port to listen on
```

# Example Client Usage

```bash
# Check server health
curl http://localhost:8080/health

# Solve a model from content
curl -X POST http://localhost:8080/solve \\
  -H "Content-Type: application/json" \\
  -d '{"content": { ... model JSON ... }}'
```
(see the README.md for a complete minimal example).
"""
module HTTPServerModule

import Pkg

using Dates
using Dictionaries
using HTTP
using JSON3
using Oxygen
using Statistics
using StructTypes
@oxidize
using StatsBase

using ..ConstraintsModule
using ..DefaultSolverModule
using ..ReconciliationModule
using ..STANModule
using ..SerializationModule

export start_server, stop_server

# Server health tracking.
const SERVER_START_TIME = Ref{DateTime}(now()) # ms
const REQUEST_TIMES = Ref{Vector{Float64}}(Float64[])
const REQUEST_COUNT = Ref{Int}(0)
# Store only the last N request times.
const MAX_REQUEST_TIMES = 100

"Track request processing time."
function track_request_time(start_time)
  elapsed = (time() - start_time) # s
  REQUEST_COUNT[] += 1

  # Store processing time, keeping only the most recent ones.
  times = REQUEST_TIMES[]
  push!(times, elapsed)
  if length(times) > MAX_REQUEST_TIMES
    popfirst!(times)
  end
  REQUEST_TIMES[] = times
end

"Struct holding either the model's content."
struct Model
  "Actual payload with the model content."
  content::Dict{String, Any}
end

StructTypes.StructType(::Type{Model}) = StructTypes.Struct()

const CORS_HEADERS = [
  "Access-Control-Allow-Origin" => "*",
  "Access-Control-Allow-Headers" => "*",
  "Access-Control-Allow-Methods" => "POST, GET, OPTIONS",
]

"CORS Middleware to handle preflight requests and add CORS headers."
function CorsMiddleware(handler)
  return function (req::HTTP.Request)
    if HTTP.method(req) == "OPTIONS"
      # Handle pre-flight requests/
      return HTTP.Response(200, CORS_HEADERS)
    else
      # Handle actual request and add CORS headers to the response.
      try
        response = handler(req)
        # Add CORS headers to the actual response
        push!(response.headers, CORS_HEADERS...)

        return response
      catch e
        # If the handler throws an error, return a 500 response
        println("Error in handler: ", e)
        # Optionally include stacktrace: showerror(stdout, e, catch_backtrace())
        error_response = HTTP.Response(
          500,
          CORS_HEADERS,
          body = JSON3.write(Dict("error" => "Internal Server Error", "details" => sprint(showerror, e))),
        )
        push!(error_response.headers, "Content-Type" => "application/json")
        return error_response
      end
    end
  end
end

# Define route for solving the model.
@post "/solve" function (req::HTTP.Request, json_model::Json{Model})
  start_time = time()

  # Extract model from JSON payload.
  model = json_model.payload

  # Load model data as a Julia dictionary.
  content = model.content

  # Load the problem as a Julia struct.
  rec = load_model(content)

  # Parse which solver should be used.
  solver = load_solver(content)

  # Solve the problem. For better user experience:
  # - Ensure that if the basic constraints are not satisfied, the result status is infeasible.
  # - Replace variable values with their upper (or lower) bound when numerical error is present, eg. return `0.0` instead of `-1.3e-9`.
  sol = solve(rec, solver; force_basic_constraints = true, round_results = true)

  # Return solution as JSON.
  result = to_dict(sol)

  track_request_time(start_time)
  result
end

# Define route for health check.
@get "/health" function (req::HTTP.Request)
  # Calculate request statistics if any
  req_times = REQUEST_TIMES[]
  avg_request_time = isempty(req_times) ? 0.0 : round(mean(req_times), sigdigits = 2) # s
  max_request_time = isempty(req_times) ? 0.0 : round(maximum(req_times), sigdigits = 2) # s

  # Calculate uptime
  uptime_seconds = Dates.value(now() - SERVER_START_TIME[]) / 1000

  # Return health information as Dict (will be converted to JSON).
  Dict(
    "uptime (sec)" => round(uptime_seconds, sigdigits = 2),
    "requests" => Dict(
      "total" => REQUEST_COUNT[],
      "average time (sec)" => avg_request_time,
      "max time (sec)" => max_request_time,
    ),
  )
end

"Start the web server with CORS middleware applied."
function start_server(;
  host = haskey(ENV, "SEFA_SOLVER_HTTP_HOST") ? ENV["SEFA_SOLVER_HTTP_HOST"] : "0.0.0.0",
  port = haskey(ENV, "SEFA_SOLVER_HTTP_PORT") ? parse(Int, ENV["SEFA_SOLVER_HTTP_PORT"]) : 8080,
  async = false,
)
  # Reset start time when server starts.
  SERVER_START_TIME[] = now() # ms

  serve(
    middleware = [CorsMiddleware],
    async = async,
    host = host,
    port = port,
  )
end

"Stop the web server."
stop_server() = terminate()

end
