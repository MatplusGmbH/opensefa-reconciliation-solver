# [HTTP API](@id http_api)

OpenSEFA.jl provides a RESTful HTTP API for solving data reconciliation problems. This enables integration with web applications, remote clients, and other programming languages.

## Starting the Server

To start the HTTP server, run:

```bash
$ julia --project -e "using OpenSEFA; start_server()"
```

If successful, you should see output similar to:

```
   ____
  / __ \_  ____  ______ ____  ____
 / / / / |/_/ / / / __ `/ _ \/ __ \
/ /_/ />  </ /_/ / /_/ /  __/ / / /
\____/_/|_|\__, /\__, /\___/_/ /_/
          /____//____/

[ Info: 📦 Version 1.7.1 (2025-02-06)
[ Info: ✅ Started server: http://0.0.0.0:8080
[ Info: 📖 Documentation: http://0.0.0.0:8080/docs
[ Info: 📊 Metrics: http://0.0.0.0:8080/docs/metrics
[ Info: Listening on: 0.0.0.0:8080, thread id: 1
```

!!! note "Server Configuration"
    By default, the server listens on `0.0.0.0:8080`. You can customize the host and port in the server configuration.

## API Endpoints

### Health Check

**Endpoint:** `GET /health`

Check the server's operational status, including uptime and request statistics.

**Example Request:**

```bash
curl http://localhost:8080/health
```

**Example Response:**

```json
{
  "uptime (seconds)": 45.68,
  "requests": {
    "total": 12,
    "max time (seconds)": 2.34,
    "average time (seconds)": 0.87
  }
}
```

**Response Fields:**
- `uptime (seconds)`: Time since server started
- `requests.total`: Total number of requests processed
- `requests.max time (seconds)`: Maximum request processing time
- `requests.average time (seconds)`: Average request processing time

### Solve Problem

**Endpoint:** `POST /solve`

Solve a data reconciliation problem specified either by file path or JSON content.

#### Method 1: Solve from File Path

Send the path to a JSON file containing the problem definition.

**Example Request:**

```bash
curl -X POST http://localhost:8080/solve \
  -H "Content-Type: application/json" \
  -d '{
    "filepath": "/path/to/model/seven_flows.json"
  }'
```

**Example Response:**

```json
{
  "status": "success",
  "solution": {
    "objective_value": 0.2959,
    "variables": {
      "F1_mass_N1": 1.12,
      "F2_mass_N1": 3.45,
      "F3_mass_N1": 2.33,
      ...
    },
    "max_equation_discrepancy": 1.026e-9,
    "max_fixed_discrepancy": 0.0
  }
}
```

#### Method 2: Solve from JSON Content

Send the complete problem definition in the request body.

**Example Request:**

```bash
curl -X POST http://localhost:8080/solve \
  -H "Content-Type: application/json" \
  -d '{
    "content": {
      "name": "Simple Example",
      "variables": [
        {
          "name": "F1_mass_N1",
          "type": "MeasuredVariable",
          "value": 1.0,
          "uncertainty": 0.1,
          "is_transfer_coefficient": false
        },
        {
          "name": "F2_mass_N1",
          "type": "UnmeasuredVariable",
          "is_transfer_coefficient": false
        },
        {
          "name": "F3_mass_N1",
          "type": "MeasuredVariable",
          "value": 2.0,
          "uncertainty": 0.3,
          "is_transfer_coefficient": false
        }
      ],
      "equations": [
        {
          "constant_term": {"value": 0},
          "linear_terms": [
            {"name": "F1_mass_N1", "factor": 1},
            {"name": "F3_mass_N1", "factor": 1},
            {"name": "F2_mass_N1", "factor": -1}
          ],
          "bilinear_terms": []
        }
      ],
      "solver": "JuMP"
    }
  }'
```

**Response Format:**

Success response (HTTP 200):
```json
{
  "status": "success",
  "solution": {
    "name": "Problem name",
    "objective_value": 0.123,
    "variables": {
      "var1": 1.23,
      "var2": 4.56,
      ...
    },
    "max_equation_discrepancy": 1.0e-9,
    "max_fixed_discrepancy": 0.0,
    "num_equations": 5,
    "num_variables": 8
  }
}
```

Error response (HTTP 400/500):
```json
{
  "status": "error",
  "message": "Error description",
  "details": "Additional error details (optional)"
}
```

## JSON Problem Format

The JSON format for defining problems follows this schema:

```json
{
  "name": "Problem Name",
  "variables": [
    {
      "name": "variable_name",
      "type": "MeasuredVariable|UnmeasuredVariable|FixedVariable",
      "value": 1.0,                    // Required for Measured/Fixed
      "uncertainty": 0.1,               // Required for Measured
      "is_transfer_coefficient": false  // Optional, default: false
    }
  ],
  "equations": [
    {
      "constant_term": {"value": 0.0},
      "linear_terms": [
        {"name": "variable_name", "factor": 1.0}
      ],
      "bilinear_terms": [               // Optional
        {
          "row_variable": "var1",
          "col_variable": "var2",
          "factor": 0.5
        }
      ]
    }
  ],
  "solver": "JuMP|Default|NonlinearSolve"  // Optional, default: "Default"
}
```

### Variable Types

- **MeasuredVariable**: Known value with uncertainty (requires `value` and `uncertainty`)
- **UnmeasuredVariable**: Unknown value to be determined by optimization
- **FixedVariable**: Known value without uncertainty (requires `value`)

### Equation Structure

Each equation represents a constraint of the form:

```
constant_term + sum(linear_terms) + sum(bilinear_terms) = 0
```

- **linear_terms**: Terms of the form `factor * variable`
- **bilinear_terms**: Terms of the form `factor * var1 * var2`

## Integration with OpenSEFA editor

OpenSEFA.jl can be integrated with the OpenSEFA editor. If you are interested in the editor source code, please [contact](mailto:contact@matplus.eu) [Matplus GmbH](https://www.matplus.eu/). A valid JointJS+ license must be held.

### Setup

1. Switch to the OpenSEFA editor project
2. Create a `config.js` file in the project root:

```javascript
window.__ENV__ = {
  EDA_HOST: "http://0.0.0.0:8080",
  SOLVER_TOKEN: ""
};
```

3. Build the editor: `npm run build`

### Workflow

1. Start the OpenSEFA.jl HTTP server
2. Launch the editor: `npm run dev`
3. Create or load a model in the editor
4. Click the "Solve" button to send a request to the Julia backend
5. View results in the editor interface

!!! tip "CORS and Security"
    When deploying in production, ensure proper CORS configuration and authentication tokens are set.

## Client Examples

### Python Client

```python
import requests
import json

# Define the problem
problem = {
    "content": {
        "name": "Python Example",
        "variables": [
            {"name": "F1", "type": "MeasuredVariable", "value": 1.0, "uncertainty": 0.1},
            {"name": "F2", "type": "UnmeasuredVariable"},
            {"name": "F3", "type": "MeasuredVariable", "value": 2.0, "uncertainty": 0.3}
        ],
        "equations": [
            {
                "constant_term": {"value": 0},
                "linear_terms": [
                    {"name": "F1", "factor": 1},
                    {"name": "F3", "factor": 1},
                    {"name": "F2", "factor": -1}
                ],
                "bilinear_terms": []
            }
        ]
    }
}

# Send request
response = requests.post(
    "http://localhost:8080/solve",
    headers={"Content-Type": "application/json"},
    data=json.dumps(problem)
)

# Parse response
if response.status_code == 200:
    result = response.json()
    print(f"Objective: {result['solution']['objective_value']}")
    print(f"Variables: {result['solution']['variables']}")
else:
    print(f"Error: {response.json()['message']}")
```

### JavaScript Client

```javascript
async function solveProblem() {
  const problem = {
    content: {
      name: "JavaScript Example",
      variables: [
        { name: "F1", type: "MeasuredVariable", value: 1.0, uncertainty: 0.1 },
        { name: "F2", type: "UnmeasuredVariable" },
        { name: "F3", type: "MeasuredVariable", value: 2.0, uncertainty: 0.3 }
      ],
      equations: [
        {
          constant_term: { value: 0 },
          linear_terms: [
            { name: "F1", factor: 1 },
            { name: "F3", factor: 1 },
            { name: "F2", factor: -1 }
          ],
          bilinear_terms: []
        }
      ]
    }
  };

  try {
    const response = await fetch("http://localhost:8080/solve", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(problem)
    });

    const result = await response.json();

    if (response.ok) {
      console.log("Objective:", result.solution.objective_value);
      console.log("Variables:", result.solution.variables);
    } else {
      console.error("Error:", result.message);
    }
  } catch (error) {
    console.error("Request failed:", error);
  }
}

solveProblem();
```

## Error Handling

The API returns appropriate HTTP status codes:

- **200 OK**: Request successful, solution found
- **400 Bad Request**: Invalid request format or problem specification
- **500 Internal Server Error**: Solver failed or internal error

Error responses include a `message` field with details:

```json
{
  "status": "error",
  "message": "Invalid variable type: UnknownType"
}
```

## Performance Considerations

- The server processes requests synchronously
- For large models, requests may take several seconds to minutes
- Consider implementing client-side timeouts (recommended: 60-300 seconds)
- Monitor the `/health` endpoint for server performance metrics

## Next Steps

- Review [Examples](@ref examples_page) for problem formulations
- Check the API Reference for programmatic Julia usage
- See [Getting Started](@ref getting_started) for installation details
