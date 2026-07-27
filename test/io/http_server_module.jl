using Test
using HTTP
using JSON3
using OpenSEFA
using OpenSEFA.HTTPServerModule
using OpenSEFA.ReconciliationModule

const TEST_PORT = 18080

@testset "HTTPServerModule" begin
  start_server(port = TEST_PORT, async = true)
  sleep(1)  # Give server time to start.

  @testset "Health endpoint returns uptime and request stats." begin
    resp = HTTP.get("http://localhost:$TEST_PORT/health")
    @test resp.status == 200
    data = JSON3.read(resp.body)
    @test haskey(data, Symbol("uptime (sec)"))
    @test haskey(data, :requests)
    @test haskey(data[:requests], :total)
  end

  @testset "Solve endpoint returns a valid solution." begin
    payload = JSON3.write(Dict(
      "content" => Dict(
        "name" => "Test",
        "variables" => [
          Dict("name" => "F1", "type" => "MeasuredVariable", "value" => 1.0, "uncertainty" => 0.1, "is_transfer_coefficient" => false),
          Dict("name" => "F2", "type" => "UnmeasuredVariable", "is_transfer_coefficient" => false),
          Dict("name" => "F3", "type" => "MeasuredVariable", "value" => 2.0, "uncertainty" => 0.3, "is_transfer_coefficient" => false),
        ],
        "equations" => [Dict(
          "constant_term" => Dict("value" => 0),
          "linear_terms" => [
            Dict("name" => "F1", "factor" => 1),
            Dict("name" => "F3", "factor" => 1),
            Dict("name" => "F2", "factor" => -1),
          ],
          "bilinear_terms" => [],
        )],
        "solver" => "JuMP",
      ),
    ))

    resp = HTTP.post(
      "http://localhost:$TEST_PORT/solve",
      ["Content-Type" => "application/json"],
      payload,
    )
    @test resp.status == 200
    data = JSON3.read(resp.body)
    @test haskey(data, :status)
    @test haskey(data, :variables_value)
    # F1 + F3 = F2, so F2 ≈ 3.0.
    @test isapprox(data[:variables_value][:F2], 3.0; atol = 0.1)
  end

  @testset "CORS headers are present in response." begin
    resp = HTTP.get("http://localhost:$TEST_PORT/health")
    headers = Dict(resp.headers)
    @test get(headers, "Access-Control-Allow-Origin", "") == "*"
  end

  stop_server()
end
