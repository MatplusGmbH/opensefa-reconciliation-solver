using MaterialFlowAnalysis
using BenchmarkTools

path = joinpath(MODELS_PATH, "Cencic_2016")
strace = parse(STANTrace, path)
rec = convert(ReconciliationProblem, strace)

# Values calculated by STAN.
sol_stan = convert(ReconciliationSolution, strace)

# Different solvers.
BENCHMARK_SOLVERS = ReconciliationSolver[
  DefaultSolver(),
  DefaultSolver(; optimization_solver = MT24Solver()),
  MT24Solver(),
]
[solve(rec, solver) for solver in BENCHMARK_SOLVERS]
