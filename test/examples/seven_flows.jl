using Distributions: Normal
using OpenSEFA.ConstraintsModule
using OpenSEFA.STANModule: parse_equation
using OpenSEFA.ReconciliationModule
using OpenSEFA.ReconciliationModule: add_variable!

# Example with seven mass flows and three processes from [1].
# [1] Cencic, Oliver. "Nonlinear data reconciliation in material flow analysis with software STAN."
# Sustainable Environment Research 26, no. 6 (2016): 291-298.
function seven_flows_cencic_2012()
  # Define the equations based on the flowsheet example using string parsing.
  f1 = parse_equation("m1 + m2 + m4 - m3 = 0")
  f2 = parse_equation("m3 - m4 - m5 = 0")
  f3 = parse_equation("m5 - m6 - m7 = 0")
  f4 = parse_equation("m4 - m3 * tc34 = 0")

  # Define the measured variables with their uncertainties.
  m1 = MeasuredVariable(; name = "m1", value = 100.0, uncertainty = 10.0)
  m3 = MeasuredVariable(; name = "m3", value = 300.0, uncertainty = 30.0)
  m5 = MeasuredVariable(; name = "m5", value = 160.0, uncertainty = 16.0)
  # This is a transfer coefficient.
  tc34 = MeasuredVariable(; name = "tc34", value = 0.5, uncertainty = 0.05, tc = true)

  # Define the unknown variables (not measured).
  m4 = UnmeasuredVariable(; name = "m4")
  m6 = UnmeasuredVariable(; name = "m6")
  m7 = UnmeasuredVariable(; name = "m7")

  # Define the variables which do not have uncertainty (fixed).
  m2 = FixedVariable(; name = "m2", value = 50.0)

  # Instantiate the reconciliation problem.
  rec = ReconciliationProblem(; name = "Seven flows Cencic2016")
  add_variable!(rec, [m1, m2, m3, m4, m5, m6, m7, tc34])
  add_equation!(rec, [f1, f2, f3, f4])

  rec
end

# Bounded variant required for global optimizers (e.g. AlpineSolver) that rely on
# partition-based relaxations: all variables need finite upper bounds.
# Upper bounds are set to 10x the respective measured values.
function seven_flows_bounded_cencic_2012()
  f1 = parse_equation("m1 + m2 + m4 - m3 = 0")
  f2 = parse_equation("m3 - m4 - m5 = 0")
  f3 = parse_equation("m5 - m6 - m7 = 0")
  f4 = parse_equation("m4 - m3 * tc34 = 0")

  m1 = MeasuredVariable(; name = "m1", value = 100.0, uncertainty = 10.0, ub = 1000.0)
  m3 = MeasuredVariable(; name = "m3", value = 300.0, uncertainty = 30.0, ub = 3000.0)
  m5 = MeasuredVariable(; name = "m5", value = 160.0, uncertainty = 16.0, ub = 1600.0)
  tc34 = MeasuredVariable(; name = "tc34", value = 0.5, uncertainty = 0.05, tc = true)

  m4 = UnmeasuredVariable(; name = "m4", ub = 3000.0)
  m6 = UnmeasuredVariable(; name = "m6", ub = 1600.0)
  m7 = UnmeasuredVariable(; name = "m7", ub = 1600.0)

  m2 = FixedVariable(; name = "m2", value = 50.0)

  rec = ReconciliationProblem(; name = "Seven flows Cencic2016 (bounded)")
  add_variable!(rec, [m1, m2, m3, m4, m5, m6, m7, tc34])
  add_equation!(rec, [f1, f2, f3, f4])
  rec
end

function random_seven_flows_cencic_2012(n::Int, factor::Float64; rng)
  f1 = parse_equation("m1 + m2 + m4 - m3 = 0")
  f2 = parse_equation("m3 - m4 - m5 = 0")
  f4 = parse_equation("m4 - m3 * tc34 = 0")
  f3 = parse_equation("m5 - m6 - m7 = 0")

  map(1:n) do i
    m1_val0 = 100.0
    m3_val0 = 300.0
    m5_val0 = 160.0
    if i == 1
      m1_val = m1_val0
      m3_val = m3_val0
      m5_val = m5_val0
    else
      m1_val = rand(rng, Normal(m1_val0, 10.0 * factor))
      m3_val = rand(rng, Normal(m3_val0, 30.0 * factor))
      m5_val = rand(rng, Normal(m5_val0, 16.0 * factor))
    end
    m1 = MeasuredVariable(; name = "m1", value = m1_val, uncertainty = 10.0 * factor)
    m3 = MeasuredVariable(; name = "m3", value = m3_val, uncertainty = 30.0 * factor)
    m5 = MeasuredVariable(; name = "m5", value = m5_val, uncertainty = 16.0 * factor)
    tc34 = MeasuredVariable(; name = "tc34", value = 0.5, uncertainty = 0.05 * factor, tc = true)

    m4 = UnmeasuredVariable(; name = "m4")
    m6 = UnmeasuredVariable(; name = "m6")
    m7 = FixedVariable(; name = "m7", value = 0.0)
    m2 = FixedVariable(; name = "m2", value = 50.0)

    rec = ReconciliationProblem(; name = "Seven flows Cencic2016")
    add_variable!(rec, [m1, m2, m3, m4, m5, m6, m7, tc34])
    add_equation!(rec, [f1, f2, f3, f4])
    rec
  end
end
