
using ForwardDiff
using TaylorSeries
using IntervalArithmetic

f(x) = [sin(x[1]) + cos(x[2]), exp(x[1] + x[2])]

taylor_expand(f, [0.0, 0.0]; order=3)
# result
# 1.0 + 1.0 x₁ - 0.5 x₂² - 0.16666666666666666 x₁³ + 𝒪(‖x‖⁴)
# 1.0 + 1.0 x₁ + 1.0 x₂ + 0.5 x₁² + 1.0 x₁ x₂ + 0.5 x₂² + 0.16666666666666666 x₁³ + 0.5 x₁² x₂ + 0.5 x₁ x₂² + 0.16666666666666666 x₂³ + 𝒪(‖x‖⁴)

# taylor_expand does not work when x0 is an interval
#taylor_expand(f, [interval(0, 0.2), interval(0, 0.1)]; order=3)
