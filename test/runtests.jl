#
# Regression test
# (each julia script of test can also be run independently)
#

tests = [
  "unary.jl",
  "optimize_rc.jl"
]

for test in tests
  include(test)
end
