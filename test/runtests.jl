#
# Regression test
# (each julia script of test can also be run independently)
#

tests = [
  "unary.jl",
  "tm_optim.jl",
  "tm_binary.jl"
]

for test in tests
  @time include(test)
end
