#
# Regression test
# (each julia script of test can also be run independently)
#

cd(Pkg.dir("TPT"))

tests = [
  "unary.jl",
  "optimize_rc.jl"
]

for test in tests
  script = joinpath("test", test)
  include(script)
end
