#
# Regression test
# (each julia script of test can also be run independently)
#

tests = [
  # Unit tests
  "unit_ahs.jl",
  "unit_wca.jl",
  "unit_lwca.jl",
  "unit_nfe.jl",
  "unit_nfetb.jl",
  "unit_hhtb.jl",

  # Investigational calculations
  "tm_optim.jl",
  "tm_binary.jl",
  "tm_energy.jl"
  "tm_hhtb.jl"
]

for test in tests
  @time include(test)
end
