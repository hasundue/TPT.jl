type TPTSystem{Tr <: RefSystem, Tp <: PertSystem}
  m::Int # number of components
  c::Vector{Float64} # conposition (mole fraction)
  label::Vector{AbstractString} # labels of components
  ref::Tr # reference system
  pert::Tp # perturbing system
end

function TPTSystem(m::Int, c::Vector{Float64}, ref::RefSystem;
                   label=AbstractString[])
  return TPTSystem(m, c, label, ref, NoPert())
end

function TPTSystem(m::Int, c::Vector{Float64}, ref::RefSystem, pert::PertSystem;
                   label=AbstractString[])
  return TPTSystem(m, c, label, ref, pert)
end

function prdf(sys::TPTSystem)
  prdf(sys.ref)
end

function psf(sys::TPTSystem)
  psf(sys.pert)
end
