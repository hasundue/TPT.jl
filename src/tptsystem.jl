type TPTSystem{Tr <: RefSystem, Tp <: PertSystem}
  ref::Tr # reference system
  pert::Tp # perturbing system
end

function prdf(sys::TPTSystem)
  if typeof(sys.ref) <: IndependentRefSystem
    prdf(sys.ref)
  else
    prdf(sys.ref, sys.pert)
  end
end

function psf(sys::TPTSystem)
  if typeof(sys.ref) <: IndependentRefSystem
    psf(sys.ref)
  else
    psf(sys.ref, sys.pert)
  end
end

function pairpotential(sys::TPTSystem)
  pairpotential(sys.pert)
end

function freeenergy(sys::TPTSystem)
  if typeof(sys.ref) <: IndependentRefSystem
    F₀, rep0 = freeenergy(sys.ref)
  else
    F₀, rep0 = freeenergy(sys.ref, sys.pert)
  end
  F₁, rep1 = freeenergy(sys.ref, sys.pert)
  return F₀ + F₁, rep0, rep1
end
