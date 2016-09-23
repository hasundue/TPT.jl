type TPTSystem{Tr <: ReferenceSystem, Tp <: Perturbation}
  ref::Tr # reference system
  pert::Tp # perturbing system
end

function prdf(sys::TPTSystem)
  if typeof(sys.ref) <: IndependentReferenceSystem
    prdf(sys.ref)
  else
    prdf(sys.ref, sys.pert)
  end
end

function psf(sys::TPTSystem)
  if typeof(sys.ref) <: IndependentReferenceSystem
    psf(sys.ref)
  else
    psf(sys.ref, sys.pert)
  end
end

function pairpotential(sys::TPTSystem)
  pairpotential(sys.pert)
end

function freeenergy(sys::TPTSystem)
  if typeof(sys.ref) <: IndependentReferenceSystem
    F₀, rep0 = freeenergy(sys.ref)
  else
    F₀, rep0 = freeenergy(sys.ref, sys.pert)
  end
  F₁, rep1 = freeenergy(sys.ref, sys.pert)
  return F₀ + F₁, rep0, rep1
end
