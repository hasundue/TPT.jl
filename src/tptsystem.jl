type TPTSystem{Tr <: RefSystem, Tp <: PertSystem}
  ref::Tr # reference system
  pert::Tp # perturbing system
end

function prdf(sys::TPTSystem)
  prdf(sys.ref)
end

function psf(sys::TPTSystem)
  psf(sys.pert)
end
