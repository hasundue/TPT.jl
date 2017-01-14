abstract ReferenceSystem
abstract IndependentReferenceSystem <: ReferenceSystem
abstract DependentReferenceSystem <: ReferenceSystem

abstract Perturbation
abstract NFEPerturbation <: Perturbation
abstract TBPerturbation <: Perturbation

abstract PseudoPotential

typealias Pairpotential Function
