module custom_structures

export ParticleInstance1D, ParticleInstance2D, MarkedParticleInstance, AbstractParticleInstance, AbstractMarkedParticleInstance

using DifferentialEquations: OrdinaryDiffEq.ODEIntegrator
using DocStringExtensions

using Architectures: AbstractParticleInstance, AbstractMarkedParticleInstance

# ParticleInstance is the Stucture that carries each particle.
mutable struct ParticleInstance2D <: AbstractParticleInstance
        position_ij::Tuple{Int, Int}
        position_xy::Tuple{Float64, Float64}
        ODEIntegrator::ODEIntegrator
        boundary :: Bool
        on::Bool
end

mutable struct ParticleInstance1D <: AbstractParticleInstance
        position_ij::Int
        position_xy::Float64
        ODEIntegrator::ODEIntegrator
        boundary::Bool
        on::Bool
end

# Debugging ParticleInstance
mutable struct MarkedParticleInstance <: AbstractMarkedParticleInstance
        Particle::AbstractParticleInstance
        time :: Float64
        state :: Vector{Any}
        errorReturnCode
end

Base.copy(s::ParticleInstance1D) = ParticleInstance1D(s.position_ij, s.position_xy, s.ODEIntegrator, s.boundary, s.on)
Base.copy(s::ParticleInstance2D) = ParticleInstance2D(s.position_ij, s.position_xy, s.ODEIntegrator, s.boundary, s.on)



end