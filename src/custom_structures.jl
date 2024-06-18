module custom_structures

export MarkedParticleInstance, AbstractParticleInstance, AbstractMarkedParticleInstance, ParticleInstance, wni

# for debuggind
export ParticleInstance1D, ParticleInstance2D, ParticleInstance2DLayer

using DifferentialEquations: OrdinaryDiffEq.ODEIntegrator
using DocStringExtensions
using StaticArrays

using ..Architectures: AbstractParticleInstance, AbstractMarkedParticleInstance


# ParticleInstance is the Stucture that carries each particle.

# ParticleInstance for multi-layer model
mutable struct ParticleInstance2DLayer <: AbstractParticleInstance
        position_ij::Tuple{Int,Int,Int}
        position_xy::Tuple{Float64,Float64}
        ODEIntegrator::ODEIntegrator
        boundary::Bool
        on::Bool
end

# ParticleInstance for a single layer model
mutable struct ParticleInstance2D <: AbstractParticleInstance
        position_ij::Tuple{Int, Int}
        position_xy::Tuple{Float64, Float64}
        ODEIntegrator::ODEIntegrator
        boundary :: Bool
        on::Bool
end

# 1D model 
mutable struct ParticleInstance1D <: AbstractParticleInstance
        position_ij::Int
        position_xy::Float64
        ODEIntegrator::ODEIntegrator
        boundary::Bool
        on::Bool
end

# multiple dispatch for ParticleInstance
ParticleInstance(ij::Int, xy::Float64, ODEI::ODEIntegrator, boundary::Bool, on::Bool) = ParticleInstance1D(ij, xy, ODEI, boundary, on)

ParticleInstance(ij::Tuple{Int,Int}, xy::Tuple{Float64,Float64}, ODEI::ODEIntegrator, boundary::Bool, on::Bool) = ParticleInstance2D(ij, xy, ODEI, boundary, on)

ParticleInstance(ij::Tuple{Int,Int,Int}, xy::Tuple{Float64,Float64}, ODEI::ODEIntegrator, boundary::Bool, on::Bool) = ParticleInstance2DLayer(ij, xy, ODEI, boundary, on)

# Debugging ParticleInstance
mutable struct MarkedParticleInstance <: AbstractMarkedParticleInstance
        Particle::AbstractParticleInstance
        time :: Float64
        state :: Vector{Any}
        errorReturnCode
end

Base.copy(s::ParticleInstance1D) = ParticleInstance1D(s.position_ij, s.position_xy, s.ODEIntegrator, s.boundary, s.on)
Base.copy(s::ParticleInstance2D) = ParticleInstance2D(s.position_ij, s.position_xy, s.ODEIntegrator, s.boundary, s.on)

# Regridding types:
"""Weights & Index (wni) FieldVector """
struct wni{TI<:SVector,TF<:SVector} <: FieldVector{4,SVector}
        xi::TI
        xw::TF
        yi::TI
        yw::TF
end


end