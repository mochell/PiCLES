module Architectures

export AbstractGrid, AbstractODESettings, AbstractParticleInstance, AbstractMarkedParticleInstance, Abstract1DModel, Abstract2DModel, AbstractModel, AbstractStore, AbstractParticleSystem
abstract type AbstractGrid end

abstract type AbstractODESettings end
abstract type AbstractParticleInstance end
abstract type AbstractMarkedParticleInstance end

"""
    AbstractModel
Abstract supertype for models.
"""
abstract type AbstractModel{TS} end

abstract type Abstract1DModel <: AbstractModel{Nothing} end
abstract type Abstract2DModel <: AbstractModel{Nothing} end

abstract type AbstractStore end

abstract type AbstractParticleSystem end

end