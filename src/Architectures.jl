module Architectures

using SharedArrays
using StaticArrays


export AbstractGrid, AbstractODESettings, AbstractParticleInstance, AbstractMarkedParticleInstance, Abstract1DModel, Abstract2DModel, AbstractModel, AbstractStore, AbstractParticleSystem, StateTypeL1, IDConstantsInstance, ScgConstantsInstance

abstract type AbstractGrid end

abstract type AbstractGridStatistics <: AbstractGrid end

abstract type CartesianGrid <: AbstractGrid end
abstract type CartesianGrid1D <: CartesianGrid end
abstract type CartesianGrid2D <: CartesianGrid end

abstract type TripolarGrid <: AbstractGrid end
# abstract type MOM6_2_3 <: TripolarGrid end


abstract type AbstractODESettings end
abstract type AbstractParticleInstance end
abstract type AbstractMarkedParticleInstance end

abstract type IDConstantsInstance end
abstract type ScgConstantsInstance end


"""
    AbstractModel
Abstract supertype for models.
"""
abstract type AbstractModel{TS} end

abstract type Abstract1DModel <: AbstractModel{Nothing} end
abstract type Abstract2DModel <: AbstractModel{Nothing} end

#All posiible types of a single-layer StateVectors
StateTypeL1 = Union{SharedArray{Float64,3},MArray}

abstract type AbstractStore end

abstract type AbstractParticleSystem end

end