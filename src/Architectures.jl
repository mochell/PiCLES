module Architectures

using SharedArrays
using StaticArrays


export AbstractGrid, AbstractODESettings, AbstractParticleInstance, AbstractMarkedParticleInstance, Abstract1DModel, Abstract2DModel, AbstractModel, AbstractStore, AbstractParticleSystem, IDConstantsInstance, ScgConstantsInstance
export StateTypeL1, StateTypeLN, AbstractState

abstract type AbstractGrid end

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

# abstract type AbstractState{TS} end


#All posiible types of a single-layer StateVectors
StateTypeL1 = Union{SharedArray{Float64,3}, MArray}#  <: AbstractState{Nothing}
StateTypeLN = Union{SharedArray{Float64,4}} #<: AbstractState{Nothing}
StateType = Union{StateTypeL1, StateTypeLN}
Position_ij_Type = Union{Int, Tuple{Int,Int}, Tuple{Int,Int,Int}}

abstract type AbstractStore end

abstract type AbstractParticleSystem end

end