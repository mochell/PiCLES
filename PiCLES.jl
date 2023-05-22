"""
Main module for 'PiCLES.jl'
"""
module PiCLES

# external modules
#using ModelingToolkit: Num, @unpack, @register_symbolic, Symbolics, @named, ODESystem
using HDF5
using DocStringExtensions
#using Printf

export 
    # run simulation
    Simulations,
    
    #Operators/Core
    Operators,

    # models
    Models,

    # Particle Systems
    ParticleSystems,

    # grids
    ParticleMesh,

    # utils
    Utils, Debugging, FetchRelations, ParticleTools, WindEmulator

    #externals


include("Architectures.jl")

include("ParticleSystems/ParticleSystems.jl")
using .ParticleSystems

include("FetchRelations.jl")
using .FetchRelations

include("Operators/Operators.jl")
include("Simulations/Simulations.jl")

#include("Grids/Grids.jl")
include("ParticleMesh.jl")
include("ParticleInCell.jl")
#include("Particles/Particles.jl")

include("Utils/Utils.jl")
using .Utils


using .Simulations
using .ParticleMesh
using .ParticleInCell

include("Models/Models.jl")
using .Models

include("visualization/Plotting.jl")
using .Plotting

end