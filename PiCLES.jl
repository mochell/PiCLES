"""
Main module for 'PiCLES.jl'
"""
module PiCLES   

export 
    # run simulation
    Simulations,
    
    # models
    WaveGrowthModels,

    # Particle Systems
    #Particles,

    # grids

    # utils
    #Debugging,FetchRelations,ParticleTools


using .Simulations
#using .Models

include("Simulations/Simulations.jl")
#include("Models/Models.jl")
#include("Grids/Grids.jl")
#include("Particles/Particles.jl")
end