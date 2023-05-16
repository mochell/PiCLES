"""
Main module for 'PiCLES.jl'
"""
module PiCLES   

export 
    # run simulation
    Simulations,
    
    # models
    WaveGrowthModels

    # Particle Systems
    #Particles,

    # grids

    # utils
    #Debugging,FetchRelations,ParticleTools



include("Simulations/Simulations.jl")
#include("Models/Models.jl")
#include("Grids/Grids.jl")
#include("Particles/Particles.jl")



using .Simulations
#using .Models
end