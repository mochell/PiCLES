module WaveGrowthModels2D

export WaveGrowth2D, init_particles!
export fields

using Architectures

using ModelingToolkit: Num, get_states

#using core_1D: MarkedParticleInstance
using ParticleMesh: OneDGrid, OneDGridNotes, TwoDGrid, TwoDGridNotes

using PiCLES.Operators.core_2D: ParticleDefaults as ParticleDefaults2D
#using core_2D: SeedParticle! as SeedParticle2D!

using SharedArrays
using Printf

import Oceananigans: fields
using Oceananigans.TimeSteppers: Clock


using PiCLES.Operators.mapping_2D
#includet("mapping_1D.jl")

function Base.show(io::IO, ow::AbstractModel)
    print(io, "WaveGrowth1D ", "\n",
        "├── grid: ", ow.grid, "\n",
        "├── layers: ", ow.layers, "\n",
        "├── clock: ", ow.clock,
        "├── State: ", size(ow.State), "\n",
        "├── ParticleCollection size: ", length(ow.ParticleCollection), "\n",
        "├── ODEs  \n",
        "|    ├── System: ", get_states(ow.ODEsystem), "\n",
        "|    ├── Defaults: ", ow.ODEdefaults, "\n",
        "|    └── Settings:  \n", ow.ODEsettings, "\n",
        "├── winds ", ow.winds, "\n",
        "├── currents ", ow.currents, "\n",
        "└── Perdiodic Boundary ", ow.periodic_boundary, "\n")
end


"""
    WaveGrowth1D{Grid, Lay, Tim, Clo, stat, PC, Ovar, Osys, Oses, Odev, bl_flag, bl, wnds, cur}

    WaveGrowth1D is a model for simulating wave growth in a 1D domain. 
    The model is based on the particle method, where each particle represents wave statistics at a given point in space.
    This structure contains all the information needed to run the model, including the grid, the particles, the ODE system, the boundary conditions, and the state of the model. 

"""
mutable struct WaveGrowth2D{Grid<:AbstractGrid,
                        Lay,
                        Tim,
                        Clo,
                        Int,
                        stat,
                        PC,
                        FPC,
                        Ovar,
                        Osys,
                        Oses,
                        Odev,
                        bl_flag,
                        bl,
                        wnds,
                        cur} <: Abstract2DModel
    grid::Grid
    layers::Lay      # number of layers used in the model, 1 is eneough
    timestepper::Tim      # silly Oceananigans
    clock::Clo
    dims::Int      # number of dimensions

    State::stat     # state of of the model at each grid point, for each layer it contains, energy, positions, group speed
    ParticleCollection::PC    # Collection (list) of Particles
    FailedCollection::FPC      # Collection (list) of Particles that failed to integrate

    ODEvars::Ovar     # list of variables in ODE system, have type Num from OrdinaryDiffEq and ModelingToolkit
    ODEsystem::Osys     # the ODE system used at each particle
    ODEsettings::Oses     # All setting needed to solve the ODEsystem
    ODEdefaults::Odev     # Dict{NUm, Float64} ODE defaults

    periodic_boundary::bl_flag # If true we use a period boundary 
    boundary::bl       # List of boundary points

    winds::wnds     # u, v, if needed u_x, u_y
    currents::cur      # u, v, currents

end




## 2D version
"""
mark_boundary(grid::TwoDGrid)
function that returns list of boundary nodes (tuples of indixes)
"""
function mark_boundary(grid::TwoDGrid)
    #get x and y coordinates
    xi = collect(range(1, stop=grid.Nx, step=1))
    yi = collect(range(1, stop=grid.Ny, step=1))

    # make boundary nodes
    a = [(yi[i], xi[j]) for i in 1:grid.Ny, j in [1, grid.Nx]]
    b = [(yi[i], xi[j]) for j in 1:grid.Nx, i in [1, grid.Nx]]
    #merge a and b
    return vec(vcat(a, b))
end


"""  
WaveGrowth2D(; grid, winds, ODEsys, ODEvars, layers, clock, ODEsets, ODEdefaults, currents, periodic_boundary, CBsets)
This is the constructor for the WaveGrowth1D model. The inputs are:
    grid             : the grid used in the model,
    winds            : the wind interpolation function here only 1D,
    ODEsys           : the ODE system used in the model,
    ODEvars          : the variables in the ODE system,
    layers           : the number of layers used in the model (default 1),
    clock            : the clock used in the model,
    ODEsets          : the ODE settings (type: ODESettings),
    ODEdefaults      : the ODE defaults (type: ParticleDefaults),
    currents         : the currents (not implimented yet),
    periodic_boundary: if true we use a periodic boundary (default true),
    CBsets           : the callback settings (not implimented yet).
"""
function WaveGrowth2D(; grid::TwoDGrid,
    winds, ODEsys, ODEvars,
    layers::Int=1,
    clock=Clock{eltype(grid)}(0, 0, 1),
    ODEsets::AbstractODESettings=nothing,  # ODE_settings
    ODEdefaults::ParticleDefaults2D=nothing,  # default_ODE_parameters
    currents=nothing,  # 
    periodic_boundary=true,
    CBsets=nothing)

    # initialize state {SharedArray} given grid and layers
    # Number of state variables 
    Nstate = 3
    if layers > 1
        State = SharedArray{Float64,4}(grid.Nx, grid.Ny, Nstate, layers)
    else
        State = SharedArray{Float64,3}(grid.Nx, grid.Ny, Nstate)
    end

    # initliaze boundary points (periodic_boundary)
    boundary = mark_boundary(grid)
    ODEdev = copy(ODEdefaults)

    ParticleCollection = []
    FailedCollection = Vector{AbstractMarkedParticleInstance}([])
    # particle initialization is  not done in the init_particle! method
    # for i in range(1,length = grid.Nx)
    #         SeedParticle!(ParticleCollection, State, i,
    #                       ODEsys, ODEdev, ODEsets,
    #                       gridnotes, winds, ODEsets.timestep, grid.Nx, boundary, periodic_boundary)
    # end        

    # check dimensions of all components:
    # check Nstate ODEvars

    # return WaveGrowth1D structure
    return WaveGrowth2D(
        grid,
        layers,
        nothing,
        clock, # ???
        2, # This is a 2D model
        State,
        ParticleCollection,
        FailedCollection, ODEvars,
        ODEsys,
        ODEsets,
        ODEdev, periodic_boundary,
        boundary, winds,
        currents)
end




"""
    fields(model::WaveGrowth)

Return a flattened `NamedTuple` of the State vector for a `WaveGrowth1D` model.
"""
fields(model::WaveGrowth2D) = (State=model.State,)
# # Oceananigans.Simulations interface
# fields(m::ContinuumIceModel) = merge(m.velocities, m.stresses)







# end of module
end