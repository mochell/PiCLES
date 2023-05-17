module WaveGrowthModels1D

export WaveGrowth1D
export fields

using Architectures
using ModelingToolkit: get_states


#using core_1D: MarkedParticleInstance
using ParticleMesh: OneDGrid, OneDGridNotes, TwoDGrid, TwoDGridNotes

using PiCLES.Operators.core_1D: ParticleDefaults as ParticleDefaults1D

using SharedArrays
using Printf

import Oceananigans: fields
using Oceananigans.TimeSteppers: Clock, tick!


using PiCLES.Operators.mapping_1D
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
mutable struct WaveGrowth1D{Grid<:AbstractGrid,
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
  bl_type,
  wnds,
  cur} <: Abstract1DModel
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
  boundary::bl        # List of boundary points
  boundary_defaults::bl_type # Dict{NUm, Float64} ODE defaults

  winds::wnds     # u, v, if needed u_x, u_y
  currents::cur      # u, v, currents

end

"""  
  WaveGrowth1D(; grid, winds, ODEsys, ODEvars, layers, clock, ODEsets, ODEdefaults, currents, periodic_boundary, CBsets)
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
      boundary_type    : the type of boundary, either "wind_sea" or "flat" (default "wind_sea"),
      CBsets           : the callback settings (not implimented yet).
"""
function WaveGrowth1D(; grid::OneDGrid, winds, ODEsys, ODEvars,
  layers::Int=1,
  clock=Clock{eltype(grid)}(0, 0, 1),
  ODEsets::AbstractODESettings=nothing,  # ODE_settings
  ODEdefaults::ParticleDefaults1D=nothing,  # default_ODE_parameters
  currents=nothing,  # 
  periodic_boundary=true,
  boundary_type="wind_sea", # or "flat"
  CBsets=nothing)

  # initialize state {SharedArray} given grid and layers
  # Number of state variables 
  Nstate = 3
  if layers > 1
    State = SharedArray{Float64,3}(grid.Nx, Nstate, layers)
  else
    State = SharedMatrix{Float64}(grid.Nx, Nstate)
  end

  # initliaze boundary points (periodic_boundary)
  if ~periodic_boundary  # if false, define boundary points here:
    boundary = [1, grid.Nx]
  else
    boundary = []
  end

  if boundary_type == "wind_sea"
    boundary_defaults = nothing
    @info "use wind_sea boundary"
  elseif boundary_type == "zero"
    @info "use zero boundary"
    boundary_defaults = copy(ParticleDefaults1D(-11.0, 1e-3, 0.0))
  elseif boundary_type == "default"
    @info "use default value boundary"
    boundary_defaults = copy(ODEdefaults)
  else
    error("boundary_type must be either 'wind_sea','zero', or 'default' ")
  end


  ODEdev = copy(ODEdefaults)
  ### seed particles
  # @printf "seed PiCLES ... \n"
  # #seed particles
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
  return WaveGrowth1D(grid, layers,
    nothing,
    clock, # ???
    1, # This is a 1D model
    State,
    ParticleCollection,
    FailedCollection,
    ODEvars,
    ODEsys,
    ODEsets,
    ODEdev,
    periodic_boundary,
    boundary,
    boundary_defaults,
    winds,
    currents)
end



"""
    fields(model::WaveGrowth1D)

Return a flattened `NamedTuple` of the State vector for a `WaveGrowth1D` model.
"""
fields(model::WaveGrowth1D) = (State=model.State,)
# # Oceananigans.Simulations interface
# fields(m::ContinuumIceModel) = merge(m.velocities, m.stresses)





# end of module
end