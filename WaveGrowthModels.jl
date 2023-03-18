module WaveGrowthModels

export WaveGrowth1D, time_step!,init_particles!
export fields

using Oceananigans.Models: AbstractModel
using Oceananigans.TimeSteppers: Clock, tick!

using ModelingToolkit: Num, get_states


using particle_waves_v3: ODESettings, init_vars_1D
using core_1D: ParticleDefaults, SeedParticle!,MarkedParticleInstance
using ParticleMesh: AbstractGrid, OneDGrid, OneDGridNotes


using SharedArrays
using Printf

import Oceananigans: fields

using mapping_1D
#includet("mapping_1D.jl")

function Base.show(io::IO, ow::AbstractModel)
  print(io, "WaveGrowth1D ", "\n",
                "├── grid: ", ow.grid, "\n",
                "├── layers: ", ow.layers, "\n",
                "├── clock: ", ow.clock,
                "├── State: ", size(ow.State), "\n",
                "├── ParticleCollection size: ", length(ow.ParticleCollection), "\n",
                "├── ODEs  \n",
                "|    ├── System: ", get_states(ow.ODEsystem) , "\n",
                "|    ├── Defaults: ", ow.ODEdefaults , "\n",
                "|    └── Settings:  \n", ow.ODEsettings , "\n",
                "├── winds ", ow.winds , "\n",
                "├── currents ", ow.currents , "\n",
                "└── Perdiodic Boundary ", ow.periodic_boundary , "\n")
end


"""
    WaveGrowth1D{Grid, Lay, Tim, Clo, stat, PC, Ovar, Osys, Oses, Odev, bl_flag, bl, wnds, cur}

    WaveGrowth1D is a model for simulating wave growth in a 1D domain. 
    The model is based on the particle method, where each particle represents wave statistics at a given point in space.
    This structure contains all the information needed to run the model, including the grid, the particles, the ODE system, the boundary conditions, and the state of the model. 

"""
mutable struct WaveGrowth1D{Grid <: OneDGrid,
                            Lay,
                            Tim,
                            Clo,
                            stat,
                            PC,
                            Ovar,
                            Osys,
                            Oses,
                            Odev,
                            bl_flag,
                            bl,
                            wnds,
                            cur} <: AbstractModel{Nothing}
        grid    :: Grid
        layers  :: Lay      # number of layers used in the model, 1 is eneough
    timestepper :: Tim      # silly Oceananigans
          clock :: Clo

          State :: stat     # state of of the model at each grid point, for each layer it contains, energy, positions, group speed
ParticleCollection :: PC    # Collection (list) of Particles

        ODEvars :: Ovar     # list of variables in ODE system, have type Num from OrdinaryDiffEq and ModelingToolkit
      ODEsystem :: Osys     # the ODE system used at each particle
    ODEsettings :: Oses     # All setting needed to solve the ODEsystem
    ODEdefaults :: Odev     # Dict{NUm, Float64} ODE defaults

periodic_boundary :: bl_flag # If true we use a period boundary 
       boundary :: bl       # List of boundary points

          winds :: wnds     # u, v, if needed u_x, u_y
       currents :: cur      # u, v, currents

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
      CBsets           : the callback settings (not implimented yet).
  """
function WaveGrowth1D(; grid, winds, ODEsys, ODEvars,
                            layers :: Int                     = 1,
                            clock                             = Clock{eltype(grid)}(0, 0, 1),
                            ODEsets  :: ODESettings           = nothing,  # ODE_settings
                            ODEdefaults :: ParticleDefaults   = nothing,  # default_ODE_parameters
                            currents                          = nothing,  # 
                            periodic_boundary                 = true,
                            CBsets                            = nothing )

        # initialize state {SharedArray} given grid and layers
        # Number of state variables 
        Nstate    = 3
        if layers > 1
          State     = SharedArray{Float64, 3}(grid.Nx, Nstate, layers)
        else
          State     = SharedMatrix{Float64}(grid.Nx, Nstate)
        end

        # initliaze boundary points (periodic_boundary)
        if ~periodic_boundary  # if false, define boundary points here:
          boundary = [1 , grid.Nx]
        end

        # define variables based on particle equationß
        #pass through from ODEvars
        #t, x, c̄_x, lne, r_g, C_α, g, C_e = particle_waves_v3.init_vars_1D()

        gridnotes = OneDGridNotes(grid)
 
        ODEdev = copy(ODEdefaults)
        ### seed particles
        # @printf "seed PiCLES ... \n"
        # #seed particles
        ParticleCollection=[]

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
                            State,
                            ParticleCollection,
                            ODEvars,
                            ODEsys,
                            ODEsets,
                            ODEdev, 
                            periodic_boundary,
                            boundary,
                            winds,
                            currents)
end



"""
    fields(model::WaveGrowth1D)

Return a flattened `NamedTuple` of the State vector for a `WaveGrowth1D` model.
"""
fields(model::WaveGrowth1D) = (State = model.State, )
# # Oceananigans.Simulations interface
# fields(m::ContinuumIceModel) = merge(m.velocities, m.stresses)



"""
time_step!(model, Δt; callbacks=nothing)

advances model by 1 time step:
1st) the model.ParticleCollection is advanced and then 
2nd) the model.State is updated.
clock is ticked by Δt

callbacks are not implimented yet

"""
function time_step!(model::WaveGrowth1D, Δt; callbacks=nothing)

  FailedCollection = Vector{MarkedParticleInstance}([])

  for a_particle in model.ParticleCollection
          #@show a_particle.position_ij
          mapping_1D.advance!(    a_particle, model.State, FailedCollection, 
                                  model.grid, model.winds , Δt , 
                                  model.ODEsettings , model.periodic_boundary )
  end

  #@printf "re-mesh"
  for a_particle in model.ParticleCollection
          mapping_1D.remesh!(a_particle, model.State, model.winds, model.clock.time, model.ODEsettings,  Δt)
  end

  tick!(model.clock, Δt)
end


"""
SeedParticle_mapper(f, p, s, b1, b2, b3, c1, c2, c3, c4, d1, d2 ) = x -> f( p, s, x, b1, b2, b3, c1, c2, c3, c4, d1, d2 )
maps to SeedParticle! function
"""
SeedParticle_mapper(f, p, s, b1, b2, b3, c1, c2, c3, c4, d1, d2 )  = x -> f( p, s, x, b1, b2, b3, c1, c2, c3, c4, d1, d2 )


"""
init_particle!(model ; defaults::Dict{Num, Float64} = nothing, verbose::Bool=false )

initialize the model.ParticleCollection based on the model.grid and the defaults. 
If defaults is nothing, then the model.ODEdev is used.
usually the initilization uses wind constitions to seed the particles.
"""
function init_particles!(model ; defaults::Dict{Num, Float64} = copy(ParticleDefaults(0.0, 1e-2, log( 4e-8 ))) , verbose::Bool=false )
         
        #defaults        = isnothing(defaults) ? model.ODEdev : defaults
        gridnotes       = OneDGridNotes(model.grid)
        if verbose
                @info "seed PiCLES ... \n"
        end

        ParticleCollection=[]
        SeedParticle_i  = SeedParticle_mapper(SeedParticle!,  ParticleCollection, model.State,
                                model.ODEsystem, defaults , model.ODEsettings,
                                gridnotes, model.winds, model.ODEsettings.timestep, model.grid.Nx, 
                                model.boundary, model.periodic_boundary  )

        map( SeedParticle_i , range(1,length = model.grid.Nx) )

        # print(defaults)
        # ParticleCollection=[]
        # for i in range(1,length = grid1d.Nx)
        #         SeedParticle!(ParticleCollection, model.State, i,
        #                         model.ODEsystem, defaults , model.ODEsettings,
        #                         gridnotes, model.winds, model.ODEsettings.timestep, model.grid.Nx,
        #                         model.boundary, model.periodic_boundary  )
        # end

        model.ParticleCollection = ParticleCollection
        nothing
end



# end of module
end