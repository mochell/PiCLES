#ENV["JULIA_INCREMENTAL_COMPILE"] = true

using Pkg
#Pkg.activate("../PiCLES/")

# using Pkg
#Pkg.activate(".")
Pkg.activate("PiCLES/")


import Plots as plt
using Setfield, IfElse

using PiCLES
using PiCLES.ParticleSystems: particle_waves_v5 as PW

import PiCLES: FetchRelations, ParticleTools
using PiCLES.Operators.core_2D: ParticleDefaults, InitParticleInstance, GetGroupVelocity
using PiCLES.Operators: TimeSteppers
using PiCLES.Simulations
using PiCLES.Operators.TimeSteppers: time_step!, movie_time_step!

using PiCLES.ParticleMesh: TwoDGrid, TwoDGridNotes, TwoDGridMesh
using PiCLES.Models.WaveGrowthModels2D

using Oceananigans.TimeSteppers: Clock, tick!
import Oceananigans: fields
using Oceananigans.Units
import Oceananigans.Utils: prettytime

using PiCLES.Architectures
using GLMakie

using PiCLES.Operators.core_2D: GetGroupVelocity, speed
using PiCLES.Plotting.movie: init_movie_2D_box_plot

using Revise

Revise.retry()
using PiCLES: ParticleInstance, ParticleInstance2D, ParticleInstance2DLayer




# debugging:
#using ProfileView

# %%
save_path = "plots/tests/T06_layers/"
mkpath(save_path)

# % Parameters
U10, V10 = 10.0, 10.0
dt_ODE_save = 30minutes
DT = 30minutes
# version 3

ODEpars, Const_ID, Const_Scg = PW.ODEParameters(r_g=0.85)


u_func(x, y, t) = U10 #* sin(t / (6 * 60 * 60 * 2π)) * sin(x / 50e3) * sin(y / 50e3)
v_func(x, y, t) = V10 #* cos(t / (6 * 60 * 60 * 2π)) * sin(x / 50e3) * sin(y / 50e3)

u(x, y, t) = u_func(x, y, t)
v(x, y, t) = v_func(x, y, t)
winds = (u=u, v=v)


grid = TwoDGrid(100e3, 11, 100e3, 11)
mesh = TwoDGridMesh(grid, skip=1);
gn = TwoDGridNotes(grid);

Revise.retry()


#ProfileView.@profview 
#ProfileView.@profview 
particle_system = PW.particle_equations(u, v, γ=Const_ID.γ, q=Const_ID.q, input=true, dissipation=true);

# define setting and standard initial conditions
WindSeamin = FetchRelations.get_minimal_windsea(U10, V10, DT);
#WindSeamin = FetchRelations.get_minimal_windsea(u(0, 0, 0), v(0, 0, 0), DT / 2)
#WindSeamin = FetchRelations.get_initial_windsea(u(0, 0, 0), v(0, 0, 0), DT/5)
lne_local = log(WindSeamin["E"])
cg_u_local = WindSeamin["cg_bar_x"]
cg_v_local = WindSeamin["cg_bar_y"]

ODE_settings = PW.ODESettings(
    Parameters=ODEpars,
    # define mininum energy threshold
    log_energy_minimum=lne_local,#log(FetchRelations.Eⱼ(0.1, DT)),
    #maximum energy threshold
    log_energy_maximum=log(27),#log(17),  # correcsponds to Hs about 16 m
    saving_step=dt_ODE_save,
    timestep=DT,
    total_time=T = 6days,
    adaptive=true,
    dt=1e-3, #60*10, 
    dtmin=1e-4, #60*5, 
    force_dtmin=true,
    callbacks=nothing,
    save_everystep=false)


# plt.heatmap(gn.x / 1e3, gn.y / 1e3, transpose(u.(mesh.x, mesh.y, 0)))
# plt.heatmap(gn.x / 1e3, gn.y / 1e3, transpose(v.(mesh.x, mesh.y, 0)))
# %% build model


Revise.retry()

wave_model = WaveGrowthModels2D.WaveGrowth2D(; grid=grid,
    winds=winds,
    ODEsys=particle_system,
    layers=2,
    ODEsets=ODE_settings,  # ODE_settings
    ODEinit_type="wind_sea",  # default_ODE_parameters
    periodic_boundary=false,
    boundary_type="same",
    #minimal_particle=FetchRelations.MinimalParticle(U10, V10, DT),
    movie=true)


size(wave_model.State)
typeof(wave_model.State) 



#particle_system
defaults = ParticleDefaults(0.1, 0.1, 0.1, 0.1, 0.1)

P1 =InitParticleInstance(particle_system, defaults, ODE_settings, (1, 2), true, true)


# 1D
ParticleInstance(1, 2.1, P1.ODEIntegrator, true, true)
# 3d
ParticleInstance2D((1, 2), (2.1, 3.1), P1.ODEIntegrator, true, true)
ParticleInstance((1, 2), (2.1, 3.1), P1.ODEIntegrator, true, true)
# 2d
ParticleInstance2DLayer((1, 2, 1), (2.1, 3.1), P1.ODEIntegrator, true, true)
ParticleInstance((1, 2, 3), (2.1, 3.1), P1.ODEIntegrator, true, true)


### build Simulation
#wave_simulation = Simulation(wave_model, Δt=10minutes, stop_time=4hours)#1hours)

wave_model.layers
typeof(wave_model.ParticleCollection)
isa(wave_model.ParticleCollection, AbstractArray)
isa(wave_model.ParticleCollection, Vector)

wave_simulation = Simulation(wave_model, Δt=10minutes, stop_time=6hour)#1hours)
initialize_simulation!(wave_simulation)


# %%
using ..Operators.mapping_2D

FailedCollection = Vector{AbstractMarkedParticleInstance}([])
model = wave_model
Δt = 10minutes
#@threads 
for a_particle in wave_model.ParticleCollection
    @info a_particle.position_ij
    mapping_2D.advance!(a_particle, model.State, FailedCollection,
        model.grid, model.winds, Δt,
        model.ODEsettings.log_energy_maximum,
        model.ODEsettings.wind_min_squared,
        model.periodic_boundary,
        model.ODEdefaults)
end



## test ParticleCollection MArray

using ..Operators.core_2D: SeedParticle!, SeedParticle

model = wave_model

print(defaults)
ParticleCollection=[]

# test if ParticleCollection is an array

#for i in range(1, length=gn.Nx)
        # SeedParticle!(ParticleCollection, model.State, i,
        #                 model.ODEsystem, defaults , model.ODEsettings,
        #                 gn, model.winds, model.ODEsettings.timestep,
        #                 model.boundary, model.periodic_boundary  )
#end

# %%

Revise.retry()

using StaticArrays
# tests
#MArray{Tuple{grid.Nx,grid.Ny,model.layers},ParticleInstance2DLayer}(undef);
#MArray{Tuple{grid.Nx,grid.Ny,model.layers},Union{ParticleInstance2DLayer,Int}}(undef);


SizedMatrix{Tuple{grid.Nx,grid.Ny,model.layers},ParticleInstance2DLayer}(undef);


#ParticleCollection2 = MArray{Tuple{grid.Nx,grid.Ny,model.layers},ParticleInstance2DLayer}(undef);
ParticleCollection3 = Array{Union{Nothing,ParticleInstance2DLayer}}(nothing, grid.Nx, grid.Ny, model.layers)
ParticleCollection = []


typeof(ParticleCollection)
typeof(ParticleCollection3)


# SeedParticle!(ParticleCollection3, model.State, (20, 2, 2),
#     model.ODEsystem, defaults, model.ODEsettings,
#     gn, model.winds, model.ODEsettings.timestep,
#     model.boundary, model.periodic_boundary)


for i in 1:model.grid.Nx, j in 1:model.grid.Ny, k in 1:model.layers

    ParticleCollection3[i,j,k]  =SeedParticle(model.State, (i,j,k),
        model.ODEsystem, defaults, model.ODEsettings,
        gn, model.winds, model.ODEsettings.timestep,
        model.boundary, model.periodic_boundary)
end

Pi = SeedParticle(model.State, (20, 2, 2),
    model.ODEsystem, defaults, model.ODEsettings,
    gn, model.winds, model.ODEsettings.timestep,
    model.boundary, model.periodic_boundary)



SeedParticle_mapper_to_array(f, s, b1, b2, b3, c1, c2, c3, d1, d2) = x -> f(s, x, b1, b2, b3, c1, c2, c3, d1, d2)

SeedParticle_i = SeedParticle_mapper_to_array(SeedParticle2D!,
    ParticleCollection, model.State,
    model.ODEsystem, defaults, model.ODEsettings,
    gridnotes, model.winds, model.ODEsettings.timestep,
    model.boundary, model.periodic_boundary)


model.State[20, 2, 2, :]




size(ParticleCollection2)
ParticleCollection3[10, 10, 1] = ParticleCollection[1] #ParticleCollection[1]


# %%




#init_state_store!(wave_simulation, save_path)
#wave_simulation.model.MovieState = wave_simulation.model.State

@time run!(wave_simulation, cash_store=true, debug=true)
#reset_simulation!(wave_simulation)
# run simulation
#ProfileView.@profview run!(wave_simulation, cash_store=true, debug=true)


istate = wave_simulation.store.store[end];
p1 = plt.heatmap(gn.x / 1e3, gn.y / 1e3, istate[:, :, 1])

