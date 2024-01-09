using Pkg

Pkg.activate("PiCLES/")

using Base.Threads
@info "Num. of threads", nthreads()

# import Plots as plt
using Setfield
using PiCLES.ParticleSystems: particle_waves_v5 as PW

import PiCLES: FetchRelations, ParticleTools
using PiCLES.Operators.core_2D: ParticleDefaults, InitParticleInstance, GetGroupVelocity
using PiCLES.Simulations
using PiCLES.Operators.TimeSteppers: time_step!, movie_time_step!

using PiCLES.ParticleMesh: TwoDGrid, TwoDGridNotes, TwoDGridMesh
using PiCLES.Models.WaveGrowthModels2D

using Oceananigans.TimeSteppers: Clock, tick!
import Oceananigans: fields
using Oceananigans.Units
import Oceananigans.Utils: prettytime

using PiCLES.Architectures
#using GLMakie
#using PiCLES.Plotting.movie: init_movie_2D_box_plot

#using ProfileView
using BenchmarkTools
#using Revise

# debugging:
#using ProfileView

@info "precompiled!"
# %%
save_path = "plots/examples/homogenous_box/"
mkpath(save_path)

# % Parameters
U10, V10 = 10.0, 10.0
dt_ODE_save = 30minutes
DT = 30minutes
# version 3
r_g0 = 0.85

# function to define constants 
Const_ID = PW.get_I_D_constant()
@set Const_ID.γ = 0.88
Const_Scg = PW.get_Scg_constants(C_alpha=-1.41, C_varphi=1.81e-5)


u(x, y, t) = U10 #* sin(t / (6 * 60 * 60 * 2π)) * sin(x / 50e3) * sin(y / 50e3)
v(x, y, t) = V10 #* cos(t / (6 * 60 * 60 * 2π)) * sin(x / 50e3) * sin(y / 50e3)

winds = (u=u, v=v)

grid = TwoDGrid(100e3, 51, 100e3, 51)
#grid = TwoDGrid(20e3, 21, 20e3, 21)
mesh = TwoDGridMesh(grid, skip=1)
gn = TwoDGridNotes(grid)
#Revise.retry()


# define variables based on particle equation
particle_system = PW.particle_equations(u, v, γ=0.88, q=Const_ID.q, input=true, dissipation=true)
#particle_equations = PW3.particle_equations_vec5(u, v, u, v, γ=0.88, q=Const_ID.q);

# define V4 parameters absed on Const NamedTuple:
default_ODE_parameters = (r_g=r_g0, C_α=Const_Scg.C_alpha,
    C_φ=Const_ID.c_β, C_e=Const_ID.C_e, g=9.81)

# define setting and standard initial conditions
WindSeamin = FetchRelations.get_minimal_windsea(U10, V10, DT)
#WindSeamin = FetchRelations.get_minimal_windsea(u(0, 0, 0), v(0, 0, 0), DT / 2)
#WindSeamin = FetchRelations.get_initial_windsea(u(0, 0, 0), v(0, 0, 0), DT/5)
lne_local = log(WindSeamin["E"])
cg_u_local = WindSeamin["cg_bar_x"]
cg_v_local = WindSeamin["cg_bar_y"]

ODE_settings = PW.ODESettings(
    Parameters=default_ODE_parameters,
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


default_particle = ParticleDefaults(lne_local, cg_u_local, cg_v_local, 0.0, 0.0)

#Revise.retry()
wave_model = WaveGrowthModels2D.WaveGrowth2D(; grid=grid,
    winds=winds,
    ODEsys=particle_system,
    ODEsets=ODE_settings,  # ODE_settings
    ODEinit_type="wind_sea",  # default_ODE_parameters
    periodic_boundary=false,
    boundary_type="same",
    #minimal_particle=FetchRelations.MinimalParticle(U10, V10, DT),
    movie=true);


#global wave_model.ParticleCollection
### build Simulation
#wave_simulation = Simulation(wave_model, Δt=10minutes, stop_time=4hours)#1hours)
wave_simulation = Simulation(wave_model, Δt=10minutes, stop_time=6hour);#1hours)
initialize_simulation!(wave_simulation)

# # %% single time steps
# import PiCLES.Operators.TimeSteppers: time_step!
# for i in 1:3
#     wave_simulation.model.State .= 0.0
#     @time time_step!(wave_simulation.model, wave_simulation.Δt, debug=false)
#     #time_step!(wave_simulation.model, wave_simulation.Δt, debug=false)
#     istate = wave_simulation.model.State
#     p1 = plt.heatmap(gn.x / 1e3, gn.y / 1e3, istate[:, :, 1])
#     display(p1)
#     sleep(0.5)
# end
# # %%


run!(wave_simulation, cash_store=false, debug=false);
reset_simulation!(wave_simulation)
wave_simulation.stop_time = 6hour
@time run!(wave_simulation, cash_store=false, debug=false);

# istate = wave_simulation.model.State;
# p1 = plt.heatmap(gn.x / 1e3, gn.y / 1e3, istate[:, :, 1])

#wave_simulation.model.ParticleCollection[10]
# run simulation
#ProfileView.@profview run!(wave_simulation, cash_store=true, debug=true)

#istate = wave_simulation.store.store[end];
#p1 = plt.heatmap(gn.x / 1e3, gn.y / 1e3, istate[:, :, 1])`
