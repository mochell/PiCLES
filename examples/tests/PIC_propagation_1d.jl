using Plots
using BenchmarkTools
using Revise

using IfElse
using Statistics

using ModelingToolkit: Num, @unpack, @register_symbolic, Symbolics, @named, ODESystem

using HDF5
using JLD2

using DocStringExtensions

# Particle Model modules
push!(LOAD_PATH, joinpath(pwd(), "code/"))
using ParticleMesh: OneDGrid, OneDGridNotes

#using core_1D: periodic_BD_single_PI!, show_pos!, periodic_condition_x
push!(LOAD_PATH, joinpath(pwd(), "code/Core"))
using core_1D: ParticleDefaults
using TimeSteppers

import ParticleInCell
import FetchRelations
using WaveGrowthModels1D

includet("../../ParticleMesh.jl")
using Debugging
using Printf

import particle_waves_v4#: particle_equations, ODESettings
PW = particle_waves_v4

using Oceananigans.TimeSteppers: Clock, tick!
import Oceananigans: fields
using Oceananigans.Units

# for Callbacks
push!(LOAD_PATH, joinpath(pwd(), "code/Simulations"))
using Simulations

using WindEmulator

"""
This load a parameter file, executes a 1D run, saves the data and run statistics.

"""

using InputOutput: Argsettings, parse_args

# %%
# Default values
save_path_base = "data/1D_gaussian/"
plot_path_base = "plots/static/"

parset = "1D_varying/"


### Boundary Conditions
periodic_boundary = false

# parametric wind forcing
U10, V = 10, 5 #m/s

r_g0 = 0.85
c_β = 4e-2
C_e0 = (2.35 / r_g0) * 2e-3 * c_β
γ = 0.7

#  rescale parameters for the right units.
T = 24 * 2 * 60 * 60 # seconds
Lx = 300 * 10e3  # km
DT = Float64(20 * 60) # seconds
Nx = 40
dt_ODE_save = 10 # 3 min

grid1d = OneDGrid(1e3, Lx - 1e3, Nx)

# create ID and save name
ID = "empl_eq_v4"
#save_path = joinpath( "data/1D_static/", parsed_args["ID"] )
save_path = save_path_base * parset * "/" * ID * "/"
plot_path = plot_path_base * parset * "/"

mkpath(save_path)
# %%
@printf "Init Forcing Field\n"
# create wind test fucntion

@register_symbolic u(x, t)

x_scale = 600e3
t_scale = (60 * 60 * 24 * 1.5)
dx = 5e3

fake_winds(x, t) = U10 + x * 0 + t * 0# 

wind_grid = IdealizedWindGrid(fake_winds, (Lx=Lx, T=T), (dx=dx, dt=DT))

# % define wind forcing globally as interp functions
interp_winds = wind_interpolator(wind_grid);

#convert to registered_symbolic need foir ODESystem
u(x, t) = interp_winds.u(x, t)

# -------------- start model definition -------------------------
# %% Load Particle equations and derive ODE system

############ try v3 again first and check ways to import module
#particle_equations = PW.particle_equations(u, γ=γ,input=true, dissipation=true, peak_shift=true )
particle_equations = PW.particle_rays()
@named particle_system = ODESystem(particle_equations)


# define variables based on particle equation
t, x, c̄_x, lne, r_g, C_α, g, C_e = PW.init_vars_1D()


# %% define storing stucture and populate inital conditions
default_ODE_parameters = Dict(
    r_g => 1 / r_g0,
    C_α => -1.41,
    g => 9.81,
    C_e => C_e0,
)


ODE_settings = PW.ODESettings(
    Parameters=default_ODE_parameters,
    # define mininum energy threshold
    log_energy_minimum=log(FetchRelations.Eⱼ(0.1, DT)),
    #maximum energy threshold
    log_energy_maximum=log(17),  # correcsponds to Hs about 16 m
    saving_step=dt_ODE_save,
    timestep=DT,
    total_time=T,
    adaptive=true,
    dt=1e-3, #60*10, 
    dtmin=1e-12, #60*5, 
    force_dtmin=true,
)


# Default values for particle
particle_defaults = ParticleDefaults(log(FetchRelations.Eⱼ(5.0, DT)), 1e-4, 0.0)

vars = PW.init_vars_1D()
wave_model = WaveGrowthModels1D.WaveGrowth1D(; grid=grid1d,
    winds=u,
    ODEsys=particle_system,
    ODEvars=vars,
    layers=1,
    ODEsets=ODE_settings,  # ODE_settings
    ODEdefaults=particle_defaults,  # default_ODE_parameters
    periodic_boundary=false,
    boundary_type="default"
)


#using storing: init_state_store!, push_state_to_storage!
#using Simulations: Simulation

Revise.retry()
# %% run simulation
wave_simulation = Simulation(wave_model, Δt=20minutes, stop_time=2hours)
initialize_simulation!(wave_simulation, particle_initials=wave_model.ODEdefaults)
# or init_particles!( wave_model, defaults= wave_model.ODEdefaults, verbose = wave_simulation.verbose )

function set_u_and_t!(integrator, u_new, t_new)
    integrator.u = u_new
    integrator.t = t_new
end


hs = 1 # meter
energy = (hs / 4)^2
c_g = 5 # m/s

for i in range(10, Integer(floor(length(wave_model.ParticleCollection) * 2 / 3)), step=1)

    PII = wave_model.ParticleCollection[i].ODEIntegrator
    ui = [log(energy), c_g, PII.u[3]]
    set_u_and_t!(PII, ui, PII.t)
    u_modified!(PII, true)
    @info i, wave_model.ParticleCollection[i].ODEIntegrator.u
end

wave_simulation.model.State[5:15, 1]
wave_simulation.model.State[5:15, 1] / wave_simulation.model.State[5:15, 2] / 2



#### 
# add flag if initial conditions are written to state vector and storage
# make rountine that write particle initial conditions to state vector and storage
# reset State vetor each time step!!
####


pp1 = [log(energy), c_g, 0.0]
state1 = core_1D.GetParticleEnergyMomentum(pp1)
state1
# using core_1D
# for i in 1:1e4
#     state1 = core_1D.GetParticleEnergyMomentum(pp1)
#     pp2 = core_1D.GetVariablesAtVertex(state1, 0.0)
#     @info pp1 -  pp2
#     pp1 = pp2
# end
# using DifferentialEquations
# step!(wave_model.ParticleCollection[1].ODEIntegrator, DT, true)
# wave_model.ParticleCollection[1].ODEIntegrator.u

# wave_model.ParticleCollection[10].ODEIntegrator.u
# step!(wave_model.ParticleCollection[10].ODEIntegrator, DT*10, true)
# wave_model.ParticleCollection[10].ODEIntegrator.u

#plot(wave_model.State[:, 2])

#using storing: init_state_store!
#init_state_store!(wave_simulation, save_path)

run!(wave_simulation, store=false, cash_store=true, debug=true)

#reset_simulation!(wave_simulation)
# wave_simulation.stop_time = Inf
# run!(wave_simulation, store=false)

@info "... finished\n"

# %% open wave_simulation.store and plot variables
## when using cash store
get_1d_data(store::Vector{Any}) = permutedims(cat(store..., dims=3), (3, 1, 2))
store_waves_data = get_1d_data(wave_simulation.store.store);
wave_x = OneDGridNotes(wave_model.grid).x;
wave_time = collect(range(wave_simulation.Δt, wave_simulation.stop_time + wave_simulation.Δt, step=wave_simulation.Δt));

# when using state store
# convert HDF5 data to Array
# store_waves_data = Array(wave_simulation.store.store["data"])
# wave_x = wave_simulation.store.store["x"][:]
# wave_time = wave_simulation.store.store["time"][:]


store_waves_energy = store_waves_data[:, :, 1];
store_waves_mx = store_waves_data[:, :, 2];
store_waves_my = store_waves_data[:, :, 3];
#print(wave_simulation.store.store["var_names"][:])

cg = store_waves_energy ./ store_waves_mx ./ 2

store_waves_energy[2, 9:12]
#contourf(wave_x / dx, wave_time / DT, store_waves_my, levels=20, colormap=:blues)
#contourf(wave_x / dx, wave_time / DT, store_waves_mx, levels=20, colormap=:blues)

#heatmap(wave_x / dx, wave_time / DT, store_waves_energy, levels=20, colormap=:blues)

heatmap(wave_x / dx, wave_time / DT, cg, levels=20, colormap=:blues)
#heatmap(wave_x / dx, wave_time / DT, store_waves_energy, levels=20, colormap=:blues)
#heatmap(wave_x / dx, wave_time / DT, store_waves_energy, levels=20, colormap=:blues)
plot!(xlabel="x (dx)", ylabel="time (DT)", title="cg") |> display


# %% plot failed particles, if there any
using ParticleTools
if ~isempty(wave_simulation.model.FailedCollection)
    @printf "plot failed particles\n"
    ParticleTools.PlotFailedParticles(wave_simulation.model.FailedCollection, Array(wave_simulation.store.store["data"]), ID, DT, dx,
        savepath=false, Npar=4)

end


### ------------ up to here ----------------------

# %% Statistics
PI = wave_simulation.model.ParticleCollection[10]
ParticleInCell.compute_weights_and_index(wave_model.grid, PI.ODEIntegrator.u[1])


store_waves_data = wave_simulation.store.store["data"]

energy = store_waves_data[end-1, :, 1] # State_collect[end-1][:,1] #GP.sel(state= 0 ).T # energy
m_x = store_waves_data[end-1, :, 2] # State_collect[end-1][:,2] #GP.sel(state= 1 ).T # energy
cg = energy ./ m_x ./ 2 # c_g


# %% Checks
@printf "\n size of collected states %s" size(store_waves_data)
@printf "\n size of each state  %s" size(store_waves_data[1, :, :])
@printf "\n length of timerange  %s" size(wave_simulation.store.store["time"])

# %% Saving Fields
@printf " write attributes to store and close \n"

mkpath(save_path)
@printf "created model run folder %s \n" save_path

# for i in 1:length(State_collect)
#         store_waves_data[i,:,:] = State_collect[i]
# end


# %% Save Winds
forcing = (u=wind_grid.u, v=nothing)
coords = (x=wind_grid.x, time=wind_grid.t)

add_winds_forcing_to_store!(wave_simulation.store, forcing, coords)
show_stored_data(wave_simulation)
close_store!(wave_simulation)


# %% Save Particles
@printf "save particle \n"
save_object(joinpath(save_path, "particles.jld2"), wave_simulation.model.ParticleCollection)


# %%
