# using Pkg; 
# Pkg.resolve()
# Pkg.instantiate()

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

import ParticleInCell
import FetchRelations

push!(LOAD_PATH, joinpath(pwd(), "code/Core"))
using WaveGrowthModels1D

includet("../ParticleMesh.jl")
using Debugging

import particle_waves_v3beta#: particle_equations, ODESettings
PW = particle_waves_v3beta


using Printf

using Oceananigans.TimeSteppers: Clock, tick!
import Oceananigans: fields
using Oceananigans.Units

# for Callbacks
#using core_1D: periodic_BD_single_PI!, show_pos!, periodic_condition_x
push!(LOAD_PATH, joinpath(pwd(), "code/Core"))
using core_1D: ParticleDefaults

push!(LOAD_PATH, joinpath(pwd(), "code/Simulations"))
using Simulations

using WindEmulator


"""
This script load a parameter file, executes a 1D run, saves the data and run statistics.

"""
#using InputOutput: Argsettings, parse_args

# %%
# Default values
save_path_base = "data/1D_gaussian/"
plot_path_base = "plots/static/"

parset = "1D_varying/"


### Boundary Conditions
periodic_boundary = false

# parametric wind forcing
U10, V = 5, 2 #m/s

r_g0 = 0.85
c_β = 4e-2
C_e0 = (2.35 / r_g0) * 2e-3 * c_β
γ = 0.7

#  rescale parameters for the right units.
T = 1.3 * 24  * 60 * 60 # seconds
Lx = 300 * 1e3  # km
DT = Float64(20 * 60) # seconds
Nx = 200
dt_ODE_save = 10 # 3 min

grid1d = OneDGrid(1e3, Lx - 1e3, Nx)
#grid1d.dx

# create ID and save name
ID = "exmpl_1D_v3"
#save_path        = joinpath( "data/1D_static/", parsed_args["ID"] )
save_path = save_path_base * parset * "/" * ID * "/"
plot_path = plot_path_base * parset * "/"


mkpath(save_path)

# %%
@printf "Init Forcing Field\n"
# create wind test fucntion

@register_symbolic u(x, t)
@register_symbolic u_x(x, t)

x_scale = 60e3
t_scale = (60 * 60 * 24 * 0.8)
dx = 20e3

Revise.retry()
fake_winds(x, t) =  U10 + x * 0 + t* 0#  
#fake_winds(x, t) = slopped_blob(x, t, U10, V, T, x_scale, t_scale; x0=3e3)# <: WindEmulationFunction2

wind_ = IdealizedWindGrid(fake_winds, (Lx=Lx, T=T), (dx=dx, dt=DT))

# % define wind forcing globally as interp functions
interp_winds = wind_interpolator(wind_)

#convert to registered_symbolic need foir ODESystem
u(x, t) = interp_winds.u(x, t)

# add wind gradient
u_x(x, t) = Interpolations.gradient(interp_winds.u, x)[1]
#u_x(x, t) = Interpolations.gradient(u_, x )[1]

#contourf(wind_.x /1e3, wind_.t / 60 /60, transpose(wind_.u))
contourf(wind_.x / dx, wind_.t / DT, transpose(wind_.u))
plot!(xlabel="x", ylabel="time") |> display

# %% -------------- start model definition -------------------------
# %% Load Particle equations and derive ODE system

particle_equations = PW.particle_equations(u, u_x, γ=γ)
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


# condition(u, t, integrator) = 0.9* u[1] > log(17)
# affect!(integrator) = terminate!(integrator)
# cb = ContinuousCallback(condition, affect!)

ODE_settings = PW.ODESettings(
        Parameters=default_ODE_parameters,
        # define mininum energy threshold
        log_energy_minimum=log(FetchRelations.Eⱼ(0.1, DT)),
        #maximum energy threshold
        log_energy_maximum=log(17),  # correcsponds to Hs about 16 m
        saving_step=dt_ODE_save,
        timestep=DT,
        total_time=T,
        callbacks=nothing,
        save_everystep=false,
        adaptive           = true,
        dt                 = 1e-3, #60*10, 
        dtmin              = 1e-5, #60*5, 
        force_dtmin        = false,
)

# Default values for particle
#particle_defaults = ParticleDefaults(0.0, 1e-2, ODE_settings.log_energy_minimum) # old version
particle_defaults = ParticleDefaults(ODE_settings.log_energy_minimum, 1e-2, 0.0)

# %%%% Define callbacks
# # terminate    = PeriodicCallback(terminate_check!, 60*12   )
# show_mean    = PeriodicCallback(show_pos!      , DT/2    )
# periodic     = PeriodicCallback(periodic_BD_single_PI! ,  DT/2    )
# cbs          = CallbackSet(periodic, show_mean )#,cb_terminate)


Revise.retry()
vars = PW.init_vars_1D()
wave_model = WaveGrowthModels1D.WaveGrowth1D(; grid=grid1d,
        winds=u,
        ODEsys=particle_system,
        ODEvars=vars,
        layers=1,
        ODEsets=ODE_settings,  # ODE_settings
        ODEdefaults=particle_defaults,  # default_ODE_parameters
        periodic_boundary=periodic_boundary)

#using storing: init_state_store!, push_state_to_storage!
#using Simulations: Simulation

# %% run simulation
wave_simulation = Simulation(wave_model, Δt=20minutes, stop_time=5hours)
initialize_simulation!(wave_simulation, defaults=wave_model.ODEdefaults)
# or init_particles!( wave_model, defaults= wave_model.ODEdefaults, verbose = wave_simulation.verbose )

# using storing: init_state_store!
init_state_store!(wave_simulation, save_path)


run!(wave_simulation, store=true)


# %% Infinity run
reset_simulation!(wave_simulation)
# wave_simulation.stop_time = Inf
#run!(wave_simulation, store = false)

@info "... finished\n"

# %%

# %% open wave_simulation.store and plot variables
# convert HDF5 data to Array
store_waves_data = Array(wave_simulation.store.store["data"])
store_waves_energy = store_waves_data[:, :, 1]
store_waves_mx = store_waves_data[:, :, 2]
store_waves_my = store_waves_data[:, :, 3]
wave_x = wave_simulation.store.store["x"][:]
wave_time = wave_simulation.store.store["time"][:]
#print(wave_simulation.store.store["var_names"][:])

cg = store_waves_energy ./ store_waves_mx ./ 2

# %
#contourf(wave_x / dx, wave_time / DT, store_waves_my, levels=20, colormap=:blues)
#contourf(wave_x / dx, wave_time / DT, store_waves_mx, levels=20, colormap=:blues)
contourf(wave_x / dx , wave_time / DT, cg, levels=20, colormap=:blues)
#heatmap(wave_x / dx, wave_time / DT, store_waves_energy, levels=20, colormap=:blues)
plot!(xlabel="x (km)", ylabel="time (hours)") |> display




# %% plot failed particles, if there any

using ParticleTools
if ~isempty(wave_simulation.model.FailedCollection)
        @printf "plot failed particles\n"
        ParticleTools.PlotFailedParticles(wave_simulation.model.FailedCollection, Array(wave_simulation.store.store["data"]), ID, DT, dx,
                savepath=false, Npar=4)

end


### ------------ up to here ----------------------
# %% plotting errors if applicaple


# %% Statistics
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
forcing = (u=wind_.u, v=nothing)
coords = (x=wind_.x, time=wind_.t)

add_winds_forcing_to_store!(wave_simulation.store, forcing, coords)
show_stored_data(wave_simulation)
close_store!(wave_simulation)


# %% Save Particles
@printf "save particle \n"
save_object(joinpath(save_path, "particles.jld2"), wave_simulation.model.ParticleCollection)

#PC2 = load_object(joinpath(save_path , "particles.jld2"))

Nstate = size(wave_model.State)[2]


# %%
@printf "save parameters \n"
using JSON3
params = Dict("ID" => ID,
        "parset" => parset,
        "elapsed_time" => wave_simulation.run_wall_time,
        "" => Dict(
                "Nx" => Nx,
                "dx" => dx,
                "Nstate" => Nstate,
                "T" => T,
                "dt" => DT
        ),
        "ODE_sampling" => Dict(
                "dt_ODE_save" => dt_ODE_save
        ),
        "forcing" => Dict(
                "U10" => U10,
                "max_alpha" => U10 / (2 * maximum([maximum(cg[.!isnan.(cg)]), 0.01]))
        )
)

open(joinpath(save_path, "parameters.json"), "w") do io
        JSON3.pretty(io, params)
end

@printf "done.\n"

# %%
