import Plots

#using PiCLES.ParticleSystems: particle_waves_v3beta as PW3
using PiCLES.ParticleSystems: particle_waves_v5 as PW
import PiCLES: FetchRelations
using Setfield, IfElse

using PiCLES.Operators.core_1D: ParticleDefaults, InitParticleValues, InitParticleInstance
using PiCLES.ParticleMesh: OneDGrid, OneDGridNotes

using DifferentialEquations

using PiCLES.Utils.ParticleTools
using Plots

using Oceananigans.Units

# % Parameters
plot_path_base = "plots/tests/plot_path_base/"
mkpath(plot_path_base)

# function to define constants 

# %%

Revise.retry()

ODEpars, Const_ID, Const_Scg = PW.ODEParameters(r_g=0.85)
@info "org. gamma gamma= " Const_ID.γ
Const_ID.γ = Const_ID.γ #/ ODEpars.r_g
@info "gamma= " Const_ID.γ


# typeof(ODEpars)
T           = 1day
DT          = 5minutes
dt_ODE_save = 30minutes
u10         = 10
# define initial conditions
WindSeamin  = FetchRelations.get_initial_windsea(u10, DT)
#WindSeamin = FetchRelations.get_minimal_windsea( u(0, 0), 5minutes)

Revise.retry()

ODE_settings = PW.ODESettings(
    Parameters=ODEpars,
    # define mininum energy threshold
    log_energy_minimum=log(WindSeamin["E"]),#log(FetchRelations.Eⱼ(0.1, DT)),
    #maximum energy threshold
    log_energy_maximum=log(17),  # correcsponds to Hs about 16 m
    saving_step=5minutes,
    timestep=DT,
    total_time=T,
    save_everystep=false,
    dt=1e-3, #60*10, 
    dtmin=1e-9, #60*5, 
)

# define model -
using PiCLES: WaveGrowthModels1D, reset_boundary!
using Oceananigans.Units
using PiCLES.Simulations
using PiCLES.Plotting

Revise.retry()

u(x, t) = x .* 0 + t * 0 + u10

# redefine model 
wave_model = WaveGrowthModels1D.WaveGrowth1D(; grid=OneDGrid(0, 30e4, 50),
    winds=u,
    ODEsys=PW.particle_equations(u, γ=Const_ID.γ, q=Const_ID.q),
    ODEsets=ODE_settings,  # ODE_settings
    ODEinit_type="wind_sea",  # ODEpars
    minimal_particle=FetchRelations.MinimalParticle(2, 0, DT), #
    periodic_boundary=false,
    boundary_type="same"  # "default" #
)

plot_path_base = "plots/tests/T02_1D_fetch/with_merge_rule/"
#mkdir(plot_path_base)
@info "experiment 1: positive winds, periodic \n"

wave_simulation = Simulation(wave_model, Δt=DT, stop_time=10hours)
initialize_simulation!(wave_simulation)
run!(wave_simulation, store=false, cash_store=true, debug=false)

#Plotting.plot_results(wave_simulation, title="$u10 m/s, periodic=" * string(wave_model.periodic_boundary))
#title!("titlde")
#savefig(joinpath(plot_path_base, "PW4_u$(u10)_per_" * string(wave_model.periodic_boundary) * ".png"))


wave_simulation.store
data = Simulations.convert_store_to_tuple(wave_simulation.store, wave_simulation)
waves_energy = data.data[:, :, 1]

using Statistics
using PiCLES: FetchRelations as FR

Revise.retry()

data.x
# %
E_pic =  mean(waves_energy[end-4:1:end,: ], dims= 1)[1:end-1]
E_pic_tilde = FR.E_tilde(E_pic, u10)
# gr(display_type=:inline)
x_tilde = FR.X_tilde(data.x, u10)[1:end-1]
plot(x_tilde, E_pic_tilde, label="PiCLES")

# fetch = range(2, 2e6, step=1e3)
# x_tilde = FR.X_tilde(fetch, u10)[1:end-1]
E_jon_tilde = FR.E_fetch_tilde(x_tilde) 

plot!(x_tilde, E_jon_tilde, label = "JONSWAP")


# single particle
xx = OneDGridNotes(wave_model.grid).x[1]
particle_defaults = ParticleDefaults(log(WindSeamin["E"]), WindSeamin["cg_bar"], 0.0)
ParticleState, particle_on = InitParticleValues(particle_defaults, xx, u(xx, 0), DT)


PI4 = InitParticleInstance(wave_model.ODEsystem, ParticleState, ODE_settings, 1, false, particle_on)

function set_u_and_t!(integrator, u_new, t_new)
    integrator.u = u_new
    integrator.t = t_new
end


step!(PI4.ODEIntegrator, T, true)

PI = ParticleTools.FormatParticleData(PI4)

plot!(FR.X_tilde(PI.x, u10), FR.E_tilde(PI.E, u10), label="Sinlge Particle")

