import Plots

#using PiCLES.ParticleSystems: particle_waves_v3beta as PW3
using PiCLES.ParticleSystems: particle_waves_v5 as PW
import PiCLES: FetchRelations
using Setfield, IfElse

using PiCLES.Operators.core_1D: ParticleDefaults, InitParticleValues, InitParticleInstance, GetVariablesAtVertex
using PiCLES.Operators.core_2D: ParticleDefaults as ParticleDefaults_2D, InitParticleValues as InitParticleValues_2D, InitParticleInstance as InitParticleInstance_2D

using PiCLES.ParticleMesh: OneDGrid, OneDGridNotes

using DifferentialEquations

using PiCLES.Utils.ParticleTools
using Plots

using Oceananigans.Units

using PiCLES: WaveGrowthModels1D, reset_boundary!
using PiCLES.Simulations
using PiCLES.Plotting

using Statistics
using PiCLES: FetchRelations as FR

# %
# % Parameters
plot_path_base = "plots/tests/plot_path_base/"
mkpath(plot_path_base)

# function to define constants 

Revise.retry()

ODEpars, Const_ID, Const_Scg = PW.ODEParameters(r_g=0.87)

@info "org. gamma gamma= " Const_ID.γ
Const_ID.γ = Const_ID.γ   #ODEpars.r_g^2
@info "gamma= " Const_ID.γ
# PW.e_T_func(Const_ID.γ, Const_ID.p, Const_ID.q, Const_ID.n; C_e=Const_ID.C_e, c_e=Const_ID.c_e, c_α=Const_ID.c_alpha)^2
# ODEpars.C_α
# ODEpars.C_φ
# ODEpars.C_e

Const_ID.c_alpha
# typeof(ODEpars)
T           = 2day
DT          = 10minutes
dt_ODE_save = 30minutes
u10         = 10.0
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
u(x, t) = x .* 0 + t * 0 + u10
# redefine model 
wave_model = WaveGrowthModels1D.WaveGrowth1D(; grid=OneDGrid(0, 1e6, 50),
    winds=u,
    ODEsys=PW.particle_equations(u, γ=Const_ID.γ, q=Const_ID.q, IDConstants=Const_ID),
    ODEsets=ODE_settings,  # ODE_settings
    ODEinit_type="wind_sea",  # ODEpars
    minimal_particle=FetchRelations.MinimalParticle(u10, 0, DT), #
    periodic_boundary=false,
    boundary_type="same"  # "default" #
)

plot_path_base = "plots/tests/T02_1D_fetch/with_merge_rule/"
#mkdir(plot_path_base)
@info "experiment 1: positive winds, periodic \n"

wave_simulation = Simulation(wave_model, Δt=DT, stop_time=T)
initialize_simulation!(wave_simulation)
run!(wave_simulation, store=false, cash_store=true, debug=false)

#Plotting.plot_results(wave_simulation, title="$u10 m/s, periodic=" * string(wave_model.periodic_boundary))
#title!("titlde")
#savefig(joinpath(plot_path_base, "PW4_u$(u10)_per_" * string(wave_model.periodic_boundary) * ".png"))

# %
squeeze(a) = dropdims(a, dims=tuple(findall(size(a) .== 1)...))


wave_simulation.store
data         = Simulations.convert_store_to_tuple(wave_simulation.store, wave_simulation)
data_slice   = squeeze(mean(data.data[end-4:1:end, :, :], dims=1))

wave_energy  = data_slice[2:end-1, 1]
wave_mx      = data_slice[2:end-1, 2]

wave_cgbar   = wave_energy ./ 2.0 ./ wave_mx
wave_omega_p = ODEpars.r_g * 9.81 ./ (2 * wave_cgbar)
wave_fp      = wave_omega_p / (2 * pi)


E_pic       = wave_energy
fp_pic      = wave_fp
E_pic_tilde = FR.E_tilde(E_pic, u10)
# gr(display_type=:inline)
x_tilde     = FR.X_tilde(data.x, u10)[2:end-1]


# fetch = range(2, 2e6, step=1e3)
# x_tilde = FR.X_tilde(fetch, u10)[1:end-1]
E_jon_tilde = FR.E_fetch_tilde(x_tilde) 

# PM 64
PM = FR.PMlimits()

pE_tilde_pic = plot(x_tilde, E_pic_tilde, label="PiCLES")
plot!(pE_tilde_pic, x_tilde, E_jon_tilde, label="JONSWAP")
plot!(pE_tilde_pic, x_tilde, x_tilde * 0 .+ PM.E_tilde, label="PW64")

pE_tilde_pic_temp = plot(x_tilde, E_pic, label="PiCLES cgmax =" * string(maximum(wave_cgbar)) )

PM.f_p_tilde * 2 * pi

pFp_tilde_pic = plot(x_tilde, x_tilde * 0 .+ PM.f_p_tilde, label="PW64", xlims=(1e3, 1e5), ylim=(0.1, 0.4 ))
plot!(pFp_tilde_pic, x_tilde, FR.f_p_tilde(fp_pic, u10), label="PiCLES cgmax =" * string(maximum(wave_cgbar))) #, ylim = (0, 20))
#pFp_tilde_pic = plot(x_tilde, FR.fₘ_from_X_tilde.(u10, x_tilde), label="PiCLES")

Revise.retry()
cg_JON = FR.c_p_fetch.(x_tilde, u10) /2  # cg

function get_fp_from_cg(cg)
    omega_p = 9.81 ./ (2 * cg)
    return omega_p / (2 * pi)
end

fp_JON = get_fp_from_cg(cg_JON)

#fp_JON = FR.fₘ_from_X_tilde.(u10, x_tilde)

9.81 ./ (4 * pi * fp_JON[end])
plot!(pFp_tilde_pic, x_tilde, FR.f_p_tilde(fp_JON, u10), label="JONSWAP")


# single particle
xx = OneDGridNotes(wave_model.grid).x[1]
particle_defaults = ParticleDefaults(log(WindSeamin["E"]), WindSeamin["cg_bar"], 0.0)
ParticleState, particle_on = InitParticleValues(particle_defaults, xx, u(xx, 0), DT)
PI4 = InitParticleInstance(wave_model.ODEsystem, ParticleState, ODE_settings, 1, false, particle_on)


step!(PI4.ODEIntegrator, T, true)
PI = ParticleTools.FormatParticleData(PI4)
plot!(pE_tilde_pic, FR.X_tilde(PI.x, u10), FR.E_tilde(PI.E, u10), label="Single Particle 1D")


plot!(pFp_tilde_pic, FR.X_tilde(PI.x, u10), FR.f_p_tilde(get_fp_from_cg(PI.cgx / ODEpars.r_g), u10), label="Single Particle 1D")

# % single particle 2D
u(x, y, t) = x .* 0 + t * 0 + u10
v(x, y,  t) = x .* 0 + t * 0 + y *0 
WindSeamin = FetchRelations.get_initial_windsea(u10, 0,  DT)
particle_defaults = ParticleDefaults_2D(log(WindSeamin["E"]), WindSeamin["cg_bar_x"], WindSeamin["cg_bar_y"], 0.0, 0.0)
ParticleState, particle_on = InitParticleValues_2D(particle_defaults, (xx, xx), (u(xx, xx, 0),0 ), DT)
ODESystem2D = PW.particle_equations(u, v, γ=Const_ID.γ, q=Const_ID.q, IDConstants=Const_ID)
PI5 = InitParticleInstance_2D(ODESystem2D, ParticleState, ODE_settings, (1, 1), false, particle_on)

step!(PI5.ODEIntegrator, T, true)
PI2 = ParticleTools.FormatParticleData(PI5)
plot!(pE_tilde_pic, FR.X_tilde(PI2.x, u10), FR.E_tilde(PI2.E, u10), label="Single Particle 2D")
plot!(pFp_tilde_pic, FR.X_tilde(PI2.x, u10), FR.f_p_tilde(get_fp_from_cg(PI2.cgx / ODEpars.r_g), u10), label="Single Particle 2D")


plot(pE_tilde_pic,pFp_tilde_pic, layout=(1, 2), size=(800, 300), title="Energy Spectra")



# E_tilde_bulk = Const_ID.c_e * x_tilde.^Const_ID.p
# plot!(x_tilde, E_tilde_bulk, label="Bulk E_tilde")