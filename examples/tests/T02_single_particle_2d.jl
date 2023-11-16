using ModelingToolkit, DifferentialEquations
using Plots
using Setfield

push!(LOAD_PATH, joinpath(pwd(), "code/"))
push!(LOAD_PATH, joinpath(pwd(), "code/Core"))

using PiCLES.ParticleSystems: particle_waves_v4 as PW4

import PiCLES: FetchRelations, ParticleTools
using PiCLES.Operators.core_2D: ParticleDefaults, InitParticleVector, InitParticleInstance
using ParticleMesh: TwoDGrid, TwoDGridNotes
using Oceananigans.Units

plot_path_base = "plots/tests/T04_2D_single_particle/"
mkpath(plot_path_base)
# %% Parameters
@register_symbolic u(x, y, t)
@register_symbolic v(x, y, t)

U10, V10 =  + 10.0, + 10.0


# version 3
r_g0 = 0.85
# function to define constants for grouwth and dissipation
Const_ID = PW4.get_I_D_constant()
#@set Const_ID.γ = 0.88
Const_Scg = PW4.get_Scg_constants()

# %%
#u(x, y, t) = 0.01 - U10 * sin(t / (6 * 60 * 60 * 2π))
#v(x, y, t) = 0.01 - V10 * cos(t / (6 * 60 * 60 * 2π))

#u(x, y, t) = - U10 + 0.01 + x * 0 + y * 0 + t *0
#v(x, y, t) = + V10 + 0.01 + x * 0 + y * 0 + t *0

u(x, y, t) =   (U10 * cos(t / (3 * 60 * 60 * 2π)) + 0.1) + x * 0 + y * 0
v(x, y, t) = - (V10 * sin(t / (3 * 60 * 60 * 2π)) + 0.1) + x * 0 + y * 0

winds = (u=u, v=v)


# define variables based on particle equation
t, x, y, c̄_x, c̄_y, lne, Δn, Δφ_p, r_g, C_α, C_φ, g, C_e = PW4.init_vars()

particle_equations = PW4.particle_equations(u, v, γ=Const_ID.γ, q=Const_ID.q)
@named particle_system = ODESystem(particle_equations)

typeof(particle_system)

# define V4 parameters absed on Const NamedTuple:
default_ODE_parameters = Dict(r_g => r_g0, C_α => Const_Scg.C_alpha,
    C_φ => Const_ID.c_β, C_e => Const_ID.C_e)

# define simple callback
condition(u, t, integrator) = 0.9 * u[1] > log(17)
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition, affect!)

# define standard initial conditions
DT = 4hours
WindSeamin = FetchRelations.get_initial_windsea(u(0, 0, 0), v(0, 0, 0), DT/2)
#WindSeamin = FetchRelations.get_minimal_windsea(u(0, 0, 0), v(0, 0, 0), 20minutes)

ODE_settings = PW4.ODESettings(
    Parameters=default_ODE_parameters,
    # define mininum energy threshold
    log_energy_minimum=log(WindSeamin["E"]),
    #maximum energy threshold
    log_energy_maximum=log(17),  # correcsponds to Hs about 16 m
    saving_step=2minutes,
    timestep=DT,
    total_time=T=6days,
    callbacks=cb,
    save_everystep=false,
    maxiters=1e4,
    adaptive = true,
    dt = 10,#60*10, 
    dtmin = 1,#60*5, 
    force_dtmin=true,
)

grid = TwoDGrid(3, 3, 3, 3)
particle_defaults = ParticleDefaults(log(WindSeamin["E"]), WindSeamin["cg_bar_x"], WindSeamin["cg_bar_y"], 0.0, 0.0)

# initialize particle given the wind conditions:
#ParticleState = InitParticleVector(copy(particle_defaults), (1, 1), TwoDGridNotes(grid), winds, DT)
ParticleState = copy(particle_defaults)
PI = InitParticleInstance(particle_system, ParticleState, ODE_settings, (0, 0), false, true)

states(particle_system)
#defaults(particle_system3)

function set_u_and_t!(integrator, u_new, t_new)
    integrator.u = u_new
    integrator.t = t_new
end

#clock_time = 0.0
for i in Base.Iterators.take(PI.ODEIntegrator, 25)
    #@info "t:", PI.ODEIntegrator.t
    #@info "u:", PI.ODEIntegrator.u
    #@info "i:", i

    @info "proposed dt", get_proposed_dt(PI.ODEIntegrator) / 60

    step!(PI.ODEIntegrator, DT, true)

    #@info "u:", PI.ODEIntegrator.u
    #clock_time += DT
    last_t = PI.ODEIntegrator.t

    ## define here the particle state at time of resetting
    #ui = [log(exp(PI.ODEIntegrator.u[1]) * 0.5), PI.ODEIntegrator.u[2] / 2, PI.ODEIntegrator.u[3] / 2, 0.0, 0.0]
    #ui = [lne_local, cg_u_local, cg_v_local, 0.0, 0.0]
    ui = PI.ODEIntegrator.u

    #ui = [PI.ODEIntegrator.u[1], PI.ODEIntegrator.u[2], PI.ODEIntegrator.u[3], 0.0, 0.0]
    WindSeamin = FetchRelations.get_initial_windsea(u(0, 0, last_t), v(0, 0, last_t), DT / 2)
    ui = [log(WindSeamin["E"]), WindSeamin["cg_bar_x"], WindSeamin["cg_bar_y"], 0.0, 0.0]

    set_u_and_t!(PI.ODEIntegrator, ui, last_t)
    # #set_u!(PI.ODEIntegrator, ui)
    #reinit!(PI.ODEIntegrator, ui, erase_sol=false, reset_dt=true, reinit_cache=true)
    #reinit!(PI3.ODEIntegrator, ui, erase_sol=false, reset_dt=true, reinit_cache=true)

    # #set_t!(PI.ODEIntegrator, last_t )
    u_modified!(PI.ODEIntegrator, true)

end

# %

PID = ParticleTools.ParticleToDataframe(PI)

gr(display_type=:inline)
# plit each row in PID and a figure

tsub = range(start=1, stop=length(PID[:, 1]), step=5)

subtitle="u=$U10 v=$V10 \n reset to windsea every $DT seconds\n"
p1 = plot(PID[tsub, 1] / (60 * 60), exp.(PID[tsub, 2]), marker=3, title=subtitle * "energy", xlabel="time (hours)", ylabel="e", label="V4") #|> display
p2 = plot(PID[tsub, 3], PID[tsub, 4], marker=3, markershape=:square, title="cg vector", xlabel="x", ylabel="y", label="V4") #|> display


# plot!(p1, PID[tsub, 1] / (60 * 60), FetchRelations.Eⱼ.(0.3 * abs(U10), PID[tsub, 1]) , marker=2, title="e", xlabel="time", ylabel="e", label="Fetch relations") #|> display

# plot!(p1, PID[tsub, 1] / (60 * 60), u.(0, 0, PID[tsub, 1]), marker=2, title="e", xlabel="time", ylabel="e", label="Fetch relations") #|> display
#plot!(p2, PID3[tsub3, 3], PID3[tsub3, 4], marker=2, title="cg vector", xlabel="x", ylabel="y", label="V3") #|> display

axlim = 10
plot!(p2, xlims=(-axlim, axlim), ylims=(-axlim, axlim))
plot!(p2, [0, 0], [-axlim, axlim], color=:black, linewidth=1, label=nothing)
plot!(p2, [-axlim, axlim], [0, 0], color=:black, linewidth=1, label=nothing)

tsubx = range(start=1, stop=length(PID[:, 1]), step=200)
time_sub = PID[tsubx, 1]
#plot quivers every qstep2
quiver!(p2, PID[tsubx, 3], PID[tsubx, 4], quiver=(u.(0, 0, time_sub) / 2, v.(0, 0, time_sub) / 2), color=:red, linewidth=2)#, label="wind")


#quiver!(p2, [0], [0], quiver=( [u.(0, 0, 0)], [v.(0, 0, 0)]), color=:red, linewidth=2, scale_units=:data, label="wind")

p3 = plot(PID[tsub, 5] / 1e3, PID[tsub, 6] / 1e3, marker=3, title="position", ylabel="postition", label="v4") #|> display

axlim = 200#1300
plot!(p3, xlims=(-axlim, axlim), ylims=(-axlim, axlim))
plot!(p3, [0, 0], [-axlim, axlim], color=:black, linewidth=1, label=nothing)
plot!(p3, [-axlim, axlim], [0, 0], color=:black, linewidth=1, label=nothing)

p4 = plot(PID[tsub, 5] / 1e3, exp.(PID[tsub, 2]), marker=3, title="e (x)", xlabel="x (km)", ylabel="e", label="V4") #|> display

plot(p1, p2, p3, p4, layout=(4, 1), legend=true, size=(600, 1600))


subtitle = "u$(U10)_v$(V10)_reset_to_windsea_dt$(DT)"
savefig(joinpath(plot_path_base, subtitle*"_continous_foreward.png"))




# # add legend ot p1
# plot!(p1, legend=:topleft)

