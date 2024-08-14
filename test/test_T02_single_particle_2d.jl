
include("dependencies_for_runtests.jl")
using DifferentialEquations
using PiCLES.ParticleSystems: particle_waves_v5 as PW
using PiCLES.Utils: Init_Standard

import PiCLES: FetchRelations, ParticleTools
using PiCLES.Operators.core_2D: InitParticleInstance
using Oceananigans.Units

# %%
u(x::Number, y::Number, t::Number) = (10.0 * cos(t / (3 * 60 * 60 * 2π)) + 0.1) + x * 0 + y * 0
v(x::Number, y::Number, t::Number) = -(10.0 * sin(t / (3 * 60 * 60 * 2π)) + 0.1) + x * 0 + y * 0

DT = 4hours
ParticleState, default_ODE_parameters, WindSeamin, Const_ID = Init_Standard(u(0.0, 0.0, 0.0), v(0.0, 0.0, 0.0), DT)
particle_system = PW.particle_equations(u, v, γ=Const_ID.γ, q=Const_ID.q)

ODE_settings = PW.ODESettings(
    Parameters=default_ODE_parameters,
    # define mininum energy threshold
    log_energy_minimum=log(WindSeamin["E"]),
    #maximum energy threshold
    log_energy_maximum=log(17),  # correcsponds to Hs about 16 m
    saving_step=2minutes,
    timestep=DT,
    total_time=T = 6days,
    save_everystep=false,
    maxiters=1e4,
    adaptive=true,
    dt=10,#60*10, 
    dtmin=1,#60*5, 
    force_dtmin=true)

PI = InitParticleInstance(particle_system, ParticleState, ODE_settings, (0, 0), false, true)


function set_u_and_t!(integrator, u_new, t_new)
    integrator.u = u_new
    integrator.t = t_new
end

function time_step_local!(PI, DT)
    "take 1 step over DT"

    step!(PI.ODEIntegrator, DT, true)
    last_t = PI.ODEIntegrator.t
    ui = PI.ODEIntegrator.u

    WindSeamin = FetchRelations.get_initial_windsea(u(0.0, 0.0, last_t), v(0.0, 0.0, last_t), DT / 2)
    ui = [log(WindSeamin["E"]), WindSeamin["cg_bar_x"], WindSeamin["cg_bar_y"], 0.0, 0.0]

    set_u_and_t!(PI.ODEIntegrator, ui, last_t)
    u_modified!(PI.ODEIntegrator, true)
    return PI
end

for i in Base.Iterators.take(PI.ODEIntegrator, 20)
    time_step_local!(PI, DT)
end


@testset "2D Single Particle reset" begin
    @testset "Assertions" begin
        PID2 = ParticleTools.FormatParticleData(PI)
        @test all(isfinite, PID2.x)
        @test all(isfinite, PID2.y)
        @test all(isfinite, PID2.cgx)
        @test all(isfinite, PID2.cgy)
        @test all(isfinite, PID2.E)
    end

    # test if PID2.E is bounded by the limits
    @testset "Energy bounds" begin
        PID2 = ParticleTools.FormatParticleData(PI)
        @test all(PID2.E .> log(WindSeamin["E"]/100))
        @test all(PID2.E .< log(17))
    end

    @testset "c_g bounds" begin
        PID2 = ParticleTools.FormatParticleData(PI)
        @test all(PID2.cgx .> -20)
        @test all(PID2.cgx .< 20)
        @test all(PID2.cgy .> -20)
        @test all(PID2.cgy .< 20)
    end
end


