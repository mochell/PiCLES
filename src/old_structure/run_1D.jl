
#@show localARGS

using Plots
using BenchmarkTools

using Interpolations
using IfElse
using ModelingToolkit, DifferentialEquations, Statistics


# Particle Model modules
push!(LOAD_PATH, joinpath(pwd(), "code/"))
using ParticleMesh: OneDGrid, OneDGridNotes
using SharedArrays

#include("./particle_waves_v1.jl")
using Revise

import ParticleInCell
import FetchRelations
#includet("ParticleInCell.jl")
#includet("FetchRelations.jl")

import particle_waves_v3
using Printf

using HDF5
using JLD2


"""
This load a parameter file, executes a 1D run, saves the data and run statistics.

"""

push!(LOAD_PATH, joinpath(pwd(), "code/"))
push!(LOAD_PATH, joinpath(pwd(), "code/Core"))
using core_1D_old: init_z0_to_State!, wrap_pos!, periodic_BD_single_PI!, show_pos!, periodic_condition_x
using core_1D_old: GetParticleEnergyMomentum, GetVariablesAtVertex, Get_u_FromShared
#using core_1D: InitParticleValues, check_boundary_point
using core_1D_old: InitParticleValues, SeedParticle!


using custom_structures: ParticleInstance1D, MarkedParticleInstance, AbstractParticleInstance

using InputOutput: Argsettings, parse_args

using ParticleTools: ParticleToDataframe

includet("Core/mapping_1D.jl")
#using mapping_1D
using Debugging


Revise.retry()

# Default values
save_path_base = "data/"
plot_path_base = "plots/static/"
parset = "1D_static/"

T = 24 * 3
Lx = 50

DT = 30 #60*30 # remeshing time stepping
dt_ODE_save = 10 # 3 min

Nx = 40

### Boundary Conditions
periodic_boundary = false

# parametric wind forcing
U10 = 20

# model parameters
r_g0 = 0.9
c_beta = 4#e-2 # 4e-2 Growth rate constant # default value from Kudr.
gamma = γ = 0.88 # dissipation wind energy input

@printf "Load passed arguments\n"
# arg_test = [    "--ID", "test1",
#                 "--Lx", "20",
#                 "--T",  "24",
#                 "--DT", "30"]

#localARGS = ["--ID", "Nx30_DT30_U15", "--Nx", "30", "--DT", "30", "--U10", "15", "--parset", "1D_static"]ß
#localARGS = ["--ID", "Nx50_DT2_U1", "--Nx", "50", "--c_beta", "2.000000", "--U10", "1", "--parset", "U10-c_beta", "--periodic"]
#localARGS = ["--ID", "Nx40_cbeta6.00_gamma1.32", "--c_beta", "6.00", "--Nx", "40", "--gamma", "1.32", "--parset", "gamma-c_beta"]
#localARGS = ["--ID", "rg0.85_cbeta6.00_gamma1.32", "--c_beta", "6.00", "--gamma", "1.32", "--rg", "0.85", "--parset", "gamma-c_beta-rg0.85"]
#localARGS = ["--ID", "rg1.25_gamma1.30_Nx40", "--Nx", "40", "--gamma", "1.32", "--rg", "1.25", "--parset", "gamma-rg0"]
#localARGS = ["--ID", "U25_DT60_Nx40", "--U10", "25", "--Nx", "40", "--DT", "60", "--parset", "U10-DT-periodic", "--periodic"]

#localARGS = ["--ID", "U20_DT3_Nx40", "--U10", "20", "--Nx", "40", "--DT", "3", "--parset", "U10-DT"]
#arg_case = ["--ID", "U25_DT60_Nx40", "--U10", "25", "--Nx", "40", "--DT", "60", "--parset", "U10-DT-periodic", "--periodic"]
#localARGS = ["--ID", "U13_DT60_Nx40", "--U10", "13", "--Nx", "40", "--DT", "60", "--parset", "U10-DT"]
localARGS = ["--ID", "gamma1.32_rg1.20_Nx50", "--gamma", "1.32", "--r_g0", "1.25", "--Nx", "50", "--parset", "gamma-rg0_2"]
@show localARGS
passed_argument = parse_args(localARGS, Argsettings)

# change here if more argument shuold be allowed to pass
#@unpack ID, parset, periodic, U10, Nx, DT = passed_argument #"U10-DT"
@unpack ID, parset, periodic, gamma, Nx, r_g0 = passed_argument #"gamma-rg_2"
#@unpack ID, parset, periodic, gamma, Nx, c_beta = passed_argument #"gamma-c_beta"

periodic_boundary = periodic
c_β = c_beta * 1e-2
C_e0 = (2.35 / r_g0) * 2e-3 * c_β
γ = Float64(gamma)
#r_g0 = rg
#
# if ~isnothing(ID)
#         #ID = "Nx"*@sprintf("%i",Nx)*"_dt"*@sprintf("%i",DT)*"_U"*@sprintf("%i",U10)
#         ID = "Nx"*@sprintf("%i",Nx)*"_cbeta"*@sprintf("%.1f",c_beta)*"_U"*@sprintf("%i",U10)
# end

#  rescale parameters for the right units.
T = T * 60 * 60 # seconds
Lx = Lx * 10e3  # km
DT = Float64(DT * 60) # seconds


# create ID and save name
#save_path = joinpath( "data/1D_static/", parsed_args["ID"] )
save_path = save_path_base * parset * "/" * ID * "/"
plot_path = plot_path_base * parset * "/"

# %%
@printf "Init Forcing Field\n"
# create wind test fucntion

@register_symbolic u(x, t)
@register_symbolic u_x(x, t)

#u_func(x, y, t) = y * 0 .+ 3 * exp(- ( ( x-25e3 + 20e3* t/T )./10e3).^2) #.+ 3 #.+ 3 * sin.(x *π/Ly/0.5 .+ y *π/Ly/0.5)
#u_func(x, t) = 3 * exp(- ( ( x-25e3  )./10e3).^2) .+ t *0 #.+ 3 #.+ 3 * sin.(x *π/Ly/0.5
u_func(x, t) = x .* 0 + U10 + t * 0 #3 * exp(- ( ( x-25e3  )./10e3).^2) .+ t *0 #.+ 3 #.+ 3 * sin.(x *π/Ly/0.5
#u_func(x, y, t) = y *0 .- x .*2/50e3 .+ 3 #.+ 3 * sin.(x *π/Ly/0.5 .+ y *π/Ly/0.5)
#u_func(x, y, t) = y *0 .- x .*0 .+ 3 #.+ 3 * sin.(x *π/Ly/0.5 .+ y *π/Ly/0.5)
#u_func(x, t) = x *0 .+ IfElse.ifelse.( x .< Lx*2.0/3.0, 3, 0.1)

# u_func(x, y, t) = y *0 .+ IfElse.ifelse.(    x .< Lx/2.0,
#                                 2 .+ 1, *sin.(x *π/Lx/2),
#                                 0 .* x + 0.1
#                                 )

# %% Define Grid, boundaries, and forcing field grid:

### Define grid
grid1d = OneDGrid(1e3, Lx - 1e3, Nx)
grid1dnotes = OneDGridNotes(grid1d)
dx = grid1d.dx

### storing matrix
Nstate = 3
Nparticle = grid1d.Nx
State = SharedMatrix{Float64}(grid1d.Nx, Nstate)

if ~periodic_boundary  # if false, define boundary points here:
        boundary = [1, Nx]
end

# provide winds in the right format
xi = range(0, Lx, step=5e3)  # note ': this is a row vector
ti = range(0, T, step=60 * 60)

# % define wind forcing globally as interp functions
u_func_gridded = [u_func(xii, tii) for xii in xi, tii in ti]
u_grid = linear_interpolation((xi, ti), u_func_gridded, extrapolation_bc=Periodic())
#u_grid = CubicSplineInterpolation( (xi, ti), u_func_gridded; bc=Line(OnGrid()), extrapolation_bc=Flat())
#u_grid = interpolate((xi, ti), u_func_gridded, Gridded(Linear()))


u(x, t) = u_grid(x, t)
# x gradient only
u_x(x, t) = Interpolations.gradient(u_grid, x)[1]

# %% Load Particle equations and derive ODE system

particle_equations0 = particle_waves_v3.particle_equations(u, u_x, γ=γ)
@named particle_system = ODESystem(particle_equations0)

# define variables based on particle equations
t, x, c̄_x, lne, r_g, C_α, g, C_e = particle_waves_v3.init_vars_1D()

# % define storing stucture and populate inital conditions

default_ODE_parameters = Dict(
        r_g => (1 / r_g0),
        C_α => -1.41,
        g => 9.81,
        C_e => C_e0,
)

### limiters:
# define mininum energy threshold
e_min_log = log(FetchRelations.Eⱼ(0.1, DT))
#e_min_log = log(FetchRelations.Eⱼ( 0.1 , 60.0 ))
e_max_log = log(17)  # correcsponds to Hs about 16 m

# Default values for particle
particle_defaults = Dict(x => 0.0, c̄_x => 1e-2, lne => e_min_log)

# %%%% Define callbacks

# terminate    = PeriodicCallback(terminate_check!, 60*12   )
show_mean = PeriodicCallback(show_pos!, DT / 2)
periodic = PeriodicCallback(periodic_BD_single_PI!, DT / 2)
cb = CallbackSet(periodic, show_mean)#,cb_terminate)

# testing
#GetParticleEnergyMomentum(z0)
#ui_x, ui_c̄_x, ui_lne  = PI.ODEIntegrator.u
#ie, imx, imy = GetParticleEnergyMomentum(PI.ODEIntegrator.u)

# testing
#GetVariablesAtVertex( [ie, imx, imy], ui_x) - PI.ODEIntegrator.u

# %% TEST 1 ODE solution
#
# z_i = copy(z0)
# z_i[x] = 20e4#grid1dnotes.x[20]
#
# ## seed particle given fetch relations
# u_init    = u(z_i[x], 0)
# z_i[c̄_x] = FetchRelations.c_g_U_tau( u_init , DT )
# z_i[lne] = log(FetchRelations.Eⱼ( u_init , DT ))
#
# #z_i[lne] = log(0.01)
#
# params_i        = copy(params0)
# problem    = ODEProblem(particle_system, z_i, (0.0,  T) , params_i)
# problem
#
# sol = solve(problem, saveat=60*5)
#
# solver_method = AutoTsit5(Rosenbrock23())
# #solver_method = Rosenbrock23()
# integrator = init(problem, solver_method , saveat =dt_ODE_save, adaptive =true, maxiters=1e3)#, reltol=1e-1)
# integrator
#
# step!(integrator, DT , true)
#
# # using PyPlot

# %%

# plt.figure()
# plt.subplot(3,1,1)
# plt.plot(sol[t], sol[c̄_x])
#
# plt.subplot(3,1,2)
# plt.plot(sol[t], sol[x].- sol[x][1])
#
# plt.subplot(3,1,3)
# plt.plot(sol[t], exp.(sol[lne]) )
#
# display(plt.gcf())
# %% seeding particles with initial states
@printf "seed particles ... \n"

# """
# InitParticleInstance(model, z_initials, pars,  ij ; cbSets=nothing)
# wrapper function to initalize a particle instance
#         inputs:
#         model           is an initlized ODESytem
#         z_initials      is the initial state of the ODESystem
#         pars            are the parameters of the ODESystem
#         ij              is the (i,j) tuple that of the initial position
#         chSet           (optional) is the set of callbacks the ODE can have
# """
# function InitParticleInstance(model, z_initials, pars,  ij, boundary_flag ; cbSets=nothing)

#         # create ODEProblem
#         problem    = ODEProblem(model, z_initials, (0.0,  T) , pars)

#         # inialize problem
#         #solver_method = Rosenbrock23()
#         #solver_method = AutoVern7(Rodas4())
#         solver_method = AutoTsit5(Rosenbrock23())
#         #solver_method =Tsit5()

#         # works best with abstol = 1e-4,reltol=1e-3,maxiters=1e4,
#         integrator = init(
#                         problem,
#                         solver_method,
#                         saveat =dt_ODE_save,
#                         abstol = 1e-4,
#                         adaptive =true,
#                         maxiters=1e4,
#                         reltol=1e-3)
#                         #abstol=1e-4)#,
#                         # callbacks= cbSets,
#                         #save_everystep=false
#         return ParticleInstance1D( ij , z_initials[x], integrator, boundary_flag )
# end


ODE_settings = particle_waves_v3.ODESettings(
        Parameters=particle_defaults,
        # define mininum energy threshold
        log_energy_minimum=e_min_log,#log(FetchRelations.Eⱼ(0.1, DT)),
        #maximum energy threshold
        log_energy_maximum=log(17),  # correcsponds to Hs about 16 m
        saving_step=dt_ODE_save,
        timestep=DT,
        total_time=T)


Revise.retry()
#seed particles
ParticleCollection = []
for i in range(1, length=Nx)
        SeedParticle!(ParticleCollection, State, i,
                particle_system, default_ODE_parameters, ODE_settings,
                grid1dnotes, u, DT, Nx, boundary, periodic_boundary)
end

# test
PI = ParticleCollection[1]
PI.ODEIntegrator
#sol = solve(PI.ODEIntegrator)
PI.ODEIntegrator.sol.u
step!(PI.ODEIntegrator, DT, true)
PI.ODEIntegrator.sol.u

#plot( PI.ODEIntegrator.sol.u[1], PI.ODEIntegrator.sol.u) 
using ParticleTools
PID = ParticleTools.ParticleToDataframe(PI)

# %%
# plit each row in PID and a figure
p1 = plot(PID[:, 2], exp.(PID[:, 4]), marker=2, title="e", xlabel="position", ylabel="e") #|> display
p2 = plot(PID[:, 2], PID[:, 3], marker=2, title="cg", ylabel="cg") #|> display
p3 = plot(PID[:, 1] / 60 / 60, PID[:, 4], marker=3, title="log e", xlabel="time", ylabel="postition") #|> display

plot(p1, p2, p3, layout=(3, 1), legend=false, size=(600, 1200))


# %%

# # #seed particles
# ParticleCollection=[]
# local z_i
# local u_init
# for i in range(1,length = Nx)
#         # there is SeedParticle! in the core module as alternative wrapper
#
#         # define initial condition
#         z_i = InitParticleValues(copy(z0), grid1dnotes, u, DT )
#         # check if point is boundary point
#         boundary_point  =  check_boundary_point(i, boundary, periodic_boundary)
#
#
#         #@show z_i, boundary_point
#
#         # add initial state to State vector
#         init_z0_to_State!(State, i,  GetParticleEnergyMomentum(z_i) )
#
#         # Push Inital condition to collection
#         push!(  ParticleCollection,
#                 InitParticleInstance(
#                                 particle_system,
#                                 z_i ,
#                                 copy(params0),
#                                 i ,
#                                 boundary_point))
#
# end



# %% single core
@printf "Start itteration ... \n"
# define state collector
State_collect = []                 # storing array
State[:, :, :] .= 0                  # reset current state # [ Energy, x momentum, not used]
push!(State_collect, copy(State)) # push initial state to storage

# collecting failed particles
FailedCollection = Vector{MarkedParticleInstance}([])

# main time stepping:
t_range = range(0.0, DT * T / DT, step=DT)

elapsed_time = @elapsed for t in t_range[2:end]

        State[:, :] .= 0
        #u_sum_m1 = mapping_1D.ShowTotalEnergyChange(ParticleCollection, u_sum_m1)

        for a_particle in ParticleCollection
                #@show a_particle.position_ij
                mapping_1D.advance!(a_particle, State, FailedCollection, grid1d, u, DT, e_min_log, e_max_log, periodic_boundary)
        end

        #u_sum_m1 = mapping_1D.ShowTotalEnergyChange(ParticleCollection, u_sum_m1)

        #@show (State[:,:,1] .- grid1d.xmin)/grid1d.dx
        #@printf "re-mesh"
        for a_particle in ParticleCollection
                mapping_1D.remesh!(a_particle, State, u, t, e_min_log, DT)
        end

        #@printf "push"
        push!(State_collect, copy(State))

        #u_sum_m1 = mapping_1D.ShowTotalEnergyChange(ParticleCollection, u_sum_m1)
        #@show (State[:,:,1] .- grid1d.xmin)/grid1d.dx
end

@printf "... finished\n"


# %% plotting errors if applicaple


if ~isempty(FailedCollection)
        @printf "plot failed particles\n"
        using ParticleTools
        mkpath(plot_path)

        ParticleTools.PlotFailedParticles(FailedCollection[1:min(4, size(FailedCollection)[1])], ID, DT, dx)
        savefig(joinpath(plot_path, "failed_ov_" * ID * ".png"))

        ParticleTools.plot_cg_version1(FailedCollection)
        savefig(joinpath(plot_path, "cg_failed_" * ID * ".png"))

        ParticleTools.plot_cg_version2(State_collect)
        savefig(joinpath(plot_path, "cg_failed2_" * ID * ".png"))

end

# %%
PI = ParticleCollection[10]
ParticleInCell.compute_weights_and_index(grid1d, PI.ODEIntegrator.u[3])

# ParticleCollection[10].ODEIntegrator.u
#
#
energy = State_collect[end-1][:, 1] #GP.sel(state= 0 ).T # energy
m_x = State_collect[end-1][:, 2] #GP.sel(state= 1 ).T # energy
cg = energy ./ m_x ./ 2 # c_g
# @printf "\n max wave age %s" U10 / (2 *  maximum( cg[.!isnan.(cg)] )  )
#


# %% Checks

@printf "\n size of collected states %s" size(State_collect)
@printf "\n size of each state  %s" size(State_collect[1])
@printf "\n length of timerange  %s" size(t_range)

# %% Saving Fields
@printf "save data \n"

mkpath(save_path)
@printf "created model run folder %s \n" save_path

rm(joinpath(save_path, "state.h5"), force=true)
file = h5open(joinpath(save_path, "state.h5"), "w")

store_waves = create_group(file, "waves")
store_waves_data = create_dataset(
        store_waves, "data", Float64,
        (length(State_collect), (grid1dnotes.Nx), (size(State_collect[1])[2]))
)#, chunk=(2,))

for i in 1:length(State_collect)
        store_waves_data[i, :, :] = State_collect[i]
end

write_attribute(store_waves, "dims", ["time", "x", "state"])
store_waves["time"] = Array(t_range)
store_waves["x"] = grid1dnotes.x
store_waves["var_names"] = ["e", "m_x", "m_y"]


store_winds = create_group(file, "wind")
store_winds_u = create_dataset(store_winds, "u", Float64, ((length(xi)), (length(ti))))#, chunk=(2,))
store_winds_v = create_dataset(store_winds, "v", Float64, ((length(xi)), (length(ti))))#, chunk=(2,))

store_winds_u[:, :] = u_func_gridded
#store_winds_v[:,:,:] = v_func_gridded

write_attribute(store_winds, "dims", ["x", "time"])
store_winds["time"] = Array(ti)
store_winds["x"] = Array(xi)
#store_winds["y"] = Array(yi)

close(file)


# %% Save Particles
# @printf "save particle \n"
# save_object(joinpath(save_path , "particles.jld2"), ParticleCollection)

#PC2 = load_object(joinpath(save_path , "particles.jld2"))

# %%
@printf "save parameters \n"
using JSON3
params = Dict("ID" => ID,
        "parset" => parset,
        "elapsed_time" => elapsed_time,
        "grid" => Dict(
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
