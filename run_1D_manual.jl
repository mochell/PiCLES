
using Plots
using BenchmarkTools

# CPU initialization
using Distributed
#addprocs(2)
#advance_workers = WorkerPool([2, 3])

@everywhere using Interpolations
@everywhere using IfElse
@everywhere using ModelingToolkit, DifferentialEquations, Statistics


# Particle Model modules
@everywhere push!(LOAD_PATH,   joinpath(pwd(), "code/")       )
@everywhere using ParticleMesh: OneDGrid, OneDGridNotes
@everywhere using SharedArrays

#include("./particle_waves_v1.jl")
using Revise

#@everywhere import ParticleInCell
includet("ParticleInCell.jl")
includet("FetchRelations.jl")

@everywhere import particle_waves_v3
@everywhere using Printf

# modules for saveing:
using HDF5
using JLD2

"""

"""
# %%
using core_1D: init_z0_to_State!, wrap_pos!, periodic_BD_single_PI!, show_pos!, periodic_condition_x
using core_1D: ParticleInstance
using core_1D: GetParticleEnergyMomentum, GetVariablesAtVertex, Get_u_FromShared
# %%
@printf "Init Forcing Field\n"

save_path_base= "data/1D_static/"

@everywhere begin

        @register_symbolic u(x, t)
        @register_symbolic u_x(x, t)

        # create fake wind data and its gradients
        T       = 60*60*24 *3
        Lx      = 50e3
        xi      = range(0,Lx,step=5e3)  # note ': this is a row vector
        ti      = range(0,T ,step=60*60)

        DT          = 60 * 20 #60*30 # remeshing time stepping
        dt_ODE_save = 60 * 10 # 3 min

        U10     = 2
        peridic_boundary = true

        #u_func(x, y, t) = y * 0 .+ 3 * exp(- ( ( x-25e3 + 20e3* t/T )./10e3).^2) #.+ 3 #.+ 3 * sin.(x *π/Ly/0.5 .+ y *π/Ly/0.5)

        #u_func(x, t) = 3 * exp(- ( ( x-25e3  )./10e3).^2) .+ t *0 #.+ 3 #.+ 3 * sin.(x *π/Ly/0.5
        u_func(x, t) = x.*0+U10 + t *0 #3 * exp(- ( ( x-25e3  )./10e3).^2) .+ t *0 #.+ 3 #.+ 3 * sin.(x *π/Ly/0.5

        #@everywhere u_func(x, y, t) = y *0 .- x .*2/50e3 .+ 3 #.+ 3 * sin.(x *π/Ly/0.5 .+ y *π/Ly/0.5)
        #@everywhere u_func(x, y, t) = y *0 .- x .*0 .+ 3 #.+ 3 * sin.(x *π/Ly/0.5 .+ y *π/Ly/0.5)

        #u_func(x, t) = x *0 .+ IfElse.ifelse.( x .< Lx*2.0/3.0, 3, 0.1)

        # @everywhere u_func(x, y, t) = y *0 .+ IfElse.ifelse.(    x .< Lx/2.0,
        #                                 2 .+ 1, *sin.(x *π/Lx/2),
        #                                 0 .* x + 0.1
        #                                 )

        # % define wind forcing globally as interp functions
        nodes = (xi, ti)

        u_func_gridded = [ u_func(xii, tii) for xii in xi, tii in ti]
        u_grid = LinearInterpolation( nodes , u_func_gridded ,extrapolation_bc=Periodic())
        #u_grid = CubicSplineInterpolation( nodes, u_func_gridded; bc=Line(OnGrid()), extrapolation_bc=Flat())
        #u_grid = interpolate(nodes, u_func_gridded, Gridded(Linear()))
        u(x, t) = u_grid(x, t)
        # only x gradient
        u_x(x, t) = Interpolations.gradient(u_grid, x )[1]

        particle_equations0 = particle_waves_v3.particle_equations(u, u_x)

        @named particle_system0 = ODESystem(particle_equations0)

        # define variables based on particle equations
        t, x, c̄_x, lne, r_g, C_α, g, C_e = particle_waves_v3.init_vars_1D()
end
# %% define storing stucture and populate inital conditions

@everywhere Nx          = 15
@everywhere grid1d      = OneDGrid(1e3, Lx-1e3, Nx)
@everywhere grid1dnotes = OneDGridNotes(grid1d)
@everywhere dx          = grid1d.dx
@everywhere Nstate      = 3


Nparticle = grid1d.Nx
State     = SharedMatrix{Float64}(grid1d.Nx, Nstate)
State.pids # check which cores can see state

# parameters
params0 = Dict(
        r_g => 0.9,
        C_α => -1.41 ,
        g   => 9.81 ,
        C_e => 2.16e-4
        )

@everywhere begin

        # define mininum energy threshold
        e_0 = log(FetchRelations.Eⱼ( 0.5 , float(DT) ))
        z0  = Dict( x => 0.0 , c̄_x=> 1e-2, lne => e_0)

end
# %%%% Define callbacks

# terminate    = PeriodicCallback(terminate_check!, 60*12   )
show_mean    = PeriodicCallback(show_pos!      , DT/2    )
periodic     = PeriodicCallback(periodic_BD_single_PI! ,  DT/2    )
cbs          = CallbackSet(periodic, show_mean )#,cb_terminate)

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
# z_i[c̄_x] = FetchRelations.c_g_U_tau( u_init , float(DT) )
# z_i[lne] = log(FetchRelations.Eⱼ( u_init , float(DT) ))
#
# #z_i[lne] = log(0.01)
#
# params_i        = copy(params0)
# problem    = ODEProblem(particle_system0, z_i, (0.0,  T) , params_i)
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
@printf "Init Particles\n"

"""
InitParticle(model, z_initials, pars,  ij ; cbSets=nothing)
wrapper function to initalize a particle instance
        inputs:
        model           is an initlized ODESytem
        z_initials      is the initial state of the ODESystem
        pars            are the parameters of the ODESystem
        ij              is the (i,j) tuple that of the initial position
        chSet           (optional) is the set of callbacks the ODE can have
"""
function InitParticle(model, z_initials, pars,  ij, boundary_flag ; cbSets=nothing)

        # create ODEProblem
        problem    = ODEProblem(model, z_initials, (0.0,  T) , pars)
        # inialize problem
        #solver_method = Rosenbrock23()
        #solver_method = AutoVern7(Rodas4())
        solver_method = AutoTsit5(Rosenbrock23())
        #solver_method =Tsit5()
        integrator = init(problem, solver_method , saveat =dt_ODE_save , abstol = 1e-3, adaptive =true, maxiters=1e3)#, reltol=1e-1)
        #abstol=1e-4)#, reltol=1e-1)

        # callbacks= cbSets,
        #save_everystep=false
        return ParticleInstance( ij , z_initials[x], integrator, boundary_flag )
end

# define boundaries
boundary = [1 , Nx]

ParticleCollection=[]
local z_i
local u_init
for i in range(1,length = Nx)

        # initalize state based on state vector
        z_i = copy(z0)
        z_i[x] = grid1dnotes.x[i]

        # take in local wind velocities
        u_init    = u(z_i[x], 0)

        # seed particle given fetch relations
        z_i[c̄_x] = FetchRelations.c_g_U_tau( abs(u_init) , float(DT) )
        z_i[lne] = log(FetchRelations.Eⱼ( abs(u_init) , float(DT) ))

        #@show e_0 > z_i[lne]
        #@show (z_i[x]-grid1d.xmin)/grid1d.dx, (z_i[y]-grid1d.ymin)/grid1d.dy

        # push z_i to shared array
        init_z0_to_State!(State, i,  GetParticleEnergyMomentum(z_i) )

        # push z_i and params_i to particle system
        # particle_system0 is an initilized ODESystem

        if ~peridic_boundary
                boundary_points = (i in boundary)
        else
                boundary_points = false
        end

        push!(  ParticleCollection,
                InitParticle(
                        particle_system0,
                        z_i ,
                        copy(params0) ,
                        i ,
                        boundary_points
                        )
                )
        #@show threadid()
end


###### Debugging functions ############
@everywhere begin


        """
                ResetParticle!(integrator)
        (debubugging function)
        Resets the integrator instance if the particle energy is nan or very high
        resets the particles position to the domain center
        resets the energy to a dfault value (e_0 is a global variable)

        """
        function ResetParticle!(integrator)
                if isnan(integrator.ODEIntegrator.u[3]) || exp(integrator.ODEIntegrator.u[3]) >= 1e-3
                        @show exp(integrator.ODEIntegrator.u[3])
                        integrator.ODEIntegrator.u[1] = integrator.ODEIntegrator.u[1] - Lx/2
                        #integrator.ODEIntegrator.u[2] = integrator.ODEIntegrator.u[2] - Ly/2
                        #integrator.ODEIntegrator.u[4] = 1e-2
                        integrator.ODEIntegrator.u[3] = e_0
                        u_modified!(integrator.ODEIntegrator,true)
                        @show "rest particle"
                end
                nothing
        end

        Lx_terminate_limit = 1


        """
                TerminateCheckSingle!(integrator)
        (debubugging function)

        """
        function TerminateCheckSingle!(integrator)
                if maximum(integrator.ODEIntegrator.u[1]) - Lx * Lx_terminate_limit >= 0 #|| maximum(exp.(integrator.u[3:N_state:end]) / e_0 ) >= 5
                        terminate!(integrator.ODEIntegrator)
                        @show "terminate"
                end
        end
end

###### remeshing routines ############


"""
        ParticleToNode!(PI::ParticleInstance, S::SharedMatrix, G::TwoDGrid)
Pushes particle values to the neighboring nodes following the ParticleInCell rules.
1.) get weights and indexes of the neighboring notes,
2.) convert the particle state to nodestate
3.) push the calculated state to the shared arrays

inputs:

PI      Particle instance
S       Shared array where particles are stored
G       (TwoDGrid) Grid that defines the nodepositions
"""
function ParticleToNode!(PI::ParticleInstance, S::SharedMatrix, G::OneDGrid)

        index_positions, weights = ParticleInCell.compute_weights_and_index_2d(G, PI.ODEIntegrator.u[1])
        #ui[1:2] .= PI.position_xy
        #@show index_positions
        u_state = GetParticleEnergyMomentum(PI.ODEIntegrator.u)
        #@show u_state, index_positions, weights
        ParticleInCell.push_to_2d_grid!(S, u_state , index_positions,  weights, G.Nx ,  peridic_boundary)
        nothing
end



"""
        NodeToParticle!(PI::ParticleInstance, S::SharedMatrix)
Pushes node value to particle:
- If Node value is smaller than a minimal value, the particle is renintialized
- If Node value is okey, it is converted to state variable and pushed to particle.
- The particle position is set to the node positions
"""
function NodeToParticle!(PI::ParticleInstance, S::SharedMatrix, ti::Number, e_0::Number)
        u_state = Get_u_FromShared( PI, S)


        last_t = PI.ODEIntegrator.t
        # test if particle is below energy threshold, or
        #      if particle is at the boundary
        if (u_state[1] < exp(e_0)) | PI.boundary # init new particle

                @show "re-init new particle"
                # @show PI.position_ij, u_state
                # z_i = copy(z0)
                # z_i[x] = PI.position_xy[1]
                # z_i[y] = PI.position_xy[2]
                # ui = z_i

                u_init       = u(PI.position_xy[1], ti)
                # derive local wave speed and direction
                cg_local     = FetchRelations.c_g_U_tau( abs(u_init) , float(DT) )
                lne_local    = log(FetchRelations.Eⱼ( abs(u_init) , float(DT) ))

                ui= [ PI.position_xy[1], cg_local, lne_local ]
                reinit!(PI.ODEIntegrator, ui , erase_sol=false, reset_dt=true)#, reinit_callbacks=true)
                #set_t!(PI.ODEIntegrator, last_t )
                u_modified!(PI.ODEIntegrator,true)

        else    # load node value to particle

                #@show "get vertex variable"
                #@show "u_state", u_state
                ui= GetVariablesAtVertex( u_state, PI.position_xy[1] )
                #@show ui
                set_u!(PI.ODEIntegrator, ui )
                #set_t!(PI.ODEIntegrator, last_t )
                u_modified!(PI.ODEIntegrator,true)
        end
        #@show ui
        nothing

end

# %%
######### Core routines for advancing and remeshing
@everywhere begin

        """
                advance!(PI::ParticleInstance, S::SharedMatrix{Float64}, G::OneDGrid, DT::Int)
        """
        function advance!(PI::ParticleInstance, S::SharedMatrix{Float64}, G::OneDGrid, DT::Int)
                #@show PI.position_ij

                add_saveat!(PI.ODEIntegrator, PI.ODEIntegrator.t )
                savevalues!(PI.ODEIntegrator)

                step!(PI.ODEIntegrator, DT , true)
                #@show PI.ODEIntegrator

                # use this for periodic boundary conditions
                #periodic_BD_single_PI!(PI, Lx) #### !!!! not sure if the boundary condition has to be set there, it might have beeen set somewhere else as well ..

                # step!(PI.ODEIntegrator, DT/2 , true)
                # PI = periodic_BD_single_PI!(PI )
                if isnan(PI.ODEIntegrator.u[1])
                        @show "position is nan"
                        @show PI
                        PI.ODEIntegrator.u = [0,0,0]
                end

                ParticleToNode!(PI, S, G)
                return PI
        end

        """
                remesh!(PI::ParticleInstance, S::SharedMatrix{Float64, 3})
                Wrapper function that does everything necessary to remesh the particles.
                - pushes the Node State to particle instance
        """
        function remesh!(PI::ParticleInstance, S::SharedMatrix{Float64}, ti::Number)
                NodeToParticle!(PI, S, ti, e_0)
                return PI
        end


        """ shows total energy of the all particles """
        function ShowTotalEnergyChange(ParticleCollection, u_sum_m1)
                u_sum = zeros(Nstate)
                for a_particle in ParticleCollection
                        u_sum +=  GetParticleEnergyMomentum(a_particle.ODEIntegrator.u)
                end
                @show u_sum_m1 - u_sum
                return u_sum
        end
end


# %% single core
@printf "Start itteration\n"
# define state collector
State_collect = []                 # storing array
State[:,:,:] .= 0                  # reset current state # [ Energy, x momentum, not used]
push!(State_collect, copy(State) ) # push initial state to storage

# main time stepping:
t_range = range(0.0, DT*T/DT, step=DT)


for t in t_range[2:end]
        #t = t_range[2]
        @show t

        #@printf "advance"
        State[:,:] .= 0
        #u_sum_m1 = ShowTotalEnergyChange(ParticleCollection, u_sum_m1)

        for a_particle in ParticleCollection
                #@show a_particle.position_ij
                #try
                advance!(a_particle, State, grid1d, DT)
                # catch err
                #          @show "err", a_particle.ODEIntegrator.u
                # end

                # advance(a_particle, State)
                # # check callback
                #show_pos!(  a_particle)
                # ResetParticle!(a_particle)
        end

        #u_sum_m1 = ShowTotalEnergyChange(ParticleCollection, u_sum_m1)

        #@show (State[:,:,1] .- grid1d.xmin)/grid1d.dx
        #@printf "re-mesh"
        for a_particle in ParticleCollection
                #@show a_particle.position_ij
                #@show "------------------------------------"
                #show_pos!(  a_particle)#, grid1d)
                remesh!(a_particle, State, t)
                #show_pos!(  a_particle)#, grid1d)
                #periodic_BD_single!(a_particle )
                #show_pos!(  a_particle, grid1d)
        end

        @printf "push"
        push!(State_collect, copy(State) )

        #u_sum_m1 = ShowTotalEnergyChange(ParticleCollection, u_sum_m1)
        #@show (State[:,:,1] .- grid1d.xmin)/grid1d.dx
end

# %% Parallel version

# @everywhere carrier(f, y, g, d) = x -> f(x,y, g, d)
# @everywhere carrier(f, y, t)    = x -> f(x,y, t)
# global ParticleCollection
#
# @time for t in t_range[2:end]
#
#         @show t
#
#         @show "advance"
#         @time ParticleCollection = fetch(pmap(carrier( advance!, State, grid1d, DT) ,advance_workers,  ParticleCollection));
#
#         @show "re-mesh"
#         @time ParticleCollection = fetch(pmap(carrier( remesh!, State, t) ,advance_workers, ParticleCollection));
#         #@show mean(State)
#
#         @show "push"
#         push!(State_collect, State )
#
# end

# %% Checks

@printf "\n size of collected states %s" size(State_collect)
@printf "\n size of each state  %s" size(State_collect[1])
@printf "\n length of timerange  %s" size(t_range)


# %% Saving Fields

@printf "save data \n"

save_path = save_path_base*"dx"* @sprintf("%i",dx)*"_dt"*@sprintf("%i",DT)*"_U"*@sprintf("%i",U10)*"/"
mkpath(save_path)

@printf "created model run folder %s" save_path

rm(joinpath(save_path , "state.h5"), force=true)
file = h5open( joinpath(save_path , "state.h5") , "w" )

store_waves      = create_group(file, "waves")
store_waves_data = create_dataset(
                        store_waves, "data", Float64,
                        (length(State_collect),  (grid1dnotes.Nx) ,  (size(State_collect[1])[2]) )
                        )#, chunk=(2,))

for i in 1:length(State_collect)
        store_waves_data[i,:,:] = State_collect[i]
end

#State_collect[1]


write_attribute(store_waves, "dims", ["time", "x", "state"])
store_waves["time"]      = Array(t_range)
store_waves["x"]         = grid1dnotes.x
store_waves["var_names"] = ["e", "m_x", "m_y"]


store_winds   = create_group(  file, "wind")
store_winds_u = create_dataset(store_winds, "u", Float64, ( (length(xi)), (length(ti))  ) )#, chunk=(2,))
store_winds_v = create_dataset(store_winds, "v", Float64, ( (length(xi)), (length(ti))  ) )#, chunk=(2,))

store_winds_u[:,:] = u_func_gridded
#store_winds_v[:,:,:] = v_func_gridded


write_attribute(store_winds, "dims", ["x", "time"])
store_winds["time"] = Array(ti)
store_winds["x"]    = Array(xi)
#store_winds["y"] = Array(yi)

#attr = create_attribute(file, "dims", ["time", "x", "y", "state"])
#HDF5.get_create_properties(state_store)

close(file)

# delete_object(parent, name)   # for groups, datasets, and datatypes
# delete_attribute(parent, name)   # for attributes

# %% Save Particles
@printf "save particle \n"
save_object(joinpath(save_path , "particles.jld2"), ParticleCollection)
#PC2 = load_object(joinpath(save_path , "particles.jld2"))

# %%
@printf "save parameters \n"
using JSON3
params = Dict("grid" => Dict(
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
                "U10" => U10
        )
        )

open( joinpath(save_path , "parameters.json") , "w") do io
    JSON3.pretty(io, params)
end

@printf "done.\n"
