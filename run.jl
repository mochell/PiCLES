
using BenchmarkTools
#using Base.Threads

# CPU initialization
using Distributed
addprocs(2)
advance_workers = WorkerPool([2])


@everywhere using Interpolations
@everywhere using IfElse
#@ redefining variables on all cores
@everywhere using ModelingToolkit, DifferentialEquations, Statistics

using Meshes

# %%

p_module = joinpath(pwd(), "code/particle_waves_v3.jl")
isfile(p_module)
@everywhere push!(LOAD_PATH,   joinpath(pwd(), "code/")   )
@everywhere push!(LOAD_PATH,   joinpath(pwd(), "analysis/")   )

#include("./particle_waves_v1.jl")
using Revise
@everywhere import particle_waves_v3
@everywhere using callbacks

@everywhere using Printf

@everywhere @register_symbolic u(x, y, t)
@everywhere @register_symbolic u_x(x, y, t)
@everywhere @register_symbolic v(x, y, t)
@everywhere @register_symbolic v_y(x, y, t)

# %% create fake wind data and its gradients

#using GLMakie
using Plots

@everywhere T = 60*60*24
@everywhere Lx, Ly = 50e3, 50e3
@everywhere xi = range(0,Lx,step=3e3)  # note ': this is a row vector
@everywhere yi = range(0,Ly,step=3e3)
@everywhere ti = range(0,T ,length=4)


@everywhere u_func(x, y, t) = y *0 .+ 3 * exp(- ( ( x-25e3 )./10e3).^2) #.+ 3 #.+ 3 * sin.(x *π/Ly/0.5 .+ y *π/Ly/0.5)
#@everywhere u_func(x, y, t) = y *0 .- x .*2/50e3 .+ 3 #.+ 3 * sin.(x *π/Ly/0.5 .+ y *π/Ly/0.5)
#@everywhere u_func(x, y, t) = y *0 .- x .*0 .+ 3 #.+ 3 * sin.(x *π/Ly/0.5 .+ y *π/Ly/0.5)

# @everywhere u_func(x, y, t) = y *0 .+ IfElse.ifelse.(    x .< Lx/2.0,
#                                 2 .+ 1, *sin.(x *π/Lx/2),
#                                 0 .* x + 0.1
#                                 )

#v_func(x, y, t) = y *0 .+ x .*0 .+ 2 * y *π/Ly/0.5 #* cos.(y *π/Ly/0.5)
@everywhere v_func(x, y, t) = y *0 .+ x .*0 .+ 1* cos.(y *π/Ly/1)
# @everywhere v_func(x, y, t) = x *0 .+ IfElse.ifelse.(    x .< Lx/2.0,
#                                 2 * cos.(y *π/Ly/0.5) ,
#                                 0 .* y + 0.1
#                                 )


@everywhere nodes = (xi, yi, ti)

@everywhere u_func_gridded = [ u_func(xii, yii, tii) for xii in xi, yii in yi, tii in ti]
@everywhere u_grid = LinearInterpolation( nodes , u_func_gridded ,extrapolation_bc=Periodic())
#u_grid = CubicSplineInterpolation( nodes, u_func_gridded; bc=Line(OnGrid()), extrapolation_bc=Flat())
#u_grid = interpolate(nodes, u_func_gridded, Gridded(Linear()))
@everywhere u(x, y, t) = u_grid(x, y, t)


@everywhere v_func_gridded = [ v_func(xii, yii, tii) for xii in xi, yii in yi, tii in ti]
@everywhere v_grid = LinearInterpolation( nodes,v_func_gridded ,extrapolation_bc=Periodic())
#v_grid = CubicSplineInterpolation( nodes, v_func_gridded; bc=Line(OnGrid()), extrapolation_bc=Periodic())

#v_grid = interpolate(nodes, v_func_gridded, Gridded(Linear()))
#v_grid = extrapolate(  interpolate(v_func_gridded, BSpline(Linear()) ) , nodes)
@everywhere v(x, y, t) = v_grid(x, y, t)#v_grid(x, y)

xmesh =  [ xii/1e3 for xii in xi, yii in yi]
ymesh =  [ yii/1e3 for xii in xi, yii in yi]

#arrows(    xi/1e3, yi/1e3 , u_func_gridded[:,:,1], v_func_gridded[:,:,1] ,arrowsize = 10, lengthscale = 2 )# , color = "black" )# |> display


#arrows(    xmesh, ymesh , u_func_gridded[:,:,1], v_func_gridded[:,:,1]  )# , color = "black" )# |> display


#arrows(    xmesh, ymesh , quiver = (u_func_gridded[:,:,1], v_func_gridded[:,:,1]) , color = "black" ) |> display
quiver(   xmesh, ymesh , quiver = (u_func_gridded[:,:,1], v_func_gridded[:,:,1]) , color = "black" ) |> display

#        angles="xy", scale_units="width", scale=100 , headwidth=2, width= 0.005, headlength=3, color = "black", alpha= 0.5)

# %%

# only x gradient
@everywhere u_x(x, y, t) = Interpolations.gradient(u_grid, x, y, t )[1]

# only y gradient
@everywhere v_y(x, y, t) = Interpolations.gradient(v_grid, x,y, t )[2]

particle_equations0 = particle_waves_v3.particle_equations(u, v, u_x, v_y)

@named particle_system0 = ODESystem(particle_equations0)

# define variables based on particle equations
@everywhere t, x, y, c̄_x, c̄_y, lne, Δn, Δφ_p, r_g, C_α, C_φ, g, C_e = particle_waves_v3.init_vars()

# ParticleInstance is the Stucture that carries each particle.
@everywhere mutable struct ParticleInstance
        position_ij :: Tuple
        position_xy :: Tuple
        ODEIntegrator::OrdinaryDiffEq.ODEIntegrator

        # dx::Float64
        # dy::Float64
        # function ParticleInstance(position_xy, ODE, dx::Float64, dy::Float64)
        #     new(  ( Int(round( position_xy[1]/dx)), Int(round( position_xy[2]/dy)) ) , position_xy  ,ODE)
        # end
end
# %% define storing stucture and populate inital conditions
@everywhere push!(LOAD_PATH,   joinpath(pwd(), "MIC/")   )
#@everywhere import mesh
@everywhere using mesh

@everywhere Nx, Ny = 30, 10
@everywhere grid2d = mesh.TwoDGrid(1e3, Lx-1e3, Nx, 10e3, Ly-10e3 , Ny)
@everywhere grid2dnotes = mesh.TwoDGridNotes(grid2d)

@everywhere using SharedArrays
@everywhere dx, dy = grid2d.dx , grid2d.dy

@everywhere Nstate    =    3
Nparticle = grid2d.Nx* grid2d.Ny
State     = SharedArray{Float64}(grid2d.Nx, grid2d.Ny, Nstate)
State.pids # check which cores can see state


@everywhere import PIC_2d_v0

# function wrap_index(pos::Int, N::Int)
#     if pos < 1
#         pos = pos + Int(N)
#         #pos =  N
#     elseif pos > N
#         pos = pos - Int(N)
#         #pos = N
#     end
#     return pos
# end

# function push_to_2d_grid_local!(grid::SharedArray{Float64, 3},
#                             charge::Vector{Float64},
#                             index_pos::Tuple{Int, Int},
#                             weights::Tuple{Float64, Float64},
#                             Nx::Int,  Ny::Int )
#
#         grid[ wrap_index(index_pos[1], Nx) , wrap_index(index_pos[2], Ny), : ] += weights[1] * weights[2] * charge
#         #grid[ index_pos[1] , index_pos[2] ] += weights[1] * weights[2] * charge
#
# end


#index_positions, weights = PIC_2d_v0.compute_weights_and_index_2d(grid2d, 12.34e3, 13.4e3)
#PIC_2d_v0.push_to_2d_grid!(State, [1.0, 1.0,1.0,1.0,1.0,1.0,1.0]  , index_positions[1], weights[1], grid2d.Nx, grid2d.Ny )

# PIC_2d_v0.push_to_2d_grid!(State, [1.0, 1.0,1.0,1.0,1.0,1.0,1.0]  , index_positions, weights, grid2d.Nx, grid2d.Ny )

# %%

# define default initial state
@everywhere e_0 = log(1e-4)
@everywhere z0  = Dict( x => 0.0, y=> 0.0, c̄_x=> 1e-2, c̄_y=> -1e-9, lne=> e_0, Δn => dx, Δφ_p => 0.0)

# parameters
params0 = Dict(
        r_g => 0.9,
        C_α => -1.41 ,
        C_φ => 1.8e-5,
        g  => 9.81 ,
        C_e => 2.16e-4)


# # copy State to shared array
# State[1, :,:] = [ z0[x] for i in range(0,length = Nparticle)  ]
#
# State[1, :,1] = [ z0[x] for i in range(0,length = Nparticle)  ]
# State[1, :,2] = [ 0 for i in y_range  ]



# @everywhere  function init_z0_to_State!(S::SharedArray, ij::Tuple{Int,Int}, zi::Vector{Float64} )
#         S[ ij[1], ij[2] , :] = [ zi[x],zi[y], zi[c̄_x],zi[c̄_y], zi[lne], zi[Δn],zi[Δφ_p]  ]#collect(values(zi))
#         nothing
# end

@everywhere  function init_z0_to_State!(S::SharedArray, ij::Tuple{Int,Int}, zi::Vector{Float64} )
        S[ ij[1], ij[2] , :] = zi
        nothing
end

@everywhere  function set_u_to_shared!(S::SharedArray, PI::ParticleInstance)
        S[ PI.position_ij[1], PI.position_ij[2] , :] = PI.ODEIntegrator.u
        nothing
end

@everywhere function set_u_to_particle!(S::SharedArray, PI::ParticleInstance)
        set_u!(PI.ODEIntegrator, get_u_from_shared(S, PI) )
        u_modified!(PI.ODEIntegrator,true)
        nothing
end

# %%%%% Define callbacks
@everywhere DT = 60*60
@everywhere dt_save = DT/3.0

@everywhere function wrap_pos!(xx::Float64, L::Float64)
        wrap_flag = false
        if xx > L
                xx = xx - L
                wrap_flag = true
        elseif xx < 0
                xx = xx + L
                wrap_flag = true
        end

        return xx, wrap_flag
end


@everywhere function periodic_BD_single_PI!(PI::ParticleInstance)

        #@printf "PI periodic condition called"
        ui = copy(PI.ODEIntegrator.u)
        ui[1], wrap_pos_PI1 = wrap_pos!(ui[1], Lx)
        ui[2], wrap_pos_PI2 = wrap_pos!(ui[2], Ly)

        if wrap_pos_PI1 || wrap_pos_PI1
                #@show wrap_pos_PI1 , wrap_pos_PI2
                #@printf "wrap pos PI"
                #@show PI.ODEIntegrator.u[1] - ui[1]
                #@show PI.ODEIntegrator.u[2] - ui[2]
                set_u!(PI.ODEIntegrator, ui )
                u_modified!(PI.ODEIntegrator,true)
        end
        nothing #return PI
end

@everywhere function periodic_BD_single!(S::SharedArray)
        S.u[1], wrap_pos_PI = wrap_pos!(S.u[1], Lx)
        S.u[2], wrap_pos_PI = wrap_pos!(S.u[2], Ly)
        nothing
end

function show_pos!(PI)
        """
        return the mean position
        """
        @show "show current position" PI.position_xy, PI.ODEIntegrator.u[1], PI.ODEIntegrator.u[2]
end

@everywhere function periodic_condition_x(u,t,integrator)
        u[1] > 0
end

@everywhere function periodic_condition_y(u,t,integrator)
        u[2] > Ly
end

@everywhere function loop_x!(integrator)
        integrator.u[1] = integrator.u[1] - Lx
        @printf "wrap pos PI x"
end

@everywhere function loop_y!(integrator)
        integrator.u[2] = integrator.u[2] - Ly
        @printf "wrap pos PI y"
end


# terminate    = PeriodicCallback(terminate_check!, 60*12   )
show_mean    = PeriodicCallback(show_pos!      , DT/2    )
# set_to_mean  = PeriodicCallback(set_to_mean!    , 60*60*6 )
periodic     = PeriodicCallback(periodic_BD_single_PI! ,  DT/2    )
#
periodic_x    = DiscreteCallback(periodic_condition_x, loop_x!)# , interp_points=0)
periodic_y    = DiscreteCallback(periodic_condition_y, loop_y!)# , interp_points=0)

cbs           = CallbackSet(periodic, show_mean )#,cb_terminate)


function init_particle(model, z_initials, pars,  ij ; cbSets=nothing)

        # create ODEProblem
        problem    = ODEProblem(model, z_initials, (0.0,  T) , pars)
        # inialize problem
        integrator = init(problem,Tsit5(), saveat =dt_save , abstol=1e-4)#, reltol=1e-1)
        # callbacks= cbSets,
        #save_everystep=false
        return ParticleInstance( ij , (z_initials[x], z_initials[y]), integrator)
end



######### Particle --> Vertex ##########
@everywhere function GetParticleEnergyMomentum(PI)
        ui_x, ui_y, ui_c̄_x, ui_c̄_y, ui_lne, ui_Δn, ui_Δφ_p  = PI.ODEIntegrator.u

        ui_e = exp(ui_lne)
        c_speed =  particle_waves_v3.speed(ui_c̄_x, ui_c̄_y)
        m_x = ui_c̄_x * ui_e  / c_speed^2 / 2
        m_y = ui_c̄_y * ui_e  / c_speed^2 / 2

        return ui_e, m_x, m_y
end

@everywhere function GetParticleEnergyMomentum(PI::ParticleInstance)
        return GetParticleEnergyMomentum(PI.ODEIntegrator.u)
end

@everywhere function GetParticleEnergyMomentum(zi::Dict)
        return GetParticleEnergyMomentum([ zi[x],zi[y], zi[c̄_x],zi[c̄_y], zi[lne], zi[Δn],zi[Δφ_p]  ])
end

@everywhere function GetParticleEnergyMomentum(z0::Vector{Float64})

        ui_x, ui_y, ui_c̄_x, ui_c̄_y, ui_lne, ui_Δn, ui_Δφ_p  = z0
        ui_e = exp(ui_lne)
        c_speed =  particle_waves_v3.speed(ui_c̄_x, ui_c̄_y)
        m_x = ui_c̄_x * ui_e  / c_speed^2 / 2
        m_y = ui_c̄_y * ui_e  / c_speed^2 / 2

        return [ui_e, m_x, m_y]
end

# testing
#GetParticleEnergyMomentum(z0)
#ui_x, ui_y, ui_c̄_x, ui_c̄_y, ui_lne, ui_Δn, ui_Δφ_p  = PI.ODEIntegrator.u
#ie, imx, imy = GetParticleEnergyMomentum(a_particle.ODEIntegrator.u)

######### Vertex --> Particle ##########

@everywhere  function GetVariablesAtVertex(i_State::Vector{Float64}, x::Float64, y::Float64, Δn::Float64, Δφ_p::Float64 )
        e, m_x, m_y = i_State
        m_amp = particle_waves_v3.speed(m_x, m_y)
        c_x = m_x * e / (2  * m_amp^2)
        c_y = m_y * e / (2  * m_amp^2)

        return [x, y,   c_x, c_y,   log(e), Δn, Δφ_p ]
end

# testing
#GetVariablesAtVertex( [ie, imx, imy], 2e3, 4e3, dx, 0.0 ) - a_particle.ODEIntegrator.u


# %% seeding particles

ParticleCollection=[]
for i in range(1,length = Nx), j in range(1,length = Ny)

        # initliza parameters
        params_i        = copy(params0)#init_particle()

        # initalize state based on state vector
        z_i = copy(z0)
        z_i[x] = grid2dnotes.x[i]
        z_i[y] = grid2dnotes.y[j]

        u_init ,v_init       = u(z_i[x], z_i[y], 0) , v(z_i[x], z_i[y], 0)
        z_i[c̄_x]  = 1e-1 * u_init / sqrt(u_init.^2 + v_init.^2)
        z_i[c̄_y]  = 1e-1 * v_init / sqrt(u_init.^2 + v_init.^2)

        #@show (z_i[x]-grid2d.xmin)/grid2d.dx, (z_i[y]-grid2d.ymin)/grid2d.dy

        #zi is set to shared array

        init_z0_to_State!(State, (i,j),  GetParticleEnergyMomentum(z_i) )
        # particle_system0 is an initilized ODESystem
        push!(ParticleCollection ,  init_particle(particle_system0, z_i , params_i, (i,j)   ))
        #@show threadid()
end


###### define remeshing routines ############

function show_pos!(integrator , G)
        """
        return the mean position
        """
        @show (integrator.ODEIntegrator.u[1] .- G.xmin)/grid2d.dx,  (integrator.ODEIntegrator.u[2] .- G.ymin)/grid2d.dy
end



Lx_terminate_limit = 1

@everywhere function reset_particle!(integrator)
        if isnan(integrator.ODEIntegrator.u[5]) || exp(integrator.ODEIntegrator.u[5]) >= 1e-3
                @show exp(integrator.ODEIntegrator.u[5])
                integrator.ODEIntegrator.u[1] = integrator.ODEIntegrator.u[1] - Lx/2
                integrator.ODEIntegrator.u[2] = integrator.ODEIntegrator.u[2] - Ly/2
                #integrator.ODEIntegrator.u[4] = 1e-2
                integrator.ODEIntegrator.u[5] = e_0
                u_modified!(integrator.ODEIntegrator,true)
                @show "rest particle"
        end
        nothing
end


@everywhere function terminate_check_single!(integrator)
        if maximum(integrator.ODEIntegrator.u[1]) - Lx * Lx_terminate_limit >= 0 #|| maximum(exp.(integrator.u[3:N_state:end]) / e_0 ) >= 5
                terminate!(integrator.ODEIntegrator)
                @show "terminate"
        end
end


#show_pos!(a_particle )
#advance(a_particle, State)

@everywhere function particle_to_cell!(PI::ParticleInstance, S::SharedArray, G::TwoDGrid)

        index_positions, weights = PIC_2d_v0.compute_weights_and_index_2d(G, PI.ODEIntegrator.u[1], PI.ODEIntegrator.u[2])
        #ui[1:2] .= PI.position_xy
        #@show index_positions
        u_state = GetParticleEnergyMomentum(PI.ODEIntegrator.u)
        #@show u_state
        PIC_2d_v0.push_to_2d_grid!(S, u_state , index_positions,  weights, G.Nx, G.Ny )
        nothing
end



@everywhere get_u_from_shared(PI::ParticleInstance, S::SharedArray, ) = S[ PI.position_ij[1], PI.position_ij[2] , :]
#u_state = get_u_from_shared(a_particle, State )


@everywhere function cell_to_particle!(PI::ParticleInstance, S::SharedArray)
        u_state = get_u_from_shared( PI, S)


        last_t = PI.ODEIntegrator.t
        if u_state[1] <= exp(e_0) # init new particle

                # @show "re-init new particle"
                # @show PI.position_ij, u_state
                # z_i = copy(z0)
                # z_i[x] = PI.position_xy[1]
                # z_i[y] = PI.position_xy[2]
                # ui = z_i
                ui= [ PI.position_xy[1], PI.position_xy[2], z0[c̄_x], z0[c̄_y], z0[lne], z0[Δn], z0[Δφ_p] ]
                reinit!(PI.ODEIntegrator, ui , erase_sol=false, reset_dt=true)#, reinit_callbacks=true)
                #set_t!(PI.ODEIntegrator, last_t )
                u_modified!(PI.ODEIntegrator,true)

        else    # load vertex value to particle

                ui= GetVariablesAtVertex( u_state, PI.position_xy[1], PI.position_xy[2], z0[Δn], z0[Δφ_p] )
                set_u!(PI.ODEIntegrator, ui )
                #set_t!(PI.ODEIntegrator, last_t )
                u_modified!(PI.ODEIntegrator,true)
        end
        #@show ui
        nothing
end


@everywhere function advance!(PI::ParticleInstance, S::SharedArray{Float64, 3}, G::TwoDGrid)
        #@show PI.position_ij

        add_saveat!(PI.ODEIntegrator, PI.ODEIntegrator.t )
        savevalues!(PI.ODEIntegrator)

        step!(PI.ODEIntegrator, DT , true)
        periodic_BD_single_PI!(PI )


        # step!(PI.ODEIntegrator, DT/2 , true)
        # PI = periodic_BD_single_PI!(PI )
        particle_to_cell!(PI, S, G)
        return PI
end

@everywhere function remesh!(PI::ParticleInstance, S::SharedArray{Float64, 3})
        cell_to_particle!(PI, S)
        return PI
end

function show_sum_change(ParticleCollection, u_sum_m1)
        u_sum = zeros(Nstate)
        for a_particle in ParticleCollection
                u_sum +=  GetParticleEnergyMomentum(a_particle.ODEIntegrator.u)
        end
        @show u_sum_m1 - u_sum
        return u_sum
end

State_collect = []

State[:,:,:] .= 0
push!(State_collect, copy(State) )


# #State[:,:,:] .= 0
u_sum_m1 = zeros(Nstate)

t_range = range(0.0, DT*10, step=DT)
for t in t_range
        @show t
        @show "advance"

        State[:,:,:] .= 0
        #u_sum_m1 = show_sum_change(ParticleCollection, u_sum_m1)

        for a_particle in ParticleCollection
                #@show a_particle.position_ij
                #try
                advance!(a_particle, State, grid2d)
                # catch err
                #         @show "err", a_particle.ODEIntegrator.u
                # end

                # advance(a_particle, State)
                # # check callback
                # # show_pos!(  a_particle)
                # reset_particle!(a_particle)
        end

        #u_sum_m1 = show_sum_change(ParticleCollection, u_sum_m1)

        #@show (State[:,:,1] .- grid2d.xmin)/grid2d.dx
        @show "re-mesh"
        for a_particle in ParticleCollection
                #@show a_particle.position_ij
                #show_pos!(  a_particle, grid2d)
                remesh!(a_particle, State)

                #periodic_BD_single!(a_particle )
                #show_pos!(  a_particle, grid2d)
        end

        push!(State_collect, copy(State) )
        #u_sum_m1 = show_sum_change(ParticleCollection, u_sum_m1)

        #@show (State[:,:,1] .- grid2d.xmin)/grid2d.dx
end
# %% Parallel version

# @everywhere carrier(f, y, g) = x -> f(x,y, g)
# @everywhere carrier(f, y) = x -> f(x,y)
#
# @time for t in range(0.0, DT*3, step=DT)
#         @show t
#         @show "advance"
#         @time ParticleCollection = fetch(pmap(carrier( advance!, State, grid2d) ,advance_workers,  ParticleCollection));
#
#         @show "re-mesh"
#         @time ParticleCollection = fetch(pmap(carrier( remesh!, State) ,advance_workers, ParticleCollection));
#         #@show mean(State)
# end
# %%
# using GLMakie
#
# aParticle  = ParticleCollection[1]
# a_matrix = aParticle.ODEIntegrator.sol.u[1]
#
# GetTimeIndex(xi, x) = argmin( abs.(x .-xi))
# #using GLMakie
# #using AbstractPlotting
#
# #Makie.inline!(true)
# #GLMakie.__init__()
# N_t     = length(aParticle.ODEIntegrator.sol.t)
# t_plot  = range(0.0, step=dt_save, length = N_t )
#
# for (tt,i) in zip(t_plot, range(1, stop=length(t_plot)+1  ) )
#
#         #for tt in aParticle.ODEIntegrator.sol.t
#
#         f = Figure()
#         ax = Axis(f[1, 1])
#
#         #aParticle  = ParticleCollection[i]
#         state_index =  GetTimeIndex(tt, t_range)
#         @show tt,i, state_index
#
#         aParticle =ParticleCollection[1]
#         for aParticle in ParticleCollection
#                 #a_matrix = mapreduce(permutedims, vcat, aParticle.ODEIntegrator.sol.u)
#                 a_matrix = aParticle.ODEIntegrator.sol.u[i]
#
#                 p_i = GLMakie.scatter!( [a_matrix[1]./1e3], [a_matrix[2]./1e3] ,markersize= 2e4* exp(a_matrix[5]) ,alpha=.2, color=:red)
#                 #p_i = plot!(a_matrix[:, 1]/1e3, a_matrix[:, 2]/1e3, linewidth = 3, markersize=8 ,alpha=.8,legend=false)
#                 #push!(plot_list, pi)
#                 #display(p_i)
#                 #@show aParticle.ODEIntegrator.u
#                 #@show aParticle.position_ij , exp(aParticle.ODEIntegrator.u[5])
#         end
#         xmesh =  [ xii/1e3 for xii in xi, yii in yi]
#         ymesh =  [ yii/1e3 for xii in xi, yii in yi]
#
#         #quiver!(    xmesh, ymesh , quiver = (u_func_gridded[:,:,1], v_func_gridded[:,:,1]), color = "black" ) |> display
#         #        angles="xy", scale_units="width", scale=100 , headwidth=2, width= 0.005, headlength=3, color = "black", alpha= 0.5)
#
#         GLMakie.arrows!(    xi/1e3, yi/1e3 , u_func_gridded[:,:,1], v_func_gridded[:,:,1] ,arrowsize = 10, lengthscale = 2 ) # |> display
#
#         hm = GLMakie.heatmap!(grid2dnotes.x/1e3, grid2dnotes.y/1e3, State_collect[state_index][:,:,1], levels = 20, color=:blues )
#
#         hm.colorrange = (0, maximum(State_collect[end][:,:,1])*1.1 )
#         #hm.levels = [0, 1e-3, 2e-3, maximum(State_collect[11][:,:,3])] #
#         #hm.levels  = [i for i in range(0, maximum(State_collect[11][:,:,3]) , length=10)]
#
#         ax.title = "time = $( tt/60/60), dt = $(aParticle.ODEIntegrator.sol.t[i]) "
#         display(f)
#
#         sleep(1/2) # refreshes the display!
#         #display()
# end
#
# %%   make Hs video with Makie

# Pi = ParticleCollection[1]
#
# t_plot = range(0.0, step=dt_save, length = N_t )
# #t_i_plot = range(1, stop=length(t_plot)+1  )
# tt  = Node(t_plot[1] ) # make dynamic variable
# i   = @lift(  Int( $tt/dt_save) +1  )
#
# state_index = @lift(  GetTimeIndex( $tt, t_range) )
#
# f = Figure()
# ax = Axis(f[1, 1])
#
# arrows!(    xi/1e3, yi/1e3 , u_func_gridded[:,:,1], v_func_gridded[:,:,1] ,arrowsize = 10, lengthscale = 2 ) # |> display
#
# data = @lift(State_collect[$state_index][:,:,1])
# hm = heatmap!(grid2dnotes.x/1e3, grid2dnotes.y/1e3, data , levels = 20, color=:blues )
#
# hm.colorrange = (0, maximum(State_collect[end][:,:,1])*1.1 )
# ax.title = "time = $( tt.val/60/60)"
# Makie.record(f, "plots/first_try/energy.mp4",  t_plot ,  framerate = 10) do tt_master
#         tt[] = tt_master
# end

# %% single particles
# #using Pyplot
# #plot()
# p1,p2 =plot(), plot()
# for a_particle in ParticleCollection
#
#         a_matrix = mapreduce(permutedims, vcat, a_particle.ODEIntegrator.sol.u)
#         time = a_particle.ODEIntegrator.sol.t/60/60
#
#
#         plot!(p1, time, sqrt.(a_matrix[:, 3].^2 + a_matrix[:, 4].^2)  ,alpha=.5,legend=false, color=:black)
#         #plot!(p1, time, a_matrix[:, 4] ,alpha=.5,legend=false)
#
#         energy = a_matrix[:, 5]
#         @show sum(isnan.(energy))
#         replace!(energy, NaN=>0)
#
#         plot!(p2, time, exp.(energy) ,alpha=.3,legend=false, color=:gray)
#
# end
#
# plot(p1, p2 , layout =   (2, 1)) |> display
#
# #3energy = a_matrix[:, 5]
# #replace!(energy, Inf=>NaN)


# %% saving files

using HDF5

save_path = "data/first_try/"
mkpath(save_path)
rm(joinpath(save_path , "state.h5"), force=true)
file = h5open( joinpath(save_path , "state.h5") , "w" )



store_waves = create_group(file, "waves")
store_waves_data = create_dataset(store_waves, "data", Float64, (length(State_collect),  (grid2dnotes.Nx), (grid2dnotes.Ny),  (size(State_collect[1])[3]) ) )#, chunk=(2,))

for i in 1:length(State_collect)
        store_waves_data[i,:,:,:] = State_collect[i]
end

#State_collect[1]

write_attribute(store_waves, "dims", ["time", "x", "y", "state"])
store_waves["time"] = Array(ti)
store_waves["x"] = grid2dnotes.x
store_waves["y"] = grid2dnotes.y
store_waves["var_names"] = ["e", "m_x", "m_y"]


store_winds = create_group(file, "wind")
store_winds_u = create_dataset(store_winds, "u", Float64, ((length(ti)), (length(xi)), (length(yi))  ) )#, chunk=(2,))
store_winds_v = create_dataset(store_winds, "v", Float64, ((length(ti)), (length(xi)), (length(yi)) ) )#, chunk=(2,))

store_winds_u[:,:,:] = u_func_gridded
store_winds_v[:,:,:] = v_func_gridded


write_attribute(store_winds, "dims", ["time", "x", "y"])
store_winds["time"] = Array(ti)
store_winds["x"] = Array(xi)
store_winds["y"] = Array(yi)


#State_collect[1]
#attr = create_attribute(file, "dims", ["time", "x", "y", "state"])
#HDF5.get_create_properties(state_store)


close(file)

# delete_object(parent, name)   # for groups, datasets, and datatypes
# delete_attribute(parent, name)   # for attributes

# %% save particles
using JLD2

save_object(joinpath(save_path , "particles.jld2"), ParticleCollection)
#PC2 = load_object(joinpath(save_path , "particles.jld2"))
