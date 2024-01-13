
using BenchmarkTools
#using Base.Threads

# CPU initialization
using Distributed
addprocs(2)
advance_workers = WorkerPool([2, 3])


@everywhere using Interpolations
@everywhere using IfElse
#@ redefining variables on all cores
@everywhere using ModelingToolkit, DifferentialEquations, Statistics

@everywhere push!(LOAD_PATH,   joinpath(pwd(), "MIC/")   )
@everywhere push!(LOAD_PATH,   joinpath(pwd(), "code/")   )
@everywhere push!(LOAD_PATH,   joinpath(pwd(), "analysis/")   )
#@everywhere import mesh
@everywhere using ParticleMesh
@everywhere using SharedArrays

#include("./particle_waves_v1.jl")
@everywhere import PIC
@everywhere import particle_waves_v3
@everywhere using Callbacks
@everywhere using Printf


@everywhere using ParticleMesh
using Revise
#using Plots
# %%


p_module = joinpath(pwd(), "code/particle_waves_v3.jl")
isfile(p_module)

@everywhere begin

        @register_symbolic u(x, y, t)
        @register_symbolic u_x(x, y, t)
        @register_symbolic v(x, y, t)
        @register_symbolic v_y(x, y, t)

        # %% create fake wind data and its gradients
        T = 60*60*24
        Lx, Ly = 50e3, 50e3
        xi = range(0,Lx,step=3e3)  # note ': this is a row vector
        yi = range(0,Ly,step=3e3)
        ti = range(0,T ,step=60*60)

        #u_func(x, y, t) = y * 0 .+ 3 * exp(- ( ( x-25e3 + 20e3* t/T )./10e3).^2) #.+ 3 #.+ 3 * sin.(x *π/Ly/0.5 .+ y *π/Ly/0.5)
        u_func(x, y, t) = y * 0 .+ 3 * exp(- ( ( x-25e3  )./10e3).^2) #.+ 3 #.+ 3 * sin.(x *π/Ly/0.5
        #@everywhere u_func(x, y, t) = y *0 .- x .*2/50e3 .+ 3 #.+ 3 * sin.(x *π/Ly/0.5 .+ y *π/Ly/0.5)
        #@everywhere u_func(x, y, t) = y *0 .- x .*0 .+ 3 #.+ 3 * sin.(x *π/Ly/0.5 .+ y *π/Ly/0.5)

        # @everywhere u_func(x, y, t) = y *0 .+ IfElse.ifelse.(    x .< Lx/2.0,
        #                                 2 .+ 1, *sin.(x *π/Lx/2),
        #                                 0 .* x + 0.1
        #                                 )

        v_func(x, y, t) = y *0 .+ x .*0 .+ 2 * y *π/Ly/0.5 * cos.(y *π/Ly/0.5)
        #v_func(x, y, t) = y * sin.(t * π ./T) * 2/Ly .+ x .*0 .+ 1* cos.( y *π/Ly/1 )
        #@everywhere v_func(x, y, t) = y *0 .+ x .*0 .+ 1* cos.(y *π/Ly/1  + t * 2 * pi /T )
        # @everywhere v_func(x, y, t) = x *0 .+ IfElse.ifelse.(    x .< Lx/2.0,
        #                                 2 * cos.(y *π/Ly/0.5) ,
        #                                 0 .* y + 0.1
        #                                 )


        # % define wind forcing globally as interp functions
        nodes = (xi, yi, ti)

        u_func_gridded = [ u_func(xii, yii, tii) for xii in xi, yii in yi, tii in ti]
        u_grid = LinearInterpolation( nodes , u_func_gridded ,extrapolation_bc=Periodic())
        #u_grid = CubicSplineInterpolation( nodes, u_func_gridded; bc=Line(OnGrid()), extrapolation_bc=Flat())
        #u_grid = interpolate(nodes, u_func_gridded, Gridded(Linear()))
        u(x, y, t) = u_grid(x, y, t)


        v_func_gridded = [ v_func(xii, yii, tii) for xii in xi, yii in yi, tii in ti]
        v_grid = LinearInterpolation( nodes,v_func_gridded ,extrapolation_bc=Periodic())
        #v_grid = CubicSplineInterpolation( nodes, v_func_gridded; bc=Line(OnGrid()), extrapolation_bc=Periodic())

        #v_grid = interpolate(nodes, v_func_gridded, Gridded(Linear()))
        #v_grid = extrapolate(  interpolate(v_func_gridded, BSpline(Linear()) ) , nodes)
        v(x, y, t) = v_grid(x, y, t)#v_grid(x, y)

end
# % plotting focring field
#using GLMakie


# create mesh
xmesh =  [ xii/1e3 for xii in xi, yii in yi]
ymesh =  [ yii/1e3 for xii in xi, yii in yi]

#arrows(    xi/1e3, yi/1e3 , u_func_gridded[:,:,1], v_func_gridded[:,:,1] ,arrowsize = 10, lengthscale = 2 )# , color = "black" )# |> display
#arrows(    xmesh, ymesh , u_func_gridded[:,:,1], v_func_gridded[:,:,1]  )# , color = "black" )# |> display


# plot 1st and last wind field
#arrows(    xmesh, ymesh , quiver = (u_func_gridded[:,:,1], v_func_gridded[:,:,1]) , color = "black" ) |> display

##############
# quiver(   xmesh, ymesh , quiver = (u_func_gridded[:,:,1], v_func_gridded[:,:,1]) , color = "black" ) |> display
#
# quiver!(   xmesh, ymesh , quiver = (u_func_gridded[:,:,Int(floor(length(ti)))], v_func_gridded[:,:,Int(floor(length(ti)/2))]) , color = "red" ) |> display
##############

#quiver!(   xmesh, ymesh , quiver = (u_func_gridded[:,:,Int(floor(length(ti)/2))], v_func_gridded[:,:,Int(floor(length(ti)/2))]) , color = "red" ) |> display
#        angles="xy", scale_units="width", scale=100 , headwidth=2, width= 0.005, headlength=3, color = "black", alpha= 0.5)

# %% define all components to solve ODE
@everywhere begin

        # only x gradient
        u_x(x, y, t) = Interpolations.gradient(u_grid, x, y, t )[1]

        # only y gradient
        v_y(x, y, t) = Interpolations.gradient(v_grid, x,y, t )[2]

        particle_equations0 = particle_waves_v3.particle_equations(u, v, u_x, v_y)

        @named particle_system0 = ODESystem(particle_equations0)

        # define variables based on particle equations
        t, x, y, c̄_x, c̄_y, lne, Δn, Δφ_p, r_g, C_α, C_φ, g, C_e = particle_waves_v3.init_vars()

        # ParticleInstance is the Stucture that carries each particle.
        mutable struct ParticleInstance
                position_ij :: Tuple
                position_xy :: Tuple
                ODEIntegrator::OrdinaryDiffEq.ODEIntegrator

                # dx::Float64
                # dy::Float64
                # function ParticleInstance(position_xy, ODE, dx::Float64, dy::Float64)
                #     new(  ( Int(round( position_xy[1]/dx)), Int(round( position_xy[2]/dy)) ) , position_xy  ,ODE)
                # end
        end

end


# %% define storing stucture and populate inital conditions


@everywhere Nx, Ny      = 30, 15
@everywhere grid2d      = ParticleMesh.TwoDGrid(1e3, Lx-1e3, Nx, 1e3, Ly-1e3 , Ny)
@everywhere grid2dnotes = ParticleMesh.TwoDGridNotes(grid2d)
@everywhere dx, dy      = grid2d.dx , grid2d.dy
@everywhere Nstate      = 3


Nparticle = grid2d.Nx* grid2d.Ny
State     = SharedArray{Float64}(grid2d.Nx, grid2d.Ny, Nstate)
State.pids # check which cores can see state


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


#index_positions, weights = PIC.compute_weights_and_index(grid2d, 12.34e3, 13.4e3)
#PIC.push_to_grid!(State, [1.0, 1.0,1.0,1.0,1.0,1.0,1.0]  , index_positions[1], weights[1], grid2d.Nx, grid2d.Ny )

# PIC.push_to_grid!(State, [1.0, 1.0,1.0,1.0,1.0,1.0,1.0]  , index_positions, weights, grid2d.Nx, grid2d.Ny )

# %%

# parameters
params0 = Dict(
        r_g => 0.9,
        C_α => -1.41 ,
        C_φ => 1.8e-5,
        g  => 9.81 ,
        C_e => 2.16e-4)

@everywhere begin

# define default initial state
e_0 = log(1e-4)
z0  = Dict( x => 0.0, y=> 0.0, c̄_x=> 1e-2, c̄_y=> -1e-9, lne=> e_0, Δn => dx, Δφ_p => 0.0)


# # copy State to shared array
# State[1, :,:] = [ z0[x] for i in range(0,length = Nparticle)  ]
#
# State[1, :,1] = [ z0[x] for i in range(0,length = Nparticle)  ]
# State[1, :,2] = [ 0 for i in y_range  ]


# @everywhere  function init_z0_to_State!(S::SharedArray, ij::Tuple{Int,Int}, zi::Vector{Float64} )
#         S[ ij[1], ij[2] , :] = [ zi[x],zi[y], zi[c̄_x],zi[c̄_y], zi[lne], zi[Δn],zi[Δφ_p]  ]#collect(values(zi))
#         nothing
# end
""" sets node state values to S at tuple position ij to zi """
function init_z0_to_State!(S::SharedArray, ij::Tuple{Int,Int}, zi::Vector{Float64} )
        S[ ij[1], ij[2] , :] = zi
        nothing
end

""" sets particle state values to S. position is taking from particle """
function set_u_to_shared!(S::SharedArray, PI::ParticleInstance)
        S[ PI.position_ij[1], PI.position_ij[2] , :] = PI.ODEIntegrator.u
        nothing
end

""" takes node state values and pushes them to particle. position is taken from particle """
function set_u_to_particle!(S::SharedArray, PI::ParticleInstance)
        set_u!(PI.ODEIntegrator, Get_u_FromShared(S, PI) )
        u_modified!(PI.ODEIntegrator,true)
        nothing
end

end
# %%%% Define callbacks

@everywhere begin

DT = 60*30
dt_save = DT/3.0

"""
wrap_pos!(xx::Float64, L::Float64)
makes periodic boundary conditions. If position xx exceeds 0 or L it is wrapped around
returns new position and if either or not the wrapping appends (bool)
"""
function wrap_pos!(xx::Float64, L::Float64)
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

"""
Checks if particle's PI position PI.ODEIntegrator.u[1,2] exceeeds the domain limits Lx, Ly
        Lx, LY are global variables
        If the PI's postions eceeeds the boundary they warpped around and the ODEIntegrator state is modified.
"""
function periodic_BD_single_PI!(PI::ParticleInstance)

        #@printf "PI periodic condition called"
        ui = copy(PI.ODEIntegrator.u)
        ui[4], wrap_pos_PI1 = wrap_pos!(ui[4], Lx)
        ui[5], wrap_pos_PI2 = wrap_pos!(ui[5], Ly)

        if wrap_pos_PI1 || wrap_pos_PI1
                #@show wrap_pos_PI1 , wrap_pos_PI2
                #@printf "wrap pos PI"
                #@show PI.ODEIntegrator.u[4] - ui[4]
                #@show PI.ODEIntegrator.u[5] - ui[5]
                set_u!(PI.ODEIntegrator, ui )
                u_modified!(PI.ODEIntegrator,true)
        end
        nothing #return PI
end


# """
# Checks if
# """
# @everywhere function periodic_BD_single!(S::SharedArray)
#         S.u[1], wrap_pos_PI = wrap_pos!(S.u[1], Lx)
#         S.u[2], wrap_pos_PI = wrap_pos!(S.u[2], Ly)
#         nothing
# end

"""
return the mean position
"""
function show_pos!(PI)
        @show "show current position" PI.position_xy, PI.ODEIntegrator.u[4], PI.ODEIntegrator.u[5]
end

# the following function are another way to have wrapping boundary conditions.

function periodic_condition_x(u,t,integrator)
        u[4] > 0
end

function periodic_condition_y(u,t,integrator)
        u[5] > Ly
end

# function loop_x!(integrator)
#         integrator.u[1] = integrator.u[1] - Lx
#         @printf "wrap pos PI x"
# end
#
# function loop_y!(integrator)
#         integrator.u[2] = integrator.u[2] - Ly
#         @printf "wrap pos PI y"
# end


end # end of everywhere

# terminate    = PeriodicCallback(terminate_check!, 60*12   )
show_mean    = PeriodicCallback(show_pos!      , DT/2    )
# set_to_mean  = PeriodicCallback(set_to_mean!    , 60*60*6 )
periodic     = PeriodicCallback(periodic_BD_single_PI! ,  DT/2    )
#
#periodic_x    = DiscreteCallback(periodic_condition_x, loop_x!)# , interp_points=0)
#periodic_y    = DiscreteCallback(periodic_condition_y, loop_y!)# , interp_points=0)

cbs           = CallbackSet(periodic, show_mean )#,cb_terminate)




######### Particle --> Node ##########
@everywhere begin

        """
        GetParticleEnergyMomentum(PI)

        """
        function GetParticleEnergyMomentum(PI)
                ui_x, ui_y, ui_c̄_x, ui_c̄_y, ui_lne, ui_Δn, ui_Δφ_p  = PI.ODEIntegrator.u

                ui_e = exp(ui_lne)
                c_speed =  particle_waves_v3.speed(ui_c̄_x, ui_c̄_y)
                m_x = ui_c̄_x * ui_e  / c_speed^2 / 2
                m_y = ui_c̄_y * ui_e  / c_speed^2 / 2

                return ui_e, m_x, m_y
        end

        function GetParticleEnergyMomentum(PI::ParticleInstance)
                return GetParticleEnergyMomentum(PI.ODEIntegrator.u)
        end

        function GetParticleEnergyMomentum(zi::Dict)
                return GetParticleEnergyMomentum([zi[lne], zi[c̄_x], zi[c̄_y], zi[x], zi[y], zi[Δn], zi[Δφ_p]])
        end

        """
        GetParticleEnergyMomentum(PI)

        """
        function GetParticleEnergyMomentum(z0::Vector{Float64})

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

        ######### node--> Particle ##########
        """
        GetVariablesAtVertex(i_State::Vector{Float64}, x::Float64, y::Float64, Δn::Float64, Δφ_p::Float64 )

        """
        function GetVariablesAtVertex(i_State::Vector{Float64}, x::Float64, y::Float64, Δn::Float64, Δφ_p::Float64 )
                e, m_x, m_y = i_State
                m_amp = particle_waves_v3.speed(m_x, m_y)
                c_x = m_x * e / (2  * m_amp^2)
                c_y = m_y * e / (2  * m_amp^2)

                return [log(e), c_x, c_y, x, y, Δn, Δφ_p]
        end
end
# testing
#GetVariablesAtVertex( [ie, imx, imy], 2e3, 4e3, dx, 0.0 ) - a_particle.ODEIntegrator.u


# %% seeding particles with initial states

"""
InitParticleInstance(model, z_initials, pars,  ij ; cbSets=nothing)
wrapper function to initalize a particle instance
        inputs:
        model           is an initlized ODESytem
        z_initials      is the initial state of the ODESystem
        pars            are the parameters of the ODESystem
        ij              is the (i,j) tuple that of the initial position
        chSet           (optional) is the set of callbacks the ODE can have
"""
function InitParticleInstance(model, z_initials, pars,  ij ; cbSets=nothing)

        # create ODEProblem
        problem    = ODEProblem(model, z_initials, (0.0,  T) , pars)
        # inialize problem
        #solver_method = Rosenbrock23()
        #solver_method = AutoVern7(Rodas4())
        solver_method = AutoTsit5(Rosenbrock23())
        #solver_method =Tsit5()
        integrator = init(problem, solver_method , saveat =dt_save , abstol = 1e-3, adaptive =true, maxiters=1e3)#, reltol=1e-1)
        #abstol=1e-4)#, reltol=1e-1)

        # callbacks= cbSets,
        #save_everystep=false
        return ParticleInstance( ij , (z_initials[x], z_initials[y]), integrator)
end


ParticleCollection=[]
for i in range(1,length = Nx), j in range(1,length = Ny)


        # initalize state based on state vector
        z_i = copy(z0)
        z_i[x] = grid2dnotes.x[i]
        z_i[y] = grid2dnotes.y[j]

        # take in local wind velocities
        u_init ,v_init       = u(z_i[x], z_i[y], 0) , v(z_i[x], z_i[y], 0)
        # derive local wave speed and direction
        # ! replace by wind sea formula in the future
        z_i[c̄_x]  = 1e-1 * u_init / sqrt(u_init.^2 + v_init.^2)
        z_i[c̄_y]  = 1e-1 * v_init / sqrt(u_init.^2 + v_init.^2)

        #@show (z_i[x]-grid2d.xmin)/grid2d.dx, (z_i[y]-grid2d.ymin)/grid2d.dy

        # push z_i to shared array
        init_z0_to_State!(State, (i,j),  GetParticleEnergyMomentum(z_i) )

        # copy initial parameters default
        params_i        = copy(params0)

        # push z_i and params_i to particle system
        # particle_system0 is an initilized ODESystem
        push!(ParticleCollection ,  InitParticleInstance(particle_system0, z_i , params_i, (i,j)   ))
        #@show threadid()
end

###### Debugging functions ############
@everywhere begin

        """ (debubugging function) returns the domains normalized position of the partivle  """
        function show_pos!(integrator , G)
                @show (integrator.ODEIntegrator.u[4] .- G.xmin)/grid2d.dx,  (integrator.ODEIntegrator.u[5] .- G.ymin)/grid2d.dy
        end


        """
                ResetParticle!(integrator)
        (debubugging function)
        Resets the integrator instance if the particle energy is nan or very high
        resets the particles position to the domain center
        resets the energy to a dfault value (e_0 is a global variable)

        """
        function ResetParticle!(integrator)
                if isnan(integrator.ODEIntegrator.u[1]) || exp(integrator.ODEIntegrator.u[1]) >= 1e-3
                        @show exp(integrator.ODEIntegrator.u[1])
                        integrator.ODEIntegrator.u[4] = integrator.ODEIntegrator.u[4] - Lx/2
                        integrator.ODEIntegrator.u[5] = integrator.ODEIntegrator.u[5] - Ly/2
                        #integrator.ODEIntegrator.u[4] = 1e-2
                        integrator.ODEIntegrator.u[1] = e_0
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
                if maximum(integrator.ODEIntegrator.u[4]) - Lx * Lx_terminate_limit >= 0 #|| maximum(exp.(integrator.u[3:N_state:end]) / e_0 ) >= 5
                        terminate!(integrator.ODEIntegrator)
                        @show "terminate"
                end
        end

end
#show_pos!(a_particle )
#advance(a_particle, State)

###### define remeshing routines ############
@everywhere begin

        """
                ParticleToNode!(PI::ParticleInstance, S::SharedArray, G::TwoDGrid)
        Pushes particle values to the neighboring nodes following the PIC rules.
        1.) get weights and indexes of the neighboring notes,
        2.) convert the particle state to nodestate
        3.) push the calculated state to the shared arrays

        inputs:

        PI      Particle instance
        S       Shared array where particles are stored
        G       (TwoDGrid) Grid that defines the nodepositions
        """
        function ParticleToNode!(PI::ParticleInstance, S::SharedArray, G::TwoDGrid)

                index_positions, weights = PIC.compute_weights_and_index(G, PI.ODEIntegrator.u[4], PI.ODEIntegrator.u[5])
                #ui[1:2] .= PI.position_xy
                #@show index_positions
                u_state = GetParticleEnergyMomentum(PI.ODEIntegrator.u)
                #@show u_state
                PIC.push_to_grid!(S, u_state , index_positions,  weights, G.Nx, G.Ny )
                nothing
        end

        """ returns node state from shared array, given paritcle index """
        Get_u_FromShared(PI::ParticleInstance, S::SharedArray, ) = S[ PI.position_ij[1], PI.position_ij[2] , :]
        #u_state = Get_u_FromShared(a_particle, State )

        """
                NodeToParticle!(PI::ParticleInstance, S::SharedArray)
        Pushes node value to particle:
        - If Node value is smaller than a minimal value, the particle is renintialized
        - If Node value is okey, it is converted to state variable and pushed to particle.
        - The particle position is set to the node positions
        """
        function NodeToParticle!(PI::ParticleInstance, S::SharedArray, ti::Number)
                u_state = Get_u_FromShared( PI, S)


                last_t = PI.ODEIntegrator.t
                if u_state[1] <= exp(e_0) # init new particle

                        # @show "re-init new particle"
                        # @show PI.position_ij, u_state
                        # z_i = copy(z0)
                        # z_i[x] = PI.position_xy[1]
                        # z_i[y] = PI.position_xy[2]
                        # ui = z_i

                        u_init ,v_init       = u(PI.position_xy[1], PI.position_xy[2], ti) , v(PI.position_xy[1], PI.position_xy[2], ti)
                        # derive local wave speed and direction
                        # ! replace by wind sea formula in the future
                        u_local = 1e-1 * u_init / sqrt(u_init.^2 + v_init.^2)
                        v_local = 1e-1 * v_init / sqrt(u_init.^2 + v_init.^2)

                        ui= [ PI.position_xy[1], PI.position_xy[2], u_local, v_local, z0[lne], z0[Δn], z0[Δφ_p] ]
                        reinit!(PI.ODEIntegrator, ui , erase_sol=false, reset_dt=true)#, reinit_callbacks=true)
                        #set_t!(PI.ODEIntegrator, last_t )
                        u_modified!(PI.ODEIntegrator,true)

                else    # load node value to particle

                        ui= GetVariablesAtVertex( u_state, PI.position_xy[1], PI.position_xy[2], z0[Δn], z0[Δφ_p] )
                        set_u!(PI.ODEIntegrator, ui )
                        #set_t!(PI.ODEIntegrator, last_t )
                        u_modified!(PI.ODEIntegrator,true)
                end
                #@show ui
                nothing
        end

        """
                advance!(PI::ParticleInstance, S::SharedArray)
        """
        function advance!(PI::ParticleInstance, S::SharedArray{Float64, 3}, G::TwoDGrid, DT::Int)
                #@show PI.position_ij

                add_saveat!(PI.ODEIntegrator, PI.ODEIntegrator.t )
                savevalues!(PI.ODEIntegrator)

                step!(PI.ODEIntegrator, DT , true)
                periodic_BD_single_PI!(PI )


                # step!(PI.ODEIntegrator, DT/2 , true)
                # PI = periodic_BD_single_PI!(PI )
                ParticleToNode!(PI, S, G)
                return PI
        end

        """
                remesh!(PI::ParticleInstance, S::SharedArray{Float64, 3})
        Wrapper function that does everything necessary to remesh the particles.
                - pushes the Node State to particle instance
        """
        function remesh!(PI::ParticleInstance, S::SharedArray{Float64, 3}, ti::Number)
                NodeToParticle!(PI, S, ti)
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

State_collect = []                 # storing array
State[:,:,:] .= 0                  # reset current state
push!(State_collect, copy(State) ) # push initial state to storage

u_sum_m1 = zeros(Nstate)

T/DT
DT
# main time stepping:
t_range = range(0.0, DT*T/DT, step=DT)



# %% single core
for t in t_range[2:end]

        @show t

        @show "advance"
        State[:,:,:] .= 0
        #u_sum_m1 = ShowTotalEnergyChange(ParticleCollection, u_sum_m1)

        for a_particle in ParticleCollection
                #@show a_particle.position_ij
                #try
                advance!(a_particle, State, grid2d, DT)
                # catch err
                #         @show "err", a_particle.ODEIntegrator.u
                # end

                # advance(a_particle, State)
                # # check callback
                # # show_pos!(  a_particle)
                # ResetParticle!(a_particle)
        end

        #u_sum_m1 = ShowTotalEnergyChange(ParticleCollection, u_sum_m1)

        #@show (State[:,:,1] .- grid2d.xmin)/grid2d.dx
        @show "re-mesh"
        for a_particle in ParticleCollection
                #@show a_particle.position_ij
                #show_pos!(  a_particle, grid2d)
                remesh!(a_particle, State, t)

                #periodic_BD_single!(a_particle )
                #show_pos!(  a_particle, grid2d)
        end

        @show "push"
        push!(State_collect, copy(State) )

        #u_sum_m1 = ShowTotalEnergyChange(ParticleCollection, u_sum_m1)
        #@show (State[:,:,1] .- grid2d.xmin)/grid2d.dx
end



# %% Parallel version

@everywhere carrier(f, y, g, d) = x -> f(x,y, g, d)
@everywhere carrier(f, y, t)    = x -> f(x,y, t)
global ParticleCollection

@time for t in t_range[2:end]

        @show t

        @show "advance"
        @time ParticleCollection = fetch(pmap(carrier( advance!, State, grid2d, DT) ,advance_workers,  ParticleCollection));

        @show "re-mesh"
        @time ParticleCollection = fetch(pmap(carrier( remesh!, State, t) ,advance_workers, ParticleCollection));
        #@show mean(State)

        @show "push"
        push!(State_collect, State )

end


# %%

heatmap( State_collect[2][:,:,1]' , levels = 20, color=:blues ) |> display


heatmap( State_collect[6][:,:,1]' , levels = 20, color=:blues ) |> display

# %%
heatmap( State_collect[10][:,:,1]' , levels = 20, color=:blues ) |> display




# %%

@printf "size of collected states %s" size(State_collect)
@printf "size of each state  %s" size(State_collect[1])
@printf "length of timerange  %s" size(t_range)


# %% saving files

@show "save data"
using HDF5

save_path = "data/first_try/"
mkpath(save_path)
rm(joinpath(save_path , "state.h5"), force=true)
file = h5open( joinpath(save_path , "state.h5") , "w" )

store_waves      = create_group(file, "waves")
store_waves_data = create_dataset(store_waves, "data", Float64, (length(State_collect),  (grid2dnotes.Nx), (grid2dnotes.Ny),  (size(State_collect[1])[3]) ) )#, chunk=(2,))

for i in 1:length(State_collect)
        store_waves_data[i,:,:,:] = State_collect[i]
end

#State_collect[1]


write_attribute(store_waves, "dims", ["time", "x", "y", "state"])
store_waves["time"] = Array(t_range)
store_waves["x"] = grid2dnotes.x
store_waves["y"] = grid2dnotes.y
store_waves["var_names"] = ["e", "m_x", "m_y"]


store_winds   = create_group(  file, "wind")
store_winds_u = create_dataset(store_winds, "u", Float64, ( (length(xi)), (length(yi)), (length(ti))  ) )#, chunk=(2,))
store_winds_v = create_dataset(store_winds, "v", Float64, ( (length(xi)), (length(yi)), (length(ti))  ) )#, chunk=(2,))

store_winds_u[:,:,:] = u_func_gridded
store_winds_v[:,:,:] = v_func_gridded


write_attribute(store_winds, "dims", ["x", "y", "time"])
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
@show "save particle"
using JLD2

save_object(joinpath(save_path , "particles.jld2"), ParticleCollection)
#PC2 = load_object(joinpath(save_path , "particles.jld2"))
