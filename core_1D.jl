module core_1D

using SharedArrays
using DifferentialEquations

#push!(LOAD_PATH,   joinpath(pwd(), "code/")       )
using ParticleMesh

using particle_waves_v3: init_vars_1D
t, x, c̄_x, lne, r_g, C_α, g, C_e = init_vars_1D()

using ParticleInCell

using ParticleMesh: OneDGrid, OneDGridNotes
using ModelingToolkit: Num,  ODESystem

using FetchRelations

# initialize:
export init_z0_to_State

# Callbacks
export wrap_pos!, periodic_BD_single_PI!, show_pos!, periodic_condition_x

# Particle Node interaction
export GetParticleEnergyMomentum, GetVariablesAtVertex, Get_u_FromShared

export ParticleInstance

export InitParticleState,check_boundary_point


# ParticleInstance is the Stucture that carries each particle.
mutable struct ParticleInstance
        position_ij :: Int
        position_xy :: Float64
        ODEIntegrator::OrdinaryDiffEq.ODEIntegrator
        boundary :: Bool
end

# Debugging ParticleInstance
mutable struct MarkedParticleInstance
        Particle :: ParticleInstance
        time :: Float64
        state :: Vector{Any}
        errorReturnCode
end


Base.copy(s::ParticleInstance) = ParticleInstance(s.position_ij, s.position_xy, s.ODEIntegrator, s.boundary)



#### initialize #####
""" sets node state values to S at tuple position ij to zi """
function init_z0_to_State!(S::SharedMatrix, ij::Int, zi::Vector{Float64} )
        S[ ij[1], :] = zi
        nothing
end

""" sets particle state values to S. position is taking from particle """
function set_u_to_shared!(S::SharedMatrix, PI::ParticleInstance)
        S[ PI.position_ij[1], PI.position_ij[2] , :] = PI.ODEIntegrator.u
        nothing
end

""" takes node state values and pushes them to particle. position is taken from particle """
function set_u_to_particle!(S::SharedMatrix, PI::ParticleInstance)
        set_u!(PI.ODEIntegrator, Get_u_FromShared(S, PI) )
        u_modified!(PI.ODEIntegrator,true)
        nothing
end

##### Callbacks ######

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
function periodic_BD_single_PI!(PI::ParticleInstance, Lx::Float64)

        #@printf "PI periodic condition called"
        ui = copy(PI.ODEIntegrator.u)
        ui[1], wrap_pos_PI1 = wrap_pos!(ui[1], Lx)
        #ui[2], wrap_pos_PI2 = wrap_pos!(ui[2], Ly)

        if wrap_pos_PI1 #|| wrap_pos_PI1
                #@show wrap_pos_PI1 , wrap_pos_PI2
                #@printf "wrap pos PI"
                #@show PI.ODEIntegrator.u[1] - ui[1]
                #@show PI.ODEIntegrator.u[2] - ui[2]
                set_u!(PI.ODEIntegrator, ui )
                u_modified!(PI.ODEIntegrator,true)
        end
        nothing #return PI
end

# """
# Checks if
# """
# @everywhere function periodic_BD_single!(S::SharedMatrix)
#         S.u[1], wrap_pos_PI = wrap_pos!(S.u[1], Lx)
#         S.u[2], wrap_pos_PI = wrap_pos!(S.u[2], Ly)
#         nothing
# end

"""
return the mean position
"""
function show_pos!(PI)
        @show "show current position" PI.position_ij, PI.position_xy/dx, PI.ODEIntegrator.u[1]/dx, PI.ODEIntegrator.u[2]
end


""" (debubugging function) returns the domains normalized position of the partivle  """
function show_pos!(integrator , G)
        @show (integrator.ODEIntegrator.u[1] .- G.xmin)/grid1d.dx
end

# the following function are another way to have wrapping boundary conditions.
function periodic_condition_x(u,t,integrator)
        u[1] > 0
end

# function loop_x!(integrator)
#         integrator.u[1] = integrator.u[1] - Lx
#         @printf "wrap pos PI x"
# end
######### Particle --> Node ##########
"""
GetParticleEnergyMomentum(PI)

"""
function GetParticleEnergyMomentum(PI)
        ui_x, ui_c̄_x, ui_lne = PI.ODEIntegrator.u

        ui_e = exp(ui_lne)
        c_speed =  particle_waves_v3.speed(ui_c̄_x, ui_c̄_y)
        m_x = ui_e  / ui_c̄_x / 2
        #m_y = ui_c̄_y * ui_e  / c_speed^2 / 2

        return ui_e, m_x, 0
end

function GetParticleEnergyMomentum(PI::ParticleInstance)
        return GetParticleEnergyMomentum(PI.ODEIntegrator.u)
end

function GetParticleEnergyMomentum(zi::Dict)
        return GetParticleEnergyMomentum([ zi[x], zi[c̄_x], zi[lne]  ])
end

"""
GetParticleEnergyMomentum(PI)

"""
function GetParticleEnergyMomentum(z0::Vector{Float64})

        ui_x, ui_c̄_x, ui_lne = z0
        ui_e = exp(ui_lne)
        #c_speed =  ui_c̄_x
        m_x = ui_e  / ui_c̄_x / 2
        #m_y = ui_c̄_y * ui_e  / c_speed^2 / 2

        return [ui_e, m_x, 0]
end

######### node--> Particle ##########
"""
GetVariablesAtVertex(i_State::Vector{Float64}, x::Float64)

"""
function GetVariablesAtVertex(i_State::Vector{Float64}, x::Float64)
        e, m_x, m_y = i_State
        #m_amp = particle_waves_v3.speed(m_x, m_y)
        c_x = e / 2  / m_x
        #c_y = m_y * e / (2  * m_amp^2)

        return [x, c_x, log(e)]
end

""" returns node state from shared array, given paritcle index """
Get_u_FromShared(PI::ParticleInstance, S::SharedMatrix, ) = S[ PI.position_ij[1], :]






###### seed particles #####


"""
InitParticleState(defaults:: Dict{Num, Float64}, i::Int64, gridnote::OneDGridNotes, u, DT)

Find initial conditions for particle. Used at the beginning of the experiment.
        inputs:
        defaults        Dict with variables of the state vector, these will be replaced in this function
        i               index of the grid point
        gridnote        grid to determine the position in the grid
        u               interp. fucntion with wind values
        DT              time step of model, used to determine fetch laws
"""
function InitParticleState(
        defaults:: Dict{Num, Float64},
        i::Int64,
        gridnote::OneDGridNotes,
        u, DT  )

        # initalize state based on state vector
        defaults[x] = gridnote.x[i]
        # take in local wind velocities
        u_init    = u(defaults[x], 0)

        # seed particle given fetch relations
        defaults[c̄_x] = FetchRelations.c_g_U_tau( abs(u_init) , DT )
        defaults[lne] = log(FetchRelations.Eⱼ( abs(u_init) , DT ))
        #@show defaults
        return defaults
end

"""
check_boundary_point(i, boundary, periodic_boundary)
checks where ever or not point is a boundary point.
returns Bool
"""
function check_boundary_point(i, boundary, periodic_boundary)
    return periodic_boundary ? false : (i in boundary)
end


"""
SeedParticle!(ParticleCollection ::Vector{Any}, State::SharedMatrix, i::Int64,
                particle_system::ODESystem, particle_defaults::Dict{Num, Float64}, ODE_defaults::Dict{Num, Float64},
                GridNotes, winds, DT:: Float64, Nx:: Int, boundary::Vector{Int}, periodic_boundary::Bool)

Seed Pickles to ParticleColletion and State
"""
function SeedParticle!(
        InitParticleInstance,
        ParticleCollection ::Vector{Any},
        State::SharedMatrix,
        i::Int64,
        particle_system::ODESystem,
        particle_defaults::Dict{Num, Float64},
        ODE_defaults::Dict{Num, Float64},
        GridNotes, # ad type of grid note
        winds,     # interp winds
        DT:: Float64,
        Nx:: Int,
        boundary::Vector{Int},
        periodic_boundary::Bool)

        # define initial condition
        z_i             = InitParticleState(particle_defaults, i, GridNotes, winds, DT)
        # check if point is boundary point
        boundary_point  =  check_boundary_point(i, boundary, periodic_boundary)

        # add initial state to State vector
        init_z0_to_State!(State, i,  GetParticleEnergyMomentum(z_i) )

        # Push Inital condition to collection
        push!(  ParticleCollection,
                InitParticleInstance(
                                particle_system,
                                z_i ,
                                ODE_defaults ,
                                i ,
                                boundary_point))
        nothing
end

# end of module
end
