module core_1D

using DifferentialEquations: OrdinaryDiffEq.ODEProblem, init
using ModelingToolkit: ODESystem

using SharedArrays
using ModelingToolkit: Num
using DocStringExtensions
# Particle-Node interaction
export GetParticleEnergyMomentum, GetVariablesAtVertex, Get_u_FromShared, ParticleDefaults
export InitParticleState, ResetParticleState


using FetchRelations

using Architectures: AbstractParticleInstance, AbstractMarkedParticleInstance
using ParticleMesh: OneDGrid, OneDGridNotes

using ..particle_waves_v3: init_vars_1D
t, x, c̄_x, lne, r_g, C_α, g, C_e = init_vars_1D()

using custom_structures: ParticleInstance1D, ParticleInstance2D, MarkedParticleInstance

export InitParticleInstance, init_z0_to_State!
include("initialize.jl")

# Callbacks
export wrap_pos!, periodic_BD_single_PI!, show_pos!, periodic_condition_x
include("utils.jl")


### particle defaults ###
"""
ParticleDefaults
Structure holds default particles
# Fields  
$(DocStringExtensions.FIELDS)
"""
struct ParticleDefaults
    "log energy"
    lne::Float64
    "horizontal velocity"
    c̄_x::Float64
    "x Position"
    x::Float64
end

# """
# ParticleDefaults
# Structure holds default particles
# # Fields  
# $(DocStringExtensions.FIELDS)
# """
# struct ParticleDefaults_old
#     "x Position"
#     x::Float64
#     "horizontal velocity"
#     c̄_x::Float64
#     "log energy"
#     lne::Float64
# end

Base.copy(s::ParticleDefaults) = Dict(lne => s.lne, c̄_x => s.c̄_x, x => s.x)

speed(x::Float64, y::Float64) = sqrt(x^2 + y^2)

######### Particle --> Node ##########

#### 1D code ####
# """
# GetParticleEnergyMomentum(PI)

# """
# function GetParticleEnergyMomentum(PI)
#     #ui_x, ui_c̄_x, ui_lne = PI.ODEIntegrator.u
#     ui_lne, ui_c̄_x, _ = PI.ODEIntegrator.u

#     ui_e = exp(ui_lne)
#     #c_speed = speed(ui_c̄_x, ui_c̄_y)
#     m_x = ui_e / ui_c̄_x / 2
#     #m_y = ui_c̄_y * ui_e  / c_speed^2 / 2

#     return ui_e, m_x, 0
# end

function GetParticleEnergyMomentum(PI::AbstractParticleInstance)
    return GetParticleEnergyMomentum(PI.ODEIntegrator.u)
end

function GetParticleEnergyMomentum(zi::Dict)
    return GetParticleEnergyMomentum([zi[lne], zi[c̄_x], zi[x]])
end

function GetParticleEnergyMomentum(zi::ParticleDefaults)
    return GetParticleEnergyMomentum([zi.lne, zi.c̄_x, zi.x])
end


"""
GetParticleEnergyMomentum(u_particle::Vector{Float64,3})
z0 is the particle state vector u with lne, cx, position
returns the state vector: energy, momentum, 0 (we are just in 1D) 
"""
function GetParticleEnergyMomentum(u_particle::Vector{Float64})

    ui_lne, ui_c̄_x, _ = u_particle
    ui_e = exp(ui_lne)
    #c_speed =  ui_c̄_x
    m_x = ui_e / ui_c̄_x / 2.0
    #m_y = ui_c̄_y * ui_e  / c_speed^2 / 2

    return [ui_e, m_x, 0]
end




######### node--> Particle ##########
"""
GetVariablesAtVertex(i_State::Vector{Float64}, x::Float64)
i_State: [e, m_x, m_y] state vector at node
x: coordinates of the vertex
"""
function GetVariablesAtVertex(i_State::Vector{Float64}, x::Float64)
    e, m_x, m_y = i_State
    #m_amp = particle_waves_v3.speed(m_x, m_y)
    c_x = e / 2.0 / m_x
    #c_y = m_y * e / (2  * m_amp^2)

    return [log(e), c_x, x]
end


""" returns node state from shared array, given particle index """
Get_u_FromShared(PI::AbstractParticleInstance, S::SharedMatrix,) = S[PI.position_ij[1], :]


###### seed particles #####


"""
InitParticleInstance(model::WaveGrowth1D, z_initials, pars,  ij ; cbSets=nothing)
wrapper function to initalize a particle instance
        inputs:
        model           is an initlized ODESytem
        z_initials      is the initial state of the ODESystem
        pars            are the parameters of the ODESystem
        ij              is the (i,j) tuple that of the initial position
        chSet           (optional) is the set of callbacks the ODE can have
"""
function InitParticleInstance(model, z_initials, ODE_settings, ij, boundary_flag; cbSets=Nothing)

    # create ODEProblem
    problem = ODEProblem(model, z_initials, (0.0, ODE_settings.total_time), ODE_settings.Parameters)
    # inialize problem
    # works best with abstol = 1e-4,reltol=1e-3,maxiters=1e4,
    integrator = init(
        problem,
        ODE_settings.solver,
        saveat=ODE_settings.saving_step,
        abstol=ODE_settings.abstol,
        adaptive=ODE_settings.adaptive,
        dt=ODE_settings.dt,
        dtmin=ODE_settings.dtmin,
        force_dtmin=ODE_settings.force_dtmin,
        maxiters=ODE_settings.maxiters,
        reltol=ODE_settings.reltol,
        callback=ODE_settings.callbacks,
        save_everystep=ODE_settings.save_everystep)
    return ParticleInstance1D(ij, z_initials[x], integrator, boundary_flag)
end



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
    defaults::PP,
    i::Int64,
    gridnote::OneDGridNotes,
    u, DT) where {PP<:Union{Dict,Nothing}}
    # take in local wind velocities

    if defaults == nothing
        #@info "init particles from fetch relations: $z_i"
        particle_defaults = Dict{Num,Float64}()

        u_init = u(gridnote.x[i], 0)
        # seed particle given fetch relations
        WindSeaMin = FetchRelations.get_initial_windsea(u_init, DT) # takes u_init just for the sign.
        particle_defaults[lne] = log(WindSeaMin["E"])
        particle_defaults[c̄_x] = WindSeaMin["cg_bar"]

    else
        particle_defaults = defaults#deepcopy(defaults)
    end

    # initalize state based on state vector
    particle_defaults[x] = gridnote.x[i]

    #@show defaults
    return particle_defaults
end

"""
ResetParticleState(defaults:: Dict{Num, Float64}, PI::AbstractParticleInstance, u_rn, DT)

resets the particle state to the default values if they given, otherwise it will use the fetch relations
        inputs:
        defaults        Dict with variables of the state vector, these will be replaced in this function
        PI              ParticleInstance to be reset
        u_rn            u value at the particle position
        DT              time step of model, used to determine fetch laws
"""
function ResetParticleState(
    defaults::PP,
    PI::AbstractParticleInstance,
    u_rn::Number, DT; vector=true) where {PP<:Union{Dict,Nothing}}
    # take in local wind velocities

    if defaults == nothing # this is boundary_defaults = "wind_sea"
        particle_defaults = Dict{Num,Float64}()
        #@info "init particles from fetch relations: $z_i"
        # seed particle given fetch relations
        WindSeaMin = FetchRelations.get_initial_windsea(u_rn, DT) # takes u_init just for the sign.
        particle_defaults[c̄_x] = WindSeaMin["cg_bar"]
        particle_defaults[lne] = log(WindSeaMin["E"])
    else
        particle_defaults = defaults#deepcopy(defaults)
    end

    # initalize state based on state vector
    particle_defaults[x] = PI.position_xy[1]

    if vector
        return [particle_defaults[lne], particle_defaults[c̄_x], particle_defaults[x]]
    else
        return particle_defaults
    end
    #@show defaults
    
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
    ParticleCollection::Vector{Any},
    State::SharedMatrix,
    i::Int64,
    particle_system::ODESystem,
    particle_defaults::PP,
    ODE_settings, #particle_waves_v3.ODESettings type
    GridNotes, # ad type of grid note
    winds,     # interp winds
    DT::Float64,
    boundary::Vector{T},
    periodic_boundary::Bool) where {T<:Union{Int,Any,Nothing,Int64},PP<:Union{Dict,Nothing}}

    # define initial condition
    z_i = InitParticleState(particle_defaults, i, GridNotes, winds, DT)

    # check if point is boundary point
    boundary_point = check_boundary_point(i, boundary, periodic_boundary)

    # add initial state to State vector
    init_z0_to_State!(State, i, GetParticleEnergyMomentum(z_i))

    # Push Inital condition to collection
    push!(ParticleCollection,
        InitParticleInstance(
            particle_system,
            z_i,
            ODE_settings,
            i,
            boundary_point))
    nothing
end


# end of module
end