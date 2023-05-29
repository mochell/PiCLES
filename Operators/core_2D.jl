module core_2D

using DifferentialEquations: OrdinaryDiffEq.ODEIntegrator, OrdinaryDiffEq.ODEProblem, init
using ModelingToolkit: ODESystem

using SharedArrays
using ModelingToolkit: Num
using DocStringExtensions

# Particle-Node interaction
export GetParticleEnergyMomentum, GetVariablesAtVertex, Get_u_FromShared, ParticleDefaults, InitParticleInstance
export InitParticleState

#include("../Utils/FetchRelations.jl")
using FetchRelations

using Architectures: AbstractParticleInstance, AbstractMarkedParticleInstance
using ParticleMesh: TwoDGrid, TwoDGridNotes

using ..particle_waves_v3: init_vars
t, x, y, c̄_x, c̄_y, lne, Δn, Δφ_p, r_g, C_α, C_φ, g, C_e = init_vars()

using custom_structures: ParticleInstance1D, ParticleInstance2D, MarkedParticleInstance

export init_z0_to_State!
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
        "zonal velocity"
        c̄_x::Float64
        "meridional velocity"
        c̄_y::Float64
        "x Position"
        x::Float64
        "y Position"
        y::Float64
end

Base.copy(s::ParticleDefaults) = Dict(lne => s.lne, c̄_x => s.c̄_x, c̄_y => s.c̄_y, x => s.x, y => s.y)


######### Particle --> Node  |  2D  ##########


"""
GetParticleEnergyMomentum(PI)

"""
function GetParticleEnergyMomentum(PI)
        ui_lne, ui_c̄_x, ui_c̄_y, ui_x, ui_y = PI.ODEIntegrator.u

        ui_e = exp(ui_lne)
        c_speed = speed(ui_c̄_x, ui_c̄_y)
        m_x = ui_c̄_x * ui_e / c_speed^2 / 2
        m_y = ui_c̄_y * ui_e / c_speed^2 / 2

        return ui_e, m_x, m_y
end

function GetParticleEnergyMomentum(PI::AbstractParticleInstance)
        return GetParticleEnergyMomentum(PI.ODEIntegrator.u)
end

function GetParticleEnergyMomentum(zi::Dict)
        return GetParticleEnergyMomentum([zi[lne], zi[c̄_x], zi[c̄_y], zi[x], zi[y]])
end

function GetParticleEnergyMomentum(zi::ParticleDefaults)
        return GetParticleEnergyMomentum([zi[lne], zi[c̄_x], zi[c̄_y], zi[x], zi[y]])
end

"""
GetParticleEnergyMomentum(PI)

"""
function GetParticleEnergyMomentum(z0::Vector{Float64})

        ui_lne, ui_c̄_x, ui_c̄_y, ui_x, ui_y =z0
        ui_e = exp(ui_lne)
        c_speed = speed(ui_c̄_x, ui_c̄_y)
        m_x = ui_c̄_x * ui_e / c_speed^2 / 2
        m_y = ui_c̄_y * ui_e / c_speed^2 / 2

        return [ui_e, m_x, m_y]
end

speed(x::Float64, y::Float64) = sqrt(x^2 + y^2)

######### node--> Particle ##########

"""
GetVariablesAtVertex(i_State::Vector{Float64}, x::Float64, y::Float64, Δn::Float64, Δφ_p::Float64 )
i_state: [e, m_x, m_y] state vector at node
x, y: coordinates of the vertex

"""
function GetVariablesAtVertex(i_State::Vector{Float64}, x::Float64, y::Float64)
    e, m_x, m_y = i_State
    m_amp = speed(m_x, m_y)
    c_x = m_x * e / (2 * m_amp^2)
    c_y = m_y * e / (2 * m_amp^2)

    return [log(e), c_x, c_y, x, y]
end


""" returns node state from shared array, given paritcle index """
Get_u_FromShared(PI::AbstractParticleInstance, S::SharedArray,) = S[PI.position_ij[1], PI.position_ij[2], :]

"""
GetGroupVelocity(i_State::Vector{Float64})
i_state: [e, m_x, m_y] state vector at node
"""
function GetGroupVelocity(i_State)
        e = i_State[:,:,1]
        m_x = i_State[:,:,2]
        m_y = i_State[:,:,3]
        m_amp = speed.(m_x, m_y)
        c_x = m_x .* e ./ (2 .* m_amp.^2)
        c_y = m_y .* e ./ (2 .* m_amp.^2)

        return (c_x=c_x, c_y=c_y)
end


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
                dt = ODE_settings.dt,
                dtmin=ODE_settings.dtmin,
                force_dtmin=ODE_settings.force_dtmin,
                maxiters=ODE_settings.maxiters,
                reltol=ODE_settings.reltol,
                callback=ODE_settings.callbacks,
                save_everystep=ODE_settings.save_everystep)
        return ParticleInstance2D(ij, (z_initials[x], z_initials[y]), integrator, boundary_flag)
end



"""
InitParticleState(defaults:: Dict{Num, Float64}, ij::Tuple, gridnote::OneDGridNotes, u, v, DT)

Find initial conditions for particle. Used at the beginning of the experiment.
        inputs:
        defaults        Dict with variables of the state vector, these will be replaced in this function
        ij              index of the grid point
        gridnote        grid to determine the position in the grid
        winds           NamedTuple (u,v) with interp. functions for wind values
        DT              time step of model, used to determine fetch laws
"""
function InitParticleState(
        defaults::PP,
        ij::Tuple{Int, Int},
        xy::Tuple{Float64, Float64},
        uv::Tuple{Float64, Float64}, 
        DT) where {PP<:Union{Dict,Nothing}}

        i,j = ij
        xx, yy = xy
 
        if defaults == nothing

                #@info "init particles from fetch relations: $z_i"
                particle_defaults = Dict{Num,Float64}()
                # take in local wind velocities
                u_init, v_init = uv#winds.u(xx, yy, 0), winds.v(xx, yy, 0)

                WindSeaMin = FetchRelations.get_initial_windsea(u_init, v_init, DT)
                # seed particle given fetch relations
                particle_defaults[lne] = log(WindSeaMin["E"])
                particle_defaults[c̄_x] = WindSeaMin["cg_bar_x"]
                particle_defaults[c̄_y] = WindSeaMin["cg_bar_y"]

        else
                particle_defaults = defaults
        end
        particle_defaults[x] = xx
        particle_defaults[y] = yy


        #@show defaults
        return particle_defaults
end

"""
ResetParticleState(PI::AbstractParticleInstance, particle_defaults::Dict{Num,Float64}, ODE_settings::ODESettings, gridnote::OneDGridNotes, winds, DT)

Reset the state of a particle instance
        inputs:
        PI                      ParticleInstance
        particle_defaults       Dict with variables of the state vector, these will be replaced in this function
        ODE_settings            ODESettings type
        gridnote                grid to determine the position in the grid
        winds                   NamedTuple (u,v) with interp. functions for wind values
        DT                      time step of model, used to determine fetch laws
returns:
        dict or vector
"""
function ResetParticleState(defaults::PP,
        PI::AbstractParticleInstance,
        wind_tuple,
        DT, vector=true) where {PP<:Union{Dict,Nothing}}

        if defaults == nothing # this is boundary_defaults = "wind_sea"
                #@info "init particles from fetch relations: $z_i"
                particle_defaults = Dict{Num,Float64}()
                #xx, yy = PI.position_xy[1], PI.position_xy[2]
                # take in local wind velocities
                u_init, v_init = wind_tuple[1], wind_tuple[2] 
                #winds.u(xx, yy, 0), winds.v(xx, yy, 0)

                WindSeaMin = FetchRelations.get_initial_windsea(u_init, v_init, DT)
                # seed particle given fetch relations
                particle_defaults[lne] = log(WindSeaMin["E"])
                particle_defaults[c̄_x] = WindSeaMin["cg_bar_x"]
                particle_defaults[c̄_y] = WindSeaMin["cg_bar_y"]

        else
                particle_defaults = defaults
        end

        # initalize state based on state vector
        particle_defaults[x] = PI.position_xy[1]
        particle_defaults[y] = PI.position_xy[2]

        #@show defaults
        if vector
                return [particle_defaults[lne], particle_defaults[c̄_x], particle_defaults[c̄_y],  particle_defaults[x],  particle_defaults[y]]
        else
                return particle_defaults
        end

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
        State::SharedArray,
        ij::Tuple{Int, Int},

        particle_system::ODESystem,
        particle_defaults::PP,
        ODE_settings, #particle_waves_v3.ODESettings type

        GridNotes, # ad type of grid note
        winds,     # interp winds
        DT::Float64,

        boundary::Vector{T},
        periodic_boundary::Bool) where {T<:Union{Int,Any,Nothing,Int64},PP<:Union{Dict,Nothing}}

        
        xx, yy = GridNotes.x[ij[1]], GridNotes.y[ij[2]]
        uv = winds.u(xx, yy, 0)::Float64, winds.v(xx, yy, 0)::Float64

        # define initial condition
        z_i = InitParticleState(particle_defaults, ij, (xx,yy), uv, DT)
        # check if point is boundary point
        boundary_point = check_boundary_point(ij, boundary, periodic_boundary)
        #@info "boundary?", boundary_point

        # @info z_i
        # @info ij
        # @info GetParticleEnergyMomentum(z_i)
        #@info State 
        # add initial state to State vector
        init_z0_to_State!(State, ij, GetParticleEnergyMomentum(z_i))

        # Push Inital condition to collection
        push!(ParticleCollection,
        InitParticleInstance(
                particle_system,
                z_i,
                ODE_settings,
                ij,
                boundary_point))
        nothing
end


# end of module
end