module mapping_2D

using SharedArrays
using StaticArrays
using DifferentialEquations
using Printf

using ...ParticleMesh: TwoDGrid, TwoDGridNotes
import ...ParticleInCell as PIC

using ...FetchRelations

using ...custom_structures: ParticleInstance1D, ParticleInstance2D, MarkedParticleInstance

using ..core_2D: GetParticleEnergyMomentum, GetVariablesAtVertex, Get_u_FromShared, ResetParticleValues, ParticleDefaults


using ...Architectures: AbstractParticleInstance, AbstractMarkedParticleInstance, AbstractODESettings, StateTypeL1
###### remeshing routines ############


speed(x::Float64, y::Float64) = sqrt(x^2 + y^2)
speed_square(x::Float64, y::Float64) = x^2 + y^2

"""
        ParticleToNode!(PI::AbstractParticleInstance, S::SharedMatrix, G::TwoDGrid)
Pushes particle values to the neighboring nodes following the ParticleInCell rules.
1.) get weights and indexes of the neighboring notes,
2.) convert the particle state to nodestate
3.) push the calculated state to the shared arrays

inputs:

PI      Particle instance
S       Shared array where particles are stored
G       (TwoDGrid) Grid that defines the nodepositions
"""

function ParticleToNode!(PI::AbstractParticleInstance, S::StateTypeL1, G::TwoDGrid, periodic_boundary::Bool)

        #u[4], u[5] are the x and y positions of the particle
        #index_positions, weights = PIC.compute_weights_and_index(G, PI.ODEIntegrator.u[4], PI.ODEIntegrator.u[5])
        weights_and_index = PIC.compute_weights_and_index_mininal(G, PI.ODEIntegrator.u[4], PI.ODEIntegrator.u[5])

        #ui[1:2] .= PI.position_xy
        #@show index_positions
        u_state = GetParticleEnergyMomentum(PI.ODEIntegrator.u)
        #@show u_state

        #PIC.push_to_grid!(S, u_state , index_positions,  weights, G.Nx, G.Ny , periodic_boundary)
        PIC.push_to_grid!(S, u_state , weights_and_index, G.Nx, G.Ny , periodic_boundary)
        nothing
end

"""
        ParticleToNode!(PI::AbstractParticleInstance, S::SharedMatrix, u_state::Vector{Float64})
Pushes values u_state to the node in S of the particle origin of PI.
"""
function ParticleToNode!(PI::AbstractParticleInstance, S::SharedMatrix, u_state::Vector{Float64})

        S[PI.position_ij[1], PI.position_ij[2], :] = u_state
        nothing
end

function set_u_and_t!(integrator, u_new::CC, t_new::Number) where CC <:Union{Vector{Float64},MVector}
        integrator.u = u_new
        integrator.t = t_new
end


function reset_PI_u!(PI::AbstractParticleInstance; ui::CC) where CC<:Union{Vector{Float64},MVector{Float64}}  
        # this method keeps the correct time for time varying forcing (~may 2023)
        set_u!(PI.ODEIntegrator, ui)
        u_modified!(PI.ODEIntegrator, true)
        auto_dt_reset!(PI.ODEIntegrator)
end


function reset_PI_ut!(PI::AbstractParticleInstance; ui::CC, ti::Number) where CC <:Union{Vector{Float64},MVector}
        # this method keeps the correct time for time varying forcing (~may 2023)
        set_u_and_t!(PI.ODEIntegrator, ui, ti)
        u_modified!(PI.ODEIntegrator, true)
        auto_dt_reset!(PI.ODEIntegrator)
end

function reset_PI_t!(PI::AbstractParticleInstance; ti::Number)
        # this method keeps the correct time for time varying forcing (~may 2023)
        set_t!(PI.ODEIntegrator, ti)
        u_modified!(PI.ODEIntegrator, true)
        auto_dt_reset!(PI.ODEIntegrator)
end

######### Core routines for advancing and remeshing

"""
        advance!(PI::AbstractParticleInstance, S::SharedMatrix{Float64}, G::TwoDGrid, DT::Float64)
"""
function advance!(PI::AbstractParticleInstance,
                        S::StateTypeL1,
                        Failed::Vector{AbstractMarkedParticleInstance},
                        G::TwoDGrid,
                        winds::NamedTuple{(:u, :v)},
                        DT::Float64, 
                        log_energy_maximum::Float64,
                        wind_min_squared::Float64,
                        periodic_boundary::Bool, 
                        default_particle::PP,
                        ) where {PP<:Union{ParticleDefaults,Nothing}}
        #@show PI.position_ij

        #if ~PI.boundary # if point is not a 
        t_start  =  copy(PI.ODEIntegrator.t)
        add_saveat!(PI.ODEIntegrator, PI.ODEIntegrator.t )
        savevalues!(PI.ODEIntegrator)
        
        # advance particle
        if PI.on #& ~PI.boundary # if Particle is on and not boundary
        
                try
                        step!(PI.ODEIntegrator, DT, true)
                catch e
                        @printf "error on advancing ODE:\n"
                        print("- time after fail $(PI.ODEIntegrator.t)\n ")
                        print("- error message: $(e)\n")
                        print("- push to failed\n")
                        print("- state of particle: $(PI.ODEIntegrator.u)\n")
                        print("- winds are: $(winds.u( PI.ODEIntegrator.u[4], PI.ODEIntegrator.u[5], PI.ODEIntegrator.t))\n")
                        print("- winds are: $(winds.v( PI.ODEIntegrator.u[4], PI.ODEIntegrator.u[5], PI.ODEIntegrator.t))\n")
                        push!(Failed,
                                MarkedParticleInstance(
                                        copy(PI),
                                        copy(PI.ODEIntegrator.t),
                                        copy(PI.ODEIntegrator.u),
                                        PI.ODEIntegrator.sol.retcode
                                ))
                        return

                end
        
        elseif ~PI.on #& ~PI.boundary # particle is off, test if there was windsea

                t_end = t_start + DT
                wind_end = convert(Tuple{Float64,Float64},
                                (winds.u(PI.position_xy[1], PI.position_xy[2], t_end),
                                winds.v(PI.position_xy[1], PI.position_xy[2], t_end)))::Tuple{Float64,Float64}

                # test if winds where strong enough
                if speed_square(wind_end[1], wind_end[2]) >= wind_min_squared
                        # winds are large eneough, reinit
                        ui = ResetParticleValues(default_particle, PI, wind_end, DT)
                        reset_PI_u!(PI, ui =ui)
                        PI.on = true
                end

        else    #particle is on and boundary
                
                #@info "particle is on and boundary"
                # particle stays off or is bounaday. do not advance
                PI.on=false
                return
        end

        # # check if integration reached limits or is nan, or what ever. if so, reset
        if sum(isnan.(PI.ODEIntegrator.u[1:3])) > 0
                @info "position or Energy is nan, reset"
                @info PI.position_ij
                @show PI
                
                t_end = t_start + DT
                winds_start = convert(  Tuple{Float64,Float64},
                        (winds.u(PI.position_xy[1], PI.position_xy[2], t_end),
                        winds.v(PI.position_xy[1], PI.position_xy[2], t_end)))::Tuple{Float64,Float64}
                @show winds_start

                ui = ResetParticleValues(default_particle, PI, winds_start, DT)
                @show PI.ODEIntegrator.u
                reset_PI_u!(PI, ui=ui)

        elseif  sum(isinf.(PI.ODEIntegrator.u[1:3])) > 0
                @info "position or Energy is inf"
                @show PI

                winds_start = convert(Tuple{Float64,Float64},
                                        (winds.u(PI.position_xy[1], PI.position_xy[2], t_start),
                                        winds.v(PI.position_xy[1], PI.position_xy[2], t_start)))::Tuple{Float64,Float64}

                ui = ResetParticleValues(default_particle, PI, winds_start, DT)
                reset_PI_u!(PI, ui=ui)

        elseif PI.ODEIntegrator.u[1] > log_energy_maximum
                @info "e_max_log is reached"
                #@show PI

                winds_start = convert(Tuple{Float64,Float64},
                                        (winds.u(PI.position_xy[1], PI.position_xy[2], t_start),
                                        winds.v(PI.position_xy[1], PI.position_xy[2], t_start)))::Tuple{Float64,Float64}

                ui = ResetParticleValues(default_particle, PI, winds_start, DT)
                reset_PI_u!(PI, ui=ui)

        end

        #if PI.ODEIntegrator.u[1] > -13.0 #ODEs.log_energy_minimum # the minimum enerçy is distributed to 4 neighbouring particles
        if PI.on 
                ParticleToNode!(PI, S, G, periodic_boundary)
        end

        #return PI
end

"""
        remesh!(PI::ParticleInstance2D, S::SharedMatrix{Float64, 3})
        Wrapper function that does everything necessary to remesh the particles.
        - pushes the Node State to particle instance
"""
function remesh!(PI::ParticleInstance2D, S::StateTypeL1,
                winds::NamedTuple{(:u, :v)}, 
                ti::Number, 
                ODEs::AbstractODESettings, DT::Float64,  #
                minimal_particle::Vector{Float64}, minimal_state::Vector{Float64},
                default_particle::PP) where {PP<:Union{ParticleDefaults,Nothing}}        
                
        winds_i::Tuple{Float64,Float64} = winds.u(PI.position_xy[1], PI.position_xy[2], ti), winds.v(PI.position_xy[1], PI.position_xy[2], ti)
        
        NodeToParticle!(PI, S, 
                        winds_i, 
                        minimal_particle, 
                        minimal_state,
                        ODEs.wind_min_squared,
                        default_particle, 
                        ODEs.log_energy_minimum, 
                        DT)
        #return PI
end


"""
        NodeToParticle!(PI::AbstractParticleInstance, S::SharedMatrix)
Pushes node value to particle:
- If Node value is smaller than a minimal value, the particle is renintialized
- If Node value is okey, it is converted to state variable and pushed to particle.
- The particle position is set to the node positions
"""
function NodeToParticle!(PI::AbstractParticleInstance, S::StateTypeL1,
        wind_tuple::Tuple{Float64,Float64}, 
        minimal_particle::Vector{Float64}, 
        minimal_state::Vector{Float64},
        wind_min_squared::Float64, 
        default_particle::PP, 
        e_min_log::Number, 
        DT::Float64,) where {PP<:Union{ParticleDefaults,Nothing}}


        # load data from shared array
        u_state = Get_u_FromShared(PI, S)

        last_t = PI.ODEIntegrator.t
        # minimal_state[1] is the minmal Energy  
        # minimal_state[2] is the minmal momentum squared  
        if ~PI.boundary & (u_state[1] >= minimal_state[1]) & (speed_square(u_state[2], u_state[3]) >= minimal_state[2]) # all integrior nodes: convert note to particle values and push to ODEIntegrator

                #@show "u_state", u_state
                ui = GetVariablesAtVertex(u_state, PI.position_xy[1], PI.position_xy[2])
                #@info exp(ui[1]), ui[2], ui[4]/1e3, ui[5]/1e3              
                reset_PI_ut!(PI; ui=ui, ti=last_t)
                PI.on = true

                # this method is more robust than the set_u! method (~february 2023)
                # reinit!(PI.ODEIntegrator, ui, erase_sol=false, reset_dt=true, reinit_cache=true)#, reinit_callbacks=true)
                # set_t!(PI.ODEIntegrator, last_t )
                # u_modified!(PI.ODEIntegrator,true)
                # auto_dt_reset!(PI.ODEIntegrator)

                # # this method is more robust than the set_u! method
                # reinit!(PI.ODEIntegrator, ui, erase_sol=false, reset_dt=true, reinit_cache=true)#, reinit_callbacks=true)
                # #set_u!(PI.ODEIntegrator, ui )
                # #set_t!(PI.ODEIntegrator, last_t )
                # u_modified!(PI.ODEIntegrator,true)

                #@show PI.ODEIntegrator.t, PI.ODEIntegrator.u
                
        elseif ~PI.boundary & (speed_square(wind_tuple[1], wind_tuple[2]) >= wind_min_squared) #minimal windsea is not big enough but local winds are strong enough  #(u_state[1] < exp(e_min_log)) | PI.boundary
                # test if particle is below energy threshold, or
                #      if particle is at the boundary

                ui = ResetParticleValues(default_particle, PI, wind_tuple, DT) # returns winds sea given DT and winds
                reinit!(PI.ODEIntegrator, ui, erase_sol=false, reset_dt=true, reinit_cache=true)#, reinit_callbacks=true)
                reset_PI_t!(PI, ti=last_t)

                PI.on = true

        elseif PI.boundary & (speed_square(wind_tuple[1], wind_tuple[2]) >= wind_min_squared) # at the boundary, reset particle if winds are strong enough

                ui = ResetParticleValues(default_particle, PI, wind_tuple, DT) # returns winds sea given DT and winds
                reinit!(PI.ODEIntegrator, ui, erase_sol=false, reset_dt=true, reinit_cache=true)#, reinit_callbacks=true)
                reset_PI_t!(PI, ti=last_t)

                PI.on = true

 
        else # particle is below energy threshold & on boundary
                #PI.ODEIntegrator.u = ResetParticleValues(minimal_particle, PI, wind_tuple, DT)
                # if ~PI.boundary
                #         @info u_state
                # end
                PI.on = false
        end
        nothing

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
