module TimeSteppers

export time_step!
using Architectures
using ..mapping_1D
using ..mapping_2D

# for debugging
using Statistics

using Oceananigans.TimeSteppers: tick!

function mean_of_state(model::Abstract2DModel)
    return mean(model.State[:, :, 1])
end

function mean_of_state(model::Abstract1DModel)
    return mean(model.State[:, 1])
end


################# 1D ####################


"""
time_step!(model, Δt; callbacks=nothing)

advances model by 1 time step:
1st) the model.ParticleCollection is advanced and then 
2nd) the model.State is updated.
clock is ticked by Δt

callbacks are not implimented yet

"""
function time_step!(model::Abstract1DModel, Δt; callbacks=nothing, debug=false)

    # temporary FailedCollection to store failed particles
    FailedCollection = Vector{AbstractMarkedParticleInstance}([])

    for a_particle in model.ParticleCollection
            #@show a_particle.position_ij
            mapping_1D.advance!(    a_particle, model.State, FailedCollection, 
                                    model.grid, model.winds , Δt , 
                                    model.ODEsettings.log_energy_maximum, 
                                    model.ODEsettings.wind_min_squared,
                                    model.periodic_boundary,
                                    model.ODEdefaults)
    end
    if debug
            model.FailedCollection = FailedCollection
            @info "advanced: "
            #@info model.State[8:12, 1], model.State[8:12, 2]
            @info model.clock.time, model.ParticleCollection[10].ODEIntegrator.t
            @info model.winds(model.ParticleCollection[10].ODEIntegrator.u[3], model.ParticleCollection[10].ODEIntegrator.t)

    end

    #@printf "re-mesh"
    for a_particle in model.ParticleCollection
            mapping_1D.remesh!(     a_particle, model.State, 
                                    model.winds, model.clock.time, 
                                    model.ODEsettings, Δt,
                                    model.minimal_particle,
                                    model.minimal_state,
                                    model.ODEdefaults)
    end

    if debug
            @info "remeshed: "
            #@info model.State[8:12, 1], model.State[8:12, 2]
            @info model.clock.time, model.ParticleCollection[10].ODEIntegrator.t

    end

    tick!(model.clock, Δt)
end



################# 2D ####################

"""
time_step!(model, Δt; callbacks=nothing)

advances model by 1 time step:
1st) the model.ParticleCollection is advanced and then 
2nd) the model.State is updated.
clock is ticked by Δt

callbacks are not implimented yet

"""
function time_step!(model::Abstract2DModel, Δt::Float64; callbacks=nothing, debug=false)

    # temporary FailedCollection to store failed particles
    FailedCollection = Vector{AbstractMarkedParticleInstance}([])

    #print("mean energy before advance ", mean_of_state(model), "\n")
    if debug
        @info "before advance"
        @info maximum(model.State[:, :, 1]), maximum(model.State[:, :, 2]), maximum(model.State[:, :, 3])
        model.FailedCollection = FailedCollection
    end 

    for a_particle in model.ParticleCollection
        #@info a_particle.position_ij
        mapping_2D.advance!(    a_particle, model.State, FailedCollection,
                                model.grid, model.winds, Δt,
                                model.ODEsettings.log_energy_maximum, 
                                model.ODEsettings.wind_min_squared,
                                model.periodic_boundary,
                                model.ODEdefaults)
    end
    
    print("mean energy after advance ", mean_of_state(model), "\n")

    if debug
        @info "advanced: "
        @info maximum(model.State[:, :, 1]), maximum(model.State[:, :, 2]), maximum(model.State[:, :, 3])
        #@info model.State[8:12, 1], model.State[8:12, 2]
        @info model.clock.time, model.ParticleCollection[10].ODEIntegrator.t
        @info "winds:", model.winds.u(model.ParticleCollection[10].ODEIntegrator.u[4], model.ParticleCollection[10].ODEIntegrator.u[5], model.ParticleCollection[10].ODEIntegrator.t)
    end

    #@printf "re-mesh"
    for a_particle in model.ParticleCollection
        mapping_2D.remesh!(a_particle, model.State, 
                        model.winds, model.clock.time, 
                        model.ODEsettings, Δt,
                        model.minimal_particle, 
                        model.minimal_state,
                        model.ODEdefaults)
    end

    if debug
        @info "remeshed: "
        #@info model.State[8:12, 1], model.State[8:12, 2]
        @info maximum(model.State[:, :, 1]), maximum(model.State[:, :, 2]), maximum(model.State[:, :, 3])
        @info model.clock.time, model.ParticleCollection[10].ODEIntegrator.t

    end
    print("mean energy after remesh ", mean_of_state(model), "\n")

    tick!(model.clock, Δt)
end


"""
movie_time_step!(model, Δt; callbacks=nothing)

advances model by 1 time step:
1st) the model.ParticleCollection is advanced and then 
2nd) the model.State is updated.
clock is ticked by Δt

callbacks are not implimented yet

"""
function movie_time_step!(model::Abstract2DModel, Δt; callbacks=nothing, debug=false)

    # temporary FailedCollection to store failed particles
    FailedCollection = Vector{AbstractMarkedParticleInstance}([])

    for a_particle in model.ParticleCollection
        #@show a_particle.position_ij
        mapping_2D.advance!(a_particle, model.State, FailedCollection,
            model.grid, model.winds, Δt,
            model.ODEsettings.log_energy_maximum,
            model.ODEsettings.wind_min_squared,
            model.periodic_boundary,
            model.ODEdefaults)
    end

    model.MovieState = copy(model.State)

    if debug
        model.FailedCollection = FailedCollection
    end

    #@printf "re-mesh"
    for a_particle in model.ParticleCollection
        mapping_2D.remesh!(a_particle, model.State,
            model.winds, model.clock.time,
            model.ODEsettings, Δt,
            model.minimal_particle,
            model.minimal_state,
            model.ODEdefaults)
    end

    
    model.State .= 0.0
    tick!(model.clock, Δt)
end


end
