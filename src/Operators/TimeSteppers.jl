module TimeSteppers

export time_step!
using ...Architectures
using ..mapping_1D
using ..mapping_2D

# for debugging
using Statistics
using Base.Threads
using Printf

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

    @threads for a_particle in model.ParticleCollection
        #@info a_particle.position_ij
        mapping_2D.advance!(    a_particle, model.ParticlesAtNode, model.State, FailedCollection,
                                model.grid, model.winds, Δt,
                                model.ODEsettings.log_energy_maximum, 
                                model.ODEsettings.wind_min_squared,
                                model.periodic_boundary,
                                model.ODEdefaults)
    end
    
    for (i,j) in [(i,j) for i in 1:model.grid.Nx for j in 1:model.grid.Ny]
        weights = [model.ParticlesAtNode[i][j][k][1] for k in 1:length(model.ParticlesAtNode[i][j])]
        values = [model.ParticlesAtNode[i][j][k][2] for k in 1:length(model.ParticlesAtNode[i][j])]
        if length(weights) > 0
            model.State[i,j,4] = sum(weights .* values) / sum(weights)
        
            for k in 1:length(model.ParticlesAtNode[i][j])
                pop!(model.ParticlesAtNode[i][j])
            end
        end
    end    

    if debug
        print("mean energy after advance ", mean_of_state(model), "\n")

        @info "advanced: "
        @info maximum(model.State[:, :, 1]), maximum(model.State[:, :, 2]), maximum(model.State[:, :, 3])
        #@info model.State[8:12, 1], model.State[8:12, 2]
        @info model.clock.time, model.ParticleCollection[10].ODEIntegrator.t
        @info "winds:", model.winds.u(model.ParticleCollection[10].ODEIntegrator.u[4], model.ParticleCollection[10].ODEIntegrator.u[5], model.ParticleCollection[10].ODEIntegrator.t)
    end

    #@printf "re-mesh"
    if model.angular_spreading_type == "stochast"
        @threads for a_particle in model.ParticleCollection
            mapping_2D.remesh!(a_particle, model.State, 
                            model.winds, model.clock.time, 
                            model.ODEsettings, Δt,
                            model.minimal_particle, 
                            model.minimal_state,
                            model.ODEdefaults)
        end
    elseif model.angular_spreading_type == "geometrical"
        i = 1
        particlesToBeReset = []
        
        for _ in 1:length(model.ParticleCollection)
            pos_ij = model.ParticleCollection[i].position_ij
            big_enough = model.State[pos_ij[1], pos_ij[2],2]^2+model.State[pos_ij[1], pos_ij[2],3]^2>model.minimal_state[2]
            if model.ParticleCollection[i].boundary || ~big_enough
                i=i+1
            elseif !(model.ParticleCollection[i].position_ij in particlesToBeReset)
                push!(particlesToBeReset, model.ParticleCollection[i].position_ij)
                deleteat!(model.ParticleCollection, i)
            end
        end
        @threads for a_particle in model.ParticleCollection
            mapping_2D.remesh!(a_particle, model.State, 
                            model.winds, model.clock.time, 
                            model.ODEsettings, Δt,
                            model.minimal_particle, 
                            model.minimal_state,
                            model.ODEdefaults)
        end
        @threads for (i,j) in particlesToBeReset
            mapping_2D.remesh!(i, j, model, Δt)
        end
        @info particlesToBeReset
    elseif model.angular_spreading_type == "nonparametric"
        
    end

    if debug
        @info "remeshed: "
        #@info model.State[8:12, 1], model.State[8:12, 2]
        @info maximum(model.State[:, :, 1]), maximum(model.State[:, :, 2]), maximum(model.State[:, :, 3])
        @info model.clock.time, model.ParticleCollection[10].ODEIntegrator.t

    end
    #print("mean energy after remesh ", mean_of_state(model), "\n")

    tick!(model.clock, Δt)
    @printf(" mean energy %.6f ", mean_of_state(model))

end


#build wrapper
advance_wrapper(f, state, Fcol, grid, winds, dt, emax, windmin, boundary, defaults) = x -> f(x, state, Fcol, grid, winds, dt, emax, windmin, boundary, defaults)
remesh_wrapper(f, state, winds, time, sets, dt, minpar, minstate, defaults) = x -> f(x, state, winds, time, sets, dt, minpar, minstate, defaults)
#global ParticleCollection


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
