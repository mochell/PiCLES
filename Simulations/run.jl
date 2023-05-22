
using ModelingToolkit: Num
using PiCLES.Operators.core_1D: ParticleDefaults

using PiCLES.Operators.core_1D: SeedParticle! as SeedParticle1D!
using PiCLES.Operators.core_2D: SeedParticle! as SeedParticle2D!

using Architectures: Abstract2DModel, Abstract1DModel
using ParticleMesh: OneDGrid, OneDGridNotes, TwoDGrid, TwoDGridNotes

#using WaveGrowthModels: init_particles!
#using WaveGrowthModels2D: init_particles!
using PiCLES.Operators.TimeSteppers

using PiCLES.Operators.mapping_1D
using PiCLES.Operators.mapping_2D
using Statistics

function mean_of_state(model::Abstract2DModel)
        return mean(model.State[:, :, 1])
end

function mean_of_state(model::Abstract1DModel)
        return mean(model.State[:, 1])
end

"""
run!(sim::Simulation; store = false, pickup=false)
main method to run the Simulation sim.
Needs time_step! to be defined for the model, and push_state_to_storage! to be defined for the store.
"""
function run!(sim; store=false, pickup=false, cash_store=false, debug=false)

        start_time_step = time_ns()

        if !(sim.initialized) # execute initialization step
                initialize_simulation!(sim)
        end

        #sim.running = true
        sim.run_wall_time = 0.0

        if sim.stop_time >= sim.model.clock.time
                sim.running = true
        else
                sim.running = false
                @info "stop_time exceeded, run not executed"
        end

        if cash_store
                sim.store = CashStore([], 1)
                sim.store.iteration += 1
                push!(sim.store.store, copy(sim.model.State))
                if sim.verbose
                        @info "write inital state to cash store..."
                end
        end

        if store
                push_state_to_storage!(sim)
                sim.store.iteration += 1
                if sim.verbose
                        @info "write inital state to store..."
                end
        end


        while sim.running

                #reset State
                sim.model.State .= 0.0
                # do time step
                
                time_step!(sim.model, sim.Δt, debug=debug)
                
                if debug & (length(sim.model.FailedCollection) > 0)
                        @info "debug mode:"
                        @info "found failed particles"
                        @info "failed particles: ", length(sim.model.FailedCollection)
                        @info "break"
                        #sim.running = false
                        # break while loop
                        break
                end

                if store
                        push_state_to_storage!(sim)
                        sim.store.iteration += 1
                        if sim.verbose
                                @info "write state to store..."
                                #print("mean energy", mean(sim.store.store["data"][:, :, 1], dims=2), "\n")
                        end

                end

                if cash_store
                        push!(sim.store.store, copy(sim.model.State))
                        sim.store.iteration += 1
                        if sim.verbose
                                @info "write state to cash store..."
                                print("mean energy ", mean_of_state(sim.model), "\n")
                        end

                end
                sim.running = sim.stop_time >= sim.model.clock.time ? true : false
                @info sim.model.clock
        end

        end_time_step = time_ns()

        # Increment the wall clock
        sim.run_wall_time += 1e-9 * (end_time_step - start_time_step)

end


"""
initialize_simulation!(sim::Simulation; particle_initials::T=nothing )
initialize the simulation sim by calling init_particles! to initialize the model.ParticleCollection.
"""
function initialize_simulation!(sim::Simulation; particle_initials::T=nothing) where {T<:Union{Dict{Num,Float64},Nothing}}
        # copy(ParticleDefaults(log(4e-8), 1e-2, 0.0)))

        if sim.verbose
                @info "init particles..."
        end
        init_particles!(sim.model, defaults=particle_initials, verbose=sim.verbose)
        
        if sim.model.clock.iteration != 0
                sim.model.clock.iteration = 0
                sim.model.clock.time = 0
        end
        
        sim.initialized = true

        nothing
end


"""
reset_simulation!(sim::Simulation; particle_initials::Dict{Num, Float64} = copy(ParticleDefaults(log( 4e-8 ), 1e-2, 0.0)) )
reset the simulation sim by calling init_particles! to reinitialize the model.ParticleCollection, sets the model.clock.time, model.clock.iteration, and model.state to 0.
"""
function reset_simulation!(sim::Simulation; particle_initials::T=nothing) where {T<:Union{Dict,Nothing}}

        sim.running = false
        sim.run_wall_time = 0.0

        sim.model.clock.iteration = 0
        sim.model.clock.time = 0

        # particles
        if sim.verbose
                @info "reset time..."
                @info "re-init particles..."
        end
        init_particles!(sim.model, defaults=particle_initials, verbose=sim.verbose)

        # state
        if sim.verbose
                @info "clear state..."
        end
        sim.model.State .= 0

        sim.initialized = true

        if sim.store isa StateStore
                reset_state_store!(sim)
        end
        nothing
end


"""
SeedParticle_mapper(f, p, s, b1, b2, b3, c1, c2, c3, c4, d1, d2 ) = x -> f( p, s, x, b1, b2, b3, c1, c2, c3, c4, d1, d2 )
maps to SeedParticle! function
"""
SeedParticle_mapper(f, p, s, b1, b2, b3, c1, c2, c3, d1, d2) = x -> f(p, s, x, b1, b2, b3, c1, c2, c3, d1, d2)


"""
init_particle!(model ; defaults::Dict{Num, Float64} = nothing, verbose::Bool=false )

initialize the model.ParticleCollection based on the model.grid and the defaults. 
If defaults is nothing, then the model.ODEdev is used.
usually the initilization uses wind constitions to seed the particles.
"""
function init_particles!(model::Abstract2DModel; defaults::T=nothing, verbose::Bool=false) where {T<:Union{Dict,Nothing}}
        #defaults        = isnothing(defaults) ? model.ODEdev : defaults
        if verbose
                @info "seed PiCLES ... \n"
                @info "defaults is $(defaults)"
                if defaults isa Dict
                        @info "found particle initials, just replace position "
                else
                        @info "no particle defaults found, use windsea to seed particles"
                end
        end

        gridnotes = TwoDGridNotes(model.grid)



        ParticleCollection = []
        SeedParticle_i = SeedParticle_mapper(SeedParticle2D!,
                ParticleCollection, model.State,
                model.ODEsystem, defaults, model.ODEsettings,
                gridnotes, model.winds, model.ODEsettings.timestep,
                model.boundary, model.periodic_boundary)

        map(SeedParticle_i, [(i, j) for i in 1:model.grid.Nx, j in 1:model.grid.Nx])



        # print(defaults)
        # ParticleCollection=[]
        # for i in range(1,length = grid1d.Nx)
        #         SeedParticle!(ParticleCollection, model.State, i,
        #                         model.ODEsystem, defaults , model.ODEsettings,
        #                         gridnotes, model.winds, model.ODEsettings.timestep, model.grid.Nx,
        #                         model.boundary, model.periodic_boundary  )
        # end

        model.ParticleCollection = ParticleCollection
        nothing
end





### 1D version ###
# """
# SeedParticle_mapper(f, p, s, b1, b2, b3, c1, c2, c3, c4, d1, d2 ) = x -> f( p, s, x, b1, b2, b3, c1, c2, c3, c4, d1, d2 )
# maps to SeedParticle! function
# """
# SeedParticle_mapper(f, p, s, b1, b2, b3, c1, c2, c3, d1, d2 )  = x -> f( p, s, x, b1, b2, b3, c1, c2, c3, d1, d2 )


"""
init_particle!(model ; defaults::Dict{Num, Float64} = nothing, verbose::Bool=false )

initialize the model.ParticleCollection based on the model.grid and the defaults. 
If defaults is nothing, then the model.ODEdev is used.
usually the initilization uses wind constitions to seed the particles.
"""
function init_particles!(model::Abstract1DModel; defaults::T=nothing, verbose::Bool=false) where {T<:Union{Dict,Nothing}}
        #defaults        = isnothing(defaults) ? model.ODEdev : defaults
        if verbose
                @info "seed PiCLES ... \n"
                @info "defaults is $(defaults)"
                if defaults isa Dict
                        @info "found particle initials, just replace position "
                else
                        @info "no particle defaults found, use windsea to seed particles"
                end
        end

        gridnotes = OneDGridNotes(model.grid)

        ParticleCollection = []
        SeedParticle_i = SeedParticle_mapper(SeedParticle1D!, ParticleCollection, model.State,
                model.ODEsystem, defaults, model.ODEsettings,
                gridnotes, model.winds, model.ODEsettings.timestep,
                model.boundary, model.periodic_boundary)

        map(SeedParticle_i, range(1, length=model.grid.Nx))


        # print(defaults)
        # ParticleCollection=[]
        # for i in range(1,length = grid1d.Nx)
        #         SeedParticle!(ParticleCollection, model.State, i,
        #                         model.ODEsystem, defaults , model.ODEsettings,
        #                         gridnotes, model.winds, model.ODEsettings.timestep, model.grid.Nx,
        #                         model.boundary, model.periodic_boundary  )
        # end

        model.ParticleCollection = ParticleCollection
        nothing
end


