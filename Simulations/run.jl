


using ModelingToolkit: Num
using core_1D: ParticleDefaults

using WaveGrowthModels: init_particles!
using Statistics

"""
run!(sim::Simulation; store = false, pickup=false)
main method to run the Simulation sim.
Needs time_step! to be defined for the model, and push_state_to_storage! to be defined for the store.
"""
function run!(sim; store = false, pickup=false)

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

    while sim.running
            time_step!(sim.model, sim.Δt)
            if store
                    push_state_to_storage!(sim)
                    sim.store.iteration += 1
                    if sim.verbose
                            @info "write state to store..."
                            print("state state",  mean(sim.store.store["data"][:, :, 1], dims =2) , "\n")
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
initialize_simulation!(sim::Simulation; defaults::Dict{Num, Float64} = copy(ParticleDefaults(0.0, 1e-2, log( 4e-8 ))) )
initialize the simulation sim by calling init_particles! to initialize the model.ParticleCollection.
"""
function initialize_simulation!(sim::Simulation ; defaults::Dict{Num, Float64} = copy(ParticleDefaults(0.0, 1e-2, log( 4e-8 ))) )

    if sim.verbose
            @info "init particles..."
    end
    init_particles!(sim.model, defaults= defaults, verbose= sim.verbose )
    sim.initialized = true
    sim.initialized = true
    nothing
end


"""
reset_simulation!(sim::Simulation; defaults::Dict{Num, Float64} = copy(ParticleDefaults(0.0, 1e-2, log( 4e-8 ))) )
reset the simulation sim by calling init_particles! to reinitialize the model.ParticleCollection, sets the model.clock.time, model.clock.iteration, and model.state to 0.
"""
function reset_simulation!(sim::Simulation ; defaults::Dict{Num, Float64} = copy(ParticleDefaults(0.0, 1e-2, log( 4e-8 ))) )

    sim.running = false
    sim.run_wall_time = 0.0

    sim.model.clock.iteration= 0
    sim.model.clock.time= 0

    # particles
    if sim.verbose
            @info "reset time..."
            @info "re-init particles..."
    end
    init_particles!(sim.model, defaults= defaults , verbose= sim.verbose )

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



