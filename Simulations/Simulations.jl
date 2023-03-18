module Simulations

export Simulation
export StateStore, state_store, init_state_store!, push_state_to_storage!, AbstractStore
export add_winds_forcing_to_store!, reset_state_store!, show_stored_data, close_store!, EmptyStore, StateStore
export run!, reset_simulation!, initialize_simulation!

using WaveGrowthModels

include("storing.jl")
include("simulation.jl")
include("run.jl")

push!(LOAD_PATH,   joinpath(pwd(), "code/")       )
using ParticleMesh: OneDGrid, OneDGridNotes
using core_1D: SeedParticle!

using Oceananigans.Models: AbstractModel




end