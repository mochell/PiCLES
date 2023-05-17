module Simulations

export Simulation
export StateStore, CashStore, state_store, init_state_store!, push_state_to_storage!, AbstractStore, convert_store_to_tuple
export add_winds_forcing_to_store!, reset_state_store!, show_stored_data, close_store!, EmptyStore, StateStore
export run!, reset_simulation!, initialize_simulation!
export convert_store_to_tuple

using Architectures: AbstractStore
include("simulation.jl")
include("storing.jl")
include("run.jl")

#push!(LOAD_PATH,   joinpath(pwd(), "code/")       )
using ParticleMesh: OneDGrid, OneDGridNotes
using PiCLES.Operators.core_1D: SeedParticle!
using PiCLES.Operators.TimeSteppers


end