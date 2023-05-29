module Plotting

export movie

#push!(LOAD_PATH, joinpath(pwd(), "code/Simulations"))
# include("../Simulations/Simulations.jl")
# using .Simulations

using ..Simulations: convert_store_to_tuple, CashStore, StateStore, Simulation
using ParticleMesh: TwoDGrid, TwoDGridNotes, TwoDGridMesh

include("plotting_1D.jl")
include("movie_2D.jl")

using .movie

end