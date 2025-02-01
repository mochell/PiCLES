module Plotting

export movie, PlotState_SingleGlobe, PlotState_DoubleGlobe

# include("../Simulations/Simulations.jl")
# using .Simulations

using ..Simulations: convert_store_to_tuple, CashStore, StateStore, Simulation
using ..ParticleMesh: TwoDGrid, TwoDGridNotes, TwoDGridMesh

include("plotting_1D.jl")
include("movie_2D.jl")
include("global.jl")

using .movie

end