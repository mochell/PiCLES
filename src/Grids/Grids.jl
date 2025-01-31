module Grids

using ...Architectures: AbstractGrid, AbstractGridStatistics, TripolarGrid

export TripolarGridMOM6
export CartesianGrid
export SphericalGrid
#, TwoDCartesianGridStatistics, TwoDCartesianGridMesh

include("mask_utils.jl")
include("spherical_grid_corrections.jl")

include("TripolarGridMOM6.jl")
using .TripolarGridMOM6

include("CartesianGrid.jl")
using .CartesianGrid

include("SphericalGrid.jl")
using .SphericalGrid


end # module