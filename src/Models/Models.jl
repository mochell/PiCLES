module Models

export WaveGrowthModels1D, WaveGrowthModels2D, reset_boundary!

include("WaveGrowthModels1D.jl")
include("WaveGrowthModels2D.jl")

using .WaveGrowthModels1D
using .WaveGrowthModels2D

end