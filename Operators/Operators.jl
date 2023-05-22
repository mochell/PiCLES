module Operators

export core_1D, core_2D, custom_structures, mapping_1D, mapping_2D, TimeSteppers
export init_z0_to_State!
using SharedArrays


#include("custom_structures.jl")
#using .custom_structures

using Architectures: AbstractParticleInstance, AbstractMarkedParticleInstance

#include("custom_structures.jl")
using custom_structures: ParticleInstance1D, ParticleInstance2D, MarkedParticleInstance

#include("../Utils/FetchRelations.jl")
using FetchRelations
include("../ParticleSystems/ParticleSystems.jl")
using .ParticleSystems

include("utils.jl")

include("initialize.jl")
include("core_1D.jl")
include("core_2D.jl")

include("mapping_1D.jl")
include("mapping_2D.jl")

include("TimeSteppers.jl")

#using PiCLES.Utils.FetchRelations

using .core_1D
using .core_2D
using .mapping_1D
using .mapping_2D
using .TimeSteppers

end