
include("dependencies_for_runtests.jl")

@testset "PiCLES.jl" begin

    # temporal Fetch relations

    @testset "Particle Equations tests" begin
        include("test_ParticleEquations_components.jl")
        include("test_ParticleEquations.jl")
    end

    include("test_T02_single_particle_2d.jl")
    include("test_T03_PIC_propagation.jl")

    # structures

    # initializations

    # Model runs
    # include("test_2D_box.jl")


end