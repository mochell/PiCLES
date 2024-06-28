include("dependencies_for_runtests.jl")

using PiCLES.ParticleSystems: particle_waves_v5 as PW


@testset "Particle Equations components" begin
    @testset "speed function tests" begin
        @test PW.speed(3, 4) ≈ 5.0
        @test PW.speed(-3, -4) ≈ 5.0
        @test PW.speed(-3, 4) === 5.0
        @test PW.speed(0, 0) === 0.0
    end

    @testset "speed_and_angles function tests" begin
        c, c_x, c_y = PW.speed_and_angles(3, 4)
        @test c ≈ 5.0
        @test c_x ≈ 0.6
        @test c_y ≈ 0.8
    end

    @testset "αₚ function tests" begin
        @test PW.αₚ(1, 0, 1, 0) == 0.5
        @test PW.αₚ(-1.0, 0, 1, 0) == -0.5
        @test PW.αₚ(1, 1, 1, 1) ≈ 0.5 atol = 1e-4
        @test PW.αₚ(0.0, 0.0, 1, 1) == 0.0
        @test PW.αₚ(1, 1, 0, 0) == 0.0
        @test PW.αₚ(-1, 1, 0.0, 0.1) ≈ 5.0 atol = 1e-4
        @test PW.αₚ(1.0, 1.0, 0.01, 0.0) == 50.0
        @test PW.αₚ(-1, 1, 0.0, 1e-4) == 5000.0
        @test PW.αₚ(-1, 1, 1e-6, 0) == -50.0
    end

    @testset "α_func function tests" begin
        @test PW.α_func(10, 5) === 1.0
        @test PW.α_func(200.0, 0.01) === 500.0
    end

    @testset "sin2_a_min_b function tests" begin
        @test PW.sin2_a_min_b(1, 0, 1, 0) === 0.0
        @test PW.sin2_a_min_b(1, -1, 1, -1) === 0.0
        @test PW.sin2_a_min_b(0.5, 0.5, 0.0, 1.0) ≈ 1.0 atol = 1e-4
        @test PW.sin2_a_min_b(0.5, 0.5, 1.0, 0.0) ≈ -1.0 atol = 1e-4
    end

    @testset "e_T_func function tests" begin
        q = -1 / 4.0
        p, q, n = PW.magic_fractions(q)
        @test (p, q, n) === (0.75, -0.25, 2.0)
        @test PW.e_T_func(0.88, p, q, n) ≈ 1.861 atol = 1e-3
    end

    @testset "H_β and Δ_β function tests" begin
        @test PW.H_β(0.85, 0.75) == 0.5
        @test PW.H_β(10.0, 1.0) ≈ 1 atol = 1e-4
        @test PW.H_β(-10.0, 1.0) ≈ 0 atol = 1e-4

        @test PW.Δ_β(0.85) == -0.25
        @test PW.Δ_β(0.0) ≈ 1 atol = 1e-6
        @test PW.Δ_β(10.0) == 1.0
        @test PW.Δ_β(0.9) ≈ 0.01694 atol = 1e-6
    end

    @testset "c_g_conversions function tests" begin
        c_gp, kₚ, ωₚ = PW.c_g_conversions_vector(5, r_g=0.85)
        @test c_gp ≈ 5.8823 atol = 1e-4
        @test kₚ ≈ 0.0708 atol = 1e-4
        @test ωₚ ≈ 0.83385 atol = 1e-4
        # Specific values should be tested based on the expected output of the function
    end

    @testset "S_cg and S_dir tests" begin
        C_alpha = -1.41
        @test PW.S_cg(log(4), 1, 0.0708, C_alpha) ≈ -0.00056 atol = 1e-4

        p = 0.75
        C_varphi = 1.81e-5
        u, v, c_gp_x, c_gp_y = 1.0, 0.0, 1.0, 0.0
        Hₚ = PW.H_β(PW.αₚ(u, v, c_gp_x, c_gp_y), p)
        @test PW.S_dir(u, v, c_gp_x, c_gp_y, C_varphi, Hₚ) == 0.0

        u, v, c_gp_x, c_gp_y = 10, 1, 0.0, 0.5
        Hₚ = PW.H_β(PW.αₚ(u, v, c_gp_x, c_gp_y), p)
        @test PW.S_dir(u, v, c_gp_x, c_gp_y, C_varphi, Hₚ) ≈ 0.00020 atol = 1e-5
    end


end


## test Particle Equations
#1D
#2d
#end