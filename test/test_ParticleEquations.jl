include("dependencies_for_runtests.jl")

using PiCLES.ParticleSystems: particle_waves_v5 as PW
using PiCLES.FetchRelations

using Oceananigans.Units

#particle_system = PW.particle_equations(u_func)
@testset "Particle Equations" begin
    ## test Particle Equations
    # 1D
    @testset "1D" begin
        # Test case 1
        @testset "Test case 1" begin
                u_func(x, t) = 10.0 + x * 0 +  t * 0 

                WindSeamin = PiCLES.FetchRelations.MinimalWindsea(u_func(0, 0), 10minutes)
                ODEpars, _, _ = PW.ODEParameters()


                particle_system = PW.particle_equations(u_func)
                zi = copy([log(WindSeamin["E"]), WindSeamin["cg_bar"], 0.0])
                dz = copy([0.0, 0.0, 0.0])
                for i in 1:20
                    dz = particle_system(dz, zi, ODEpars, 0.0)
                    zi += dz
                end
                @testset "Assertions" begin
                    @test length(zi) == 3
                    @test length(dz) == 3
                    @test all(isfinite, zi)
                    @test all(isfinite, dz)
                end

                @testset "values" begin
                    @test zi[1] ≈ -80.0 atol = 1
                    @test zi[2] ≈ 0.5 atol = 0.1
                    @test zi[3] ≈ 11.0 atol = 1
                end
        end

        @testset "Test case 2" begin
            u_func(x, t) = -10.0 + x * 0 + t * 0
            # v_func(x, t) = 5.0 + x * 0 +  t * 0

            WindSeamin = PiCLES.FetchRelations.MinimalWindsea(u_func(0, 0), 10minutes)
            ODEpars, _, _ = PW.ODEParameters()

            particle_system = PW.particle_equations(u_func)
            zi = copy([log(WindSeamin["E"]), WindSeamin["cg_bar"], 0.0])
            dz = copy([0.0, 0.0, 0.0])
            for i in 1:20
                dz = particle_system(dz, zi, ODEpars, 0.0)
                zi += dz
            end
            @testset "Assertions" begin
                @test length(zi) == 3
                @test length(dz) == 3
                @test all(isfinite, zi)
                @test all(isfinite, dz)
            end

            @testset "values" begin
                @test zi[1] ≈ -80.0 atol = 1
                @test zi[2] ≈ -0.5 atol = 0.1
                @test zi[3] ≈ -11.0 atol = 1
            end
        end

        @testset "Test case 3" begin
            T = 20.0
            u_func(x, t) = -10.0 + x * 0 + cos(pi * t/T) 

            WindSeamin = PiCLES.FetchRelations.MinimalWindsea(u_func(0, 0), 5minutes)
            ODEpars, _, _ = PW.ODEParameters()

            particle_system = PW.particle_equations(u_func)
            zi = copy([log(WindSeamin["E"]), WindSeamin["cg_bar"], 0.0])
            dz = copy([0.0, 0.0, 0.0])
            for t in 1:T:20
                dz = particle_system(dz, zi, ODEpars, t)
                zi += dz
            end

            @testset "Assertions" begin
                @test length(zi) == 3
                @test length(dz) == 3
                @test all(isfinite, zi)
                @test all(isfinite, dz)
            end
        end
    end

    @testset "2D" begin
        # Test case 1
        @testset "Test case 1" begin
            u_func(x, y, t)  = 10.0 + x * 0 + y * 0 + t * 0
            v_func(x, y, t)  = 5.0 + x * 0 + y * 0 + t * 0

            WindSeamin       = PiCLES.FetchRelations.MinimalWindsea(u_func(0, 0, 0), v_func(0,0,0) , 10minutes)
            ODEpars, _, _    = PW.ODEParameters()
            particle_system  = PW.particle_equations(u_func, v_func)
            zi               = copy([log(WindSeamin["E"]), WindSeamin["cg_bar_x"], WindSeamin["cg_bar_y"], 0.0, 0.0])
            dz               = copy([0.0, 0.0, 0.0, 0.0, 0.0])

            for i in 1:100
                dz = particle_system(dz, zi, ODEpars, 0.0)
                zi += dz
            end

            @testset "Assertions" begin
                @test length(zi) == 5
                @test length(dz) == 5
                @test all(isfinite, zi)
                @test all(isfinite, dz)
            end

            @testset "values" begin
                @test zi[1] ≈ -50.0 atol = 5
                @test zi[2] ≈ 0.5 atol = 0.2
                @test zi[3] ≈ 0.2 atol =0.2
                @test zi[4] ≈ 50.0 atol = 4
                @test zi[5] ≈ 25 atol = 4
            end
        end
    end
end
