module mapping_2D

using SharedArrays
using ParticleMesh: TwoDGrid, TwoDGridNotes
using Printf

import ParticleInCell as PIC
import .FetchRelations

using custom_structures: ParticleInstance1D, ParticleInstance2D, MarkedParticleInstance
using ..core_2D: GetParticleEnergyMomentum, GetVariablesAtVertex, Get_u_FromShared

using DifferentialEquations
using ..particle_waves_v3: ODESettings

using Architectures: AbstractParticleInstance, AbstractMarkedParticleInstance, AbstractODESettings
###### remeshing routines ############

"""
        ParticleToNode!(PI::AbstractParticleInstance, S::SharedMatrix, G::TwoDGrid)
Pushes particle values to the neighboring nodes following the ParticleInCell rules.
1.) get weights and indexes of the neighboring notes,
2.) convert the particle state to nodestate
3.) push the calculated state to the shared arrays

inputs:

PI      Particle instance
S       Shared array where particles are stored
G       (TwoDGrid) Grid that defines the nodepositions
"""

function ParticleToNode!(PI::AbstractParticleInstance, S::SharedArray, G::TwoDGrid, periodic_boundary::Bool)

        #u[4], u[5] are the x and y positions of the particle
        index_positions, weights = PIC.compute_weights_and_index(G, PI.ODEIntegrator.u[4], PI.ODEIntegrator.u[5])
        #ui[1:2] .= PI.position_xy
        #@show index_positions
        u_state = GetParticleEnergyMomentum(PI.ODEIntegrator.u)
        #@show u_state
        PIC.push_to_grid!(S, u_state , index_positions,  weights, G.Nx, G.Ny , periodic_boundary)
        nothing
end

"""
        NodeToParticle!(PI::AbstractParticleInstance, S::SharedMatrix)
Pushes node value to particle:
- If Node value is smaller than a minimal value, the particle is renintialized
- If Node value is okey, it is converted to state variable and pushed to particle.
- The particle position is set to the node positions
"""
function NodeToParticle!(PI::AbstractParticleInstance, S::SharedArray, ti::Number, wind_tuple::Tuple{Number,Number}, e_min_log::Number, DT::Float64)
        u_state = Get_u_FromShared( PI, S)


        last_t = PI.ODEIntegrator.t
        # test if particle is below energy threshold, or
        #      if particle is at the boundary
        if (u_state[1] < exp(e_min_log)) | PI.boundary # init new particle

                #@show "re-init new particle"
                # @show PI.position_ij, u_state
                # z_i = copy(z0)
                # z_i[x] = PI.position_xy[1]
                # z_i[y] = PI.position_xy[2]
                # ui = z_i

                # derive local wave speed and direction
                u_init, v_init = wind_tuple
                u_abs          = sqrt(u_init .^ 2 + v_init .^ 2) # absolute value of u


                # define new particle state
                cg_x_local     = FetchRelations.c_g_U_tau( abs(u_init) , DT )
                cg_y_local     = FetchRelations.c_g_U_tau(abs(v_init), DT)
                lne_local      = log(FetchRelations.Eⱼ(abs(u_abs), DT))

                ui             = [lne_local, cg_x_local, cg_y_local, PI.position_xy[1], PI.position_xy[2]]
                reinit!(PI.ODEIntegrator, ui , erase_sol=false, reset_dt=true, reinit_cache=true)#, reinit_callbacks=true)
                #set_t!(PI.ODEIntegrator, last_t )
                u_modified!(PI.ODEIntegrator,true)

        else    # load node value to particle

                #@show "get vertex variable"
                #@show "u_state", u_state
                ui = GetVariablesAtVertex( u_state, PI.position_xy[1], PI.position_xy[2] )
                #@show ui

                # this method is more robust than the set_u! method
                reinit!(PI.ODEIntegrator, ui, erase_sol=false, reset_dt=true, reinit_cache=true)#, reinit_callbacks=true)
                #set_u!(PI.ODEIntegrator, ui )
                #set_t!(PI.ODEIntegrator, last_t )
                u_modified!(PI.ODEIntegrator,true)
        end
        #@show ui
        nothing

end


######### Core routines for advancing and remeshing

"""
        advance!(PI::AbstractParticleInstance, S::SharedMatrix{Float64}, G::TwoDGrid, DT::Float64)
"""
function advance!(PI::AbstractParticleInstance,
                        S::SharedArray,
                        Failed::Vector{AbstractMarkedParticleInstance},
                        G::TwoDGrid,
                        winds,
                        DT::Float64, ODEs:: AbstractODESettings,  periodic_boundary :: Bool)
        #@show PI.position_ij

        t_start  = copy(PI.ODEIntegrator.t)
        add_saveat!(PI.ODEIntegrator, PI.ODEIntegrator.t )
        savevalues!(PI.ODEIntegrator)
        
        # advance particle
        step!(PI.ODEIntegrator, DT , true)

        # check if integration was successful
        if check_error(PI.ODEIntegrator) != :Success
                @printf "no Success on advancing ODE:\n"
                print("- time after fail $(PI.ODEIntegrator.t)\n ")

                @printf "- push to failed\n"
                push!(Failed,
                        MarkedParticleInstance(
                                        copy(PI),
                                        copy(PI.ODEIntegrator.t),
                                        copy(PI.ODEIntegrator.u),
                                        PI.ODEIntegrator.sol.retcode
                                                ))

                # reinitialize from grid
                @printf "- adjust reltol and retry:\n"
                winds_start = winds.u(PI.position_xy[1], PI.position_xy[2], t_start), winds.v(PI.position_xy[1], PI.position_xy[2], t_start)

                NodeToParticle!(PI, S, t_start, winds_start, ODEs.log_energy_minimum, DT)
                reinit!(PI.ODEIntegrator, PI.ODEIntegrator.u , erase_sol=false, reset_dt=true, reinit_cache=true)

                #set_t!(PI.ODEIntegrator, time)
                # # change relative tolerance
                # u_modified!(PI.ODEIntegrator,true)
                # auto_dt_reset!(PI.ODEIntegrator)
                u_modified!(PI.ODEIntegrator, true)
                #

                PI.ODEIntegrator.opts.reltol = PI.ODEIntegrator.opts.reltol/10

                # try integration again
                step!(PI.ODEIntegrator, DT , true)
                #@show check_error(PI.ODEIntegrator)
                # readjust error
                PI.ODEIntegrator.opts.reltol = PI.ODEIntegrator.opts.reltol * 10
                #@show PI.ODEIntegrator.opts.reltol
                #@show "  final time ", PI.ODEIntegrator.t
                print("- 2nd try success? ", check_error(PI.ODEIntegrator) != :Success, "\n" )
                @printf "--- \n"
        end


        # use this for periodic boundary conditions
        #periodic_BD_single_PI!(PI, Lx) #### !!!! not sure if the boundary condition has to be set there, it might have beeen set somewhere else as well ..

        # step!(PI.ODEIntegrator, DT/2 , true)
        # PI = periodic_BD_single_PI!(PI )

        if isnan(PI.ODEIntegrator.u[1]) | isnan(PI.ODEIntegrator.u[2]) | isnan(PI.ODEIntegrator.u[3])
                @show "position or Energy is nan"
                @show PI
                set_u!(PI.ODEIntegrator, [0,0,0, 0,0] )
                u_modified!(PI.ODEIntegrator,true)

        elseif isinf(PI.ODEIntegrator.u[1]) | isinf(PI.ODEIntegrator.u[2]) | isinf(PI.ODEIntegrator.u[3])
                @show "position or Energy is inf"
                @show PI
                set_u!(PI.ODEIntegrator, [0,0,0, 0, 0] )
                u_modified!(PI.ODEIntegrator,true)

        elseif PI.ODEIntegrator.u[1] > ODEs.log_energy_maximum
                @show "e_max_log is reached"
                #@show PI
                set_u!(PI.ODEIntegrator, [0,0,0, 0, 0] )
                u_modified!(PI.ODEIntegrator,true)
        end

        ParticleToNode!(PI, S, G, periodic_boundary)
        return PI
end

"""
        remesh!(PI::ParticleInstance2D, S::SharedMatrix{Float64, 3})
        Wrapper function that does everything necessary to remesh the particles.
        - pushes the Node State to particle instance
"""
function remesh!(PI::ParticleInstance2D, S::SharedArray{Float64}, winds, ti::Number, ODEs::AbstractODESettings, DT::Float64)        #ui = u(PI.position_xy[1], ti)
        winds_i = winds.u(PI.position_xy[1], PI.position_xy[2], ti), winds.v(PI.position_xy[1], PI.position_xy[2], ti)
        NodeToParticle!(PI, S, ti, winds_i, ODEs.log_energy_minimum, DT)
        return PI
end




""" shows total energy of the all particles """
function ShowTotalEnergyChange(ParticleCollection, u_sum_m1)
        u_sum = zeros(Nstate)
        for a_particle in ParticleCollection
                u_sum +=  GetParticleEnergyMomentum(a_particle.ODEIntegrator.u)
        end
        @show u_sum_m1 - u_sum
        return u_sum
end

end
