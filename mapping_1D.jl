module mapping_1D

# includet("ParticleInCell.jl")
# includet("FetchRelations.jl")
using SharedArrays
using ParticleMesh: OneDGrid, OneDGridNotes
using Printf

import ParticleInCell
import FetchRelations

using core_1D: ParticleInstance, MarkedParticleInstance
using core_1D: GetParticleEnergyMomentum, GetVariablesAtVertex, Get_u_FromShared

using ModelingToolkit, DifferentialEquations
using particle_waves_v3: ODESettings
###### remeshing routines ############

"""
        ParticleToNode!(PI::ParticleInstance, S::SharedMatrix, G::TwoDGrid)
Pushes particle values to the neighboring nodes following the ParticleInCell rules.
1.) get weights and indexes of the neighboring notes,
2.) convert the particle state to nodestate
3.) push the calculated state to the shared arrays

inputs:

PI      Particle instance
S       Shared array where particles are stored
G       (TwoDGrid) Grid that defines the nodepositions
"""
function ParticleToNode!(PI::ParticleInstance, S::SharedMatrix, G::OneDGrid,  periodic_boundary :: Bool)


        index_positions, weights = ParticleInCell.compute_weights_and_index_2d(G, PI.ODEIntegrator.u[1])
        #ui[1:2] .= PI.position_xy
        #@show index_positions
        u_state = GetParticleEnergyMomentum(PI.ODEIntegrator.u)
        #@show u_state, index_positions, weights
        ParticleInCell.push_to_2d_grid!(S, u_state , index_positions,  weights, G.Nx ,  periodic_boundary)
        nothing
end


"""
        NodeToParticle!(PI::ParticleInstance, S::SharedMatrix)
Pushes node value to particle:
- If Node value is smaller than a minimal value, the particle is renintialized
- If Node value is okey, it is converted to state variable and pushed to particle.
- The particle position is set to the node positions
"""
function NodeToParticle!(PI::ParticleInstance, S::SharedMatrix, ti::Number, u_init::Number, e_min_log::Number, DT::Float64)
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

                #u_init       = u(PI.position_xy[1], ti)
                # derive local wave speed and direction
                cg_local     = FetchRelations.c_g_U_tau( abs(u_init) , DT )
                lne_local    = log(FetchRelations.Eⱼ( abs(u_init) , DT ))

                ui= [ PI.position_xy[1], cg_local, lne_local ]
                reinit!(PI.ODEIntegrator, ui , erase_sol=false, reset_dt=true, reinit_cache=true)#, reinit_callbacks=true)
                #set_t!(PI.ODEIntegrator, last_t )
                u_modified!(PI.ODEIntegrator,true)

        else    # load node value to particle

                #@show "get vertex variable"
                #@show "u_state", u_state
                ui= GetVariablesAtVertex( u_state, PI.position_xy[1] )
                #@show ui
                set_u!(PI.ODEIntegrator, ui )
                #set_t!(PI.ODEIntegrator, last_t )
                u_modified!(PI.ODEIntegrator,true)
        end
        #@show ui
        nothing

end


######### Core routines for advancing and remeshing

"""
        advance!(PI::ParticleInstance, S::SharedMatrix{Float64}, G::OneDGrid, DT::Float64)
"""
function advance!(      PI::ParticleInstance,
                        S::SharedMatrix{Float64},
                        Failed::Vector{MarkedParticleInstance},
                        G::OneDGrid,
                        u,
                        DT::Float64, ODEs:: ODESettings,  periodic_boundary :: Bool)
        #@show PI.position_ij

        t_start  = copy(PI.ODEIntegrator.t)
        add_saveat!(PI.ODEIntegrator, PI.ODEIntegrator.t )
        savevalues!(PI.ODEIntegrator)

        step!(PI.ODEIntegrator, DT , true)

        if check_error(PI.ODEIntegrator) != :Success
                @printf "no Success on advancing ODE, adjust reltol and retry:\n"

                @printf "push to failed\n"
                push!(Failed,
                        MarkedParticleInstance(
                                        copy(PI),
                                        copy(PI.ODEIntegrator.t),
                                        copy(PI.ODEIntegrator.u),
                                        PI.ODEIntegrator.sol.retcode
                                                ))

                # reinitialize from grid

                u_start = u(PI.position_xy[1], t_start)

                @show " time after fail",  PI.ODEIntegrator.t
                NodeToParticle!(PI, S, t_start, u_start, ODEs.log_energy_minimum, DT)
                reinit!(PI.ODEIntegrator, PI.ODEIntegrator.u , erase_sol=false, reset_dt=true, reinit_cache=true)

                #set_t!(PI.ODEIntegrator, time)
                # # change relative tolerance
                # u_modified!(PI.ODEIntegrator,true)
                # auto_dt_reset!(PI.ODEIntegrator)
                # # u_modified!(PI.ODEIntegrator,true)
                #
                PI.ODEIntegrator.opts.reltol = PI.ODEIntegrator.opts.reltol/10

                # try integration again
                step!(PI.ODEIntegrator, DT , true)
                #@show check_error(PI.ODEIntegrator)
                # readjust error
                PI.ODEIntegrator.opts.reltol = PI.ODEIntegrator.opts.reltol * 10
                #@show PI.ODEIntegrator.opts.reltol
                #@show "  final time ", PI.ODEIntegrator.t
                @printf "--- \n"
        end


        # use this for periodic boundary conditions
        #periodic_BD_single_PI!(PI, Lx) #### !!!! not sure if the boundary condition has to be set there, it might have beeen set somewhere else as well ..

        # step!(PI.ODEIntegrator, DT/2 , true)
        # PI = periodic_BD_single_PI!(PI )

        if isnan(PI.ODEIntegrator.u[1]) | isnan(PI.ODEIntegrator.u[2]) | isnan(PI.ODEIntegrator.u[3])
                @show "position or Energy is nan"
                @show PI
                set_u!(PI.ODEIntegrator, [0,0,0] )
                u_modified!(PI.ODEIntegrator,true)

        elseif isinf(PI.ODEIntegrator.u[1]) | isinf(PI.ODEIntegrator.u[2]) | isinf(PI.ODEIntegrator.u[3])
                @show "position or Energy is inf"
                @show PI
                set_u!(PI.ODEIntegrator, [0,0,0] )
                u_modified!(PI.ODEIntegrator,true)

        elseif PI.ODEIntegrator.u[3] > ODEs.log_energy_maximum
                @show "e_max_log is reached"
                @show PI
                set_u!(PI.ODEIntegrator, [0,0,0] )
                u_modified!(PI.ODEIntegrator,true)
        end

        ParticleToNode!(PI, S, G, periodic_boundary)
        return PI
end

"""
        remesh!(PI::ParticleInstance, S::SharedMatrix{Float64, 3})
        Wrapper function that does everything necessary to remesh the particles.
        - pushes the Node State to particle instance
"""
function remesh!(PI::ParticleInstance, S::SharedMatrix{Float64}, u, ti::Number, ODEs:: ODESettings, DT::Float64)
        ui = u(PI.position_xy[1], ti)
        NodeToParticle!(PI, S, ti, ui, ODEs.log_energy_minimum, DT)
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
