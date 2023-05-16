

#### initialize #####


""" sets node state values to S at tuple position ij to zi """
function init_z0_to_State!(S::SharedMatrix, ij::Int, zi::Vector{Float64} )
        S[ ij[1], :] = zi
        nothing
end

""" sets node state values to S at tuple position ij to zi """
function init_z0_to_State!(S::SharedArray, ij::Tuple{Int64,Int64}, zi::Vector{Float64})
        S[ij[1],ij[2], :] = zi
        nothing
end


""" sets particle state values to S. position is taking from particle """
function set_u_to_shared!(S::SharedMatrix, PI::ParticleInstance2D)
        S[ PI.position_ij[1], PI.position_ij[2] , :] = PI.ODEIntegrator.u
        nothing
end

""" takes node state values and pushes them to particle. position is taken from particle """
function set_u_to_particle!(S::SharedMatrix, PI::AbstractParticleInstance)
        set_u!(PI.ODEIntegrator, Get_u_FromShared(S, PI))
        u_modified!(PI.ODEIntegrator,true)
        nothing
end

