

#### initialize #####

using StaticArrays

### 1D

""" sets node state values to S at tuple position ij to zi """
function init_z0_to_State!(S::SharedMatrix, ij::Int, zi::TT) where TT<:Union{Vector{Float64},SVector{3,Float64}}
        S[ ij[1], :] = zi
        nothing
end

### 2D

""" sets node state values to S at tuple position ij to zi """
function init_z0_to_State!(S::StateTypeL1, ij::Tuple{Int64,Int64}, zi::TT) where {TT<:Union{Vector{Float64},SVector{3,Float64}}}
        S[ij[1],ij[2], :] = zi
        nothing
end


# not used anymore
# """ sets particle state values to S. position is taking from particle """
# function set_u_to_shared!(S::StateTypeL1, PI::ParticleInstance)
#         S[ PI.position_ij[1], PI.position_ij[2] , :] = PI.ODEIntegrator.u
#         nothing
# end

""" takes node state values and pushes them to particle. position is taken from particle """
function set_u_to_particle!(S::StateTypeL1, PI::AbstractParticleInstance)
        set_u!(PI.ODEIntegrator, Get_u_FromShared(S, PI))
        u_modified!(PI.ODEIntegrator,true)
        nothing
end


### multi-layer version

""" sets node state values to S at tuple position ij to zi """
function init_z0_to_State!(S::StateTypeLN, ij::Tuple{Int64,Int64,Int64}, zi::TT) where {TT<:Union{Vector{Float64},SVector{3,Float64}}}
        S[ij[1], ij[2], ij[3], :] = zi
        nothing
end

# """ sets particle state values to S. position is taking from particle """
# function set_u_to_shared!(S::StateTypeLN, PI::ParticleInstance)
#         S[PI.position_ij[1], PI.position_ij[2], :] = PI.ODEIntegrator.u
#         nothing
# end

# """ takes node state values and pushes them to particle. position is taken from particle """
# function set_u_to_particle!(S::StateTypeLN, PI::AbstractParticleInstance)
#         set_u!(PI.ODEIntegrator, Get_u_FromShared(S, PI))
#         u_modified!(PI.ODEIntegrator, true)
#         nothing
# end
