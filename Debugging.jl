module Debugging



###### Debugging functions ############

"""
        ResetParticle!(integrator)
(debubugging function)
Resets the integrator instance if the particle energy is nan or very high
resets the particles position to the domain center
resets the energy to a dfault value (e_min_log is a global variable)
"""
function ResetParticle!(integrator)
        if isnan(integrator.ODEIntegrator.u[3]) || exp(integrator.ODEIntegrator.u[3]) >= 1e-3
                @show exp(integrator.ODEIntegrator.u[3])
                integrator.ODEIntegrator.u[1] = integrator.ODEIntegrator.u[1] - Lx/2
                #integrator.ODEIntegrator.u[2] = integrator.ODEIntegrator.u[2] - Ly/2
                #integrator.ODEIntegrator.u[4] = 1e-2
                integrator.ODEIntegrator.u[3] = e_min_log
                u_modified!(integrator.ODEIntegrator,true)
                @show "rest particle"
        end
        nothing
end

Lx_terminate_limit = 1


"""
        TerminateCheckSingle!(integrator)
(debubugging function)

"""
function TerminateCheckSingle!(integrator)
        if maximum(integrator.ODEIntegrator.u[1]) - Lx * Lx_terminate_limit >= 0 #|| maximum(exp.(integrator.u[3:N_state:end]) / e_min_log ) >= 5
                terminate!(integrator.ODEIntegrator)
                @show "terminate"
        end
end



end
