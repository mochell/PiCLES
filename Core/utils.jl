
##### Callbacks ######

"""
wrap_pos!(xx::Float64, L::Float64)
makes periodic boundary conditions. If position xx exceeds 0 or L it is wrapped around
returns new position and if either or not the wrapping appends (bool)
"""
function wrap_pos!(xx::Float64, L::Float64)
    wrap_flag = false
    if xx > L
        xx = xx - L
        wrap_flag = true
    elseif xx < 0
        xx = xx + L
        wrap_flag = true
    end

    return xx, wrap_flag
end

"""
Checks if particle's PI position PI.ODEIntegrator.u[1,2] exceeeds the domain limits Lx, Ly
        Lx, LY are global variables
        If the PI's postions eceeeds the boundary they warpped around and the ODEIntegrator state is modified.
"""
function periodic_BD_single_PI!(PI::AbstractParticleInstance, Lx::Float64)

    #@printf "PI periodic condition called"
    ui = copy(PI.ODEIntegrator.u)
    ui[3], wrap_pos_PI1 = wrap_pos!(ui[3], Lx)
    #ui[2], wrap_pos_PI2 = wrap_pos!(ui[2], Ly)

    if wrap_pos_PI1 #|| wrap_pos_PI1
        #@show wrap_pos_PI1 , wrap_pos_PI2
        #@printf "wrap pos PI"
        #@show PI.ODEIntegrator.u[1] - ui[1]
        #@show PI.ODEIntegrator.u[2] - ui[2]
        set_u!(PI.ODEIntegrator, ui)
        u_modified!(PI.ODEIntegrator, true)
    end
    nothing #return PI
end

# """
# Checks if
# """
# @everywhere function periodic_BD_single!(S::SharedMatrix)
#         S.u[1], wrap_pos_PI = wrap_pos!(S.u[1], Lx)
#         S.u[2], wrap_pos_PI = wrap_pos!(S.u[2], Ly)
#         nothing
# end

"""
return the mean position
"""
function show_pos!(PI)
    @show "show current position" PI.position_ij, PI.position_xy / dx, PI.ODEIntegrator.u[3] / dx, PI.ODEIntegrator.u[2]
end


""" (debubugging function) returns the domains normalized position of the partivle  """
function show_pos!(integrator, G)
    @show (integrator.ODEIntegrator.u[3] .- G.xmin) / grid1d.dx
end

# the following function are another way to have wrapping boundary conditions.
function periodic_condition_x(u, t, integrator)
    u[3] > 0
end
