module WindEmulator

export WindEmulationFunction, IdealizedWindGrid, wind_interpolator
export slopped_blob, contant_winds

using Interpolations

function IdealizedWindGrid(u_func::Function, dims::NamedTuple, steps::NamedTuple)

    xi = range(0, dims.Lx, step=steps.dx)
    ti = range(0, dims.T, step=steps.dt)

    u_func_gridded = [ u_func(xii, tii) for xii in xi, tii in ti]
    return (u= u_func_gridded,  x=xi, t=ti)
end


function wind_interpolator(wind_grid)
    if haskey(wind_grid, :y)
            nodes = (wind_grid.x, wind_grid.y, wind_grid.t)
    else
            nodes = (wind_grid.x, wind_grid.t)
    end

    if haskey(wind_grid, :u)
            u_grid = linear_interpolation( nodes , wind_grid.u ,extrapolation_bc=Periodic())
            #u_grid = CubicSplineInterpolation( (xi, ti), u_func_gridded; bc=Line(OnGrid()), extrapolation_bc=Flat())
            #u_grid = interpolate((xi, ti), u_func_gridded, Gridded(Linear()))
    else
            u_grid = nothing
    end

    if haskey(wind_grid, :v)
            v_grid = linear_interpolation( nodes , wind_grid.v ,extrapolation_bc=Periodic())
    else
            v_grid = nothing
    end


    tcol= (u=u_grid, v= v_grid)
    # Create a new NamedTuple without the entries that are nothing
    return (; [p for p in pairs(tcol) if ~isnothing(p[2])]...)
end




abstract type WindEmulationFunction2 <: Function end

function slopped_blob(x, t, U10, V, T , x_scale, t_scale; x0=300e3)# <: WindEmulationFunction
    return 0.5+ U10 * ( exp(- ( ( x - (x0 + t * V)  )./x_scale).^2) .*  exp(- ( ( t- (T/2)  )./  t_scale ).^2) )#.+ 3 * sin.(x *π/Ly/0.5 .+ y *π/Ly/0.5)
end


function contant_winds(x, t, U10)# <: WindEmulationFunction
    return x.*0+U10 + t *0 
end


#u_func(x, t) = 3 * exp(- ( ( x-25e3 + 20e3* t/T )./10e3).^2) #.+ 3 #.+ 3 * sin.(x *π/Ly/0.5 .+ y *π/Ly/0.5)
#u_func(x, t) = 3 * exp(- ( ( x-25e3  )./10e3).^2) .+ t *0 #.+ 3 #.+ 3 * sin.(x *π/Ly/0.5

#u_func(x, y, t) = y *0 .- x .*2/50e3 .+ 3 #.+ 3 * sin.(x *π/Ly/0.5 .+ y *π/Ly/0.5)
#u_func(x, y, t) = y *0 .- x .*0 .+ 3 #.+ 3 * sin.(x *π/Ly/0.5 .+ y *π/Ly/0.5)
#u_func(x, t) = x *0 .+ IfElse.ifelse.( x .< Lx*2.0/3.0, 3, 0.1)

# u_func(x, y, t) = y *0 .+ IfElse.ifelse.(    x .< Lx/2.0,
#                                 2 .+ 1, *sin.(x *π/Lx/2),
#                                 0 .* x + 0.1
#                                 )



end