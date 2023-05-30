module movie


using GLMakie
using ParticleMesh: TwoDGrid, TwoDGridNotes, TwoDGridMesh
using PiCLES.Operators.core_2D: GetGroupVelocity

import Oceananigans.Utils: prettytime

function init_movie_2D_box_plot(wave_simulation)

    n = Observable(1) # for visualization
    # Ocean vorticity
    grid = wave_simulation.model.grid
    mesh = TwoDGridMesh(grid, skip=1)
    gn = TwoDGridNotes(grid)

    model_time = @lift ($n; wave_simulation.model.clock.time)
    uo, vo = wave_simulation.model.winds
    ocean_wind_u = @lift(uo.(mesh.x, mesh.y, $model_time))
    ocean_wind_v = @lift(vo.(mesh.x, mesh.y, $model_time))


    #ocean_wind = @lift ($n; sqrt.(ocean_wind_u.^2 + ocean_wind_u.^2))
    #ocean_wind = @lift($n; sqrt(vo.(mesh.x, mesh.y, $model_time)^2 + uo.(mesh.x, mesh.y, $model_time)^2))
    ocean_wind = @lift(sqrt.(vo.(mesh.x, mesh.y, $model_time).^2 + uo.(mesh.x, mesh.y, $model_time).^2))
    strength = @lift( vec(sqrt.(vo.(mesh.x, mesh.y, $model_time).^2 + uo.(mesh.x, mesh.y, $model_time).^2)))
    #ocean_wind = @lift ($n; model_time * 2)

    #@info ocean_wind

    wave_energy = @lift ($n; 4 * sqrt.(wave_simulation.model.MovieState[:,:, 1]))
    wave_momentum_x = @lift ($n; wave_simulation.model.MovieState[:, :, 2])
    wave_momentum_y = @lift ($n; wave_simulation.model.MovieState[:, :, 3])
    cx = @lift ($n; GetGroupVelocity(wave_simulation.model.MovieState).c_x)
    cy = @lift ($n; GetGroupVelocity(wave_simulation.model.MovieState).c_y)

    #c_max = @lift ($n; sqrt(GetGroupVelocity(wave_simulation.model.MovieState).c_y^2 + GetGroupVelocity(wave_simulation.model.MovieState).c_x)^2 )
    #c_max = round(maximum(sqrt(GetGroupVelocity(wave_simulation.model.MovieState).c_y^2 + GetGroupVelocity(wave_simulation.model.MovieState).c_x)^2))

    # for testing
    #we = wave_simulation.model.State[:, :, 1]

    # Make figure
    #x, y, z = nodes((Center, Center, Center), grid)
    fig = Figure(resolution=(900, 1200))

    ax_wind = Axis(fig[1, 1], aspect=1, xlabel="x (km)", ylabel="y (km)", title="Winds")
    ax_o  = Axis(fig[1, 2], aspect=1, xlabel="x (km)",  title="Hs")
    ax_mx = Axis(fig[2, 1], aspect=1, xlabel="x (km)", ylabel="y (km)", title="x momentum")
    ax_my = Axis(fig[2, 2], aspect=1, xlabel="x (km)",  title="y momentum")

    ax_cx = Axis(fig[3, 1], aspect=1, xlabel="x (km)", ylabel="y (km)", title="c_x")
    ax_cy = Axis(fig[3, 2], aspect=1, xlabel="x (km)",  title="c_y")

    #hm_i = heatmap!(ax_wind, 1e-3 * x, 1e-3 * y, ice_speed_n)
    hm_wind = heatmap!(ax_wind, 1e-3 * gn.x[1:3:end], 1e-3 * gn.y[1:3:end], ocean_wind, colormap=:jet, colorbar=true, colorrange=(-11, 11))
    # colormap options for heatmap 

    #quiver!(ax_wind, 1e-3 * gn.x, 1e-3 * gn.y, quiver=(ocean_wind_u, ocean_wind_v))#, color=:red, scale_unit=:data, label="wind")
    #strength = ocean_wind_u.^2 .+ ocean_wind_v.^2
    arrows!(ax_wind, 1e-3 * gn.x, 1e-3 * gn.y, ocean_wind_u, ocean_wind_v, arrowsize=10, lengthscale = 0.08, length=strength)#, arrowcolor=strength, linecolor=strength)
    #scatter!(ax_wind, vec(gridmesh.* 1e-3), rotations=0, markersize=20, marker='↑')


    hm_o = heatmap!(ax_o, 1e-3 * gn.x, 1e-3 * gn.y, wave_energy, colormap=:dense, colorrange=(0, 1.5))
    hm_x = heatmap!(ax_mx, 1e-3 * gn.x, 1e-3 * gn.y, wave_momentum_x, colormap=:balance, colorrange=(-0.01, 0.01))
    hm_y = heatmap!(ax_my, 1e-3 * gn.x, 1e-3 * gn.y, wave_momentum_y, colormap=:balance, colorrange=(-0.01, 0.01))

    hm_cx = heatmap!(ax_cx, 1e-3 * gn.x, 1e-3 * gn.y, cx, colormap=:balance, colorrange=(-4, 4))
    hm_cy = heatmap!(ax_cy, 1e-3 * gn.x, 1e-3 * gn.y, cy, colormap=:balance, colorrange=(-4, 4))
    #colormaps

    #colorbar(ax_my)
    cb_wind = Colorbar(fig[1, 0], hm_wind, label="winds [m/s]", tickalignmode=:right)
    cb_wind.alignmode = Mixed(right=1)
    #Colorbar(fig[1, 3], hm_o, label="Wave energy [m^2]")

    Colorbar(fig[1, 4], hm_o, label = "Wave energy [m^2]")
    Colorbar(fig[2, 4], hm_x, label = "Wave momentum x []")
    Colorbar(fig[3, 4], hm_cx, label="Group Velocity [m/s]")

    limits!(ax_wind, (1e-3 * gn.x[1], 1e-3 * gn.x[end]), (1e-3 * gn.y[1], 1e-3 * gn.y[end]))
    # sc
    #hm_o = heatmap(ax_o, 1e-3 * mesh.x, 1e-3 * mesh.y, ocean_wind_v, colormap=:redblue)
    DT = wave_simulation.model.ODEsettings.timestep
    #c_yi = GetGroupVelocity(wave_simulation.model.MovieState).c_y
    #c_xi = GetGroupVelocity(wave_simulation.model.MovieState).c_x
    #CFL = @lift ($n; round(maximum(sqrt(cx[]^2 + cx[]^2))) * DT / gn.dx)

    title = @lift ($n; "Static winds DT=$DT , dx=$(round(gn.dx)), CFL= $(round( maximum(sqrt.(cx[].^2+cy[].^2)) * DT /gn.dx; digits=3 )), \ntime=" * prettytime(wave_simulation.model.clock.time))
    Label(fig[0, :], title)
    display(fig)

    return fig, n
end


end # end of module