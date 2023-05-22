using Plots

"""
plot_results(wave_simulation::Simulation; wind_grid=nothing)
    plot the results of a simulation in a simple way
"""
function plot_results(wave_simulation; wind_grid=nothing, title="")
    output = convert_store_to_tuple(wave_simulation.store, wave_simulation)

    store_waves_energy = output.data[:, :, 1]
    store_waves_mx = output.data[:, :, 2]
    #store_waves_my = output.data[:, :, 3];
    cg = store_waves_energy ./ store_waves_mx ./ 2

    gr(display_type=:inline)

    dx = wave_simulation.model.grid.dx
    DT = wave_simulation.model.ODEsettings.timestep
    # ## make a three panel plot for energy, mx, and cg 
    p1 = heatmap(output.x / dx, output.time / DT, store_waves_energy, levels=20, colormap=:dense)
    p2 = heatmap(output.x / dx, output.time / DT, store_waves_mx, levels=20, colormap=:dense)
    p3 = heatmap(output.x / dx, output.time / DT, cg, levels=40, colormap=:dense)

    if wind_grid != nothing
        # add wind forcing to plot
        contour!(p1, wind_grid.x / dx, wind_grid.t / DT, transpose(wind_grid.u) / 10, levels=5, colormap=:black)
        contour!(p2, wind_grid.x / dx, wind_grid.t / DT, transpose(wind_grid.u) / 200, levels=5, colormap=:black)
        contour!(p3, wind_grid.x / dx, wind_grid.t / DT, transpose(wind_grid.u) / 3, levels=5, colormap=:black)
    end

    plot(p1, p2, p3, layout=(3, 1), size=(600, 1200), title=[title * "\nEnergy" "mx" "cg"], xlabel="x (dx)", ylabel="time (DT)", left_margin=10 * Plots.mm) |> display

end

