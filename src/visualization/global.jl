
using CairoMakie, GeoMakie

function OrthographicTwoMaps(fig, lons, lats, field; title="", kargs=nothing)

    ax_left =  GeoAxis(fig[:, 1]; dest="+proj=ortho +lon_0=-20 +lat_0=40", xticks=-180:30:180, yticks=-90:30:90)
    ax_right = GeoAxis(fig[:, 2]; dest="+proj=ortho +lon_0=200 +lat_0=-40", xticks=-180:30:180, yticks=-90:30:90)

    if kargs == nothing
        kargs = (shading=NoShading, colormap=:dense)#, colorrange=range)
    end
    #sp = GeoMakie.surface!(ax, lons, lats, wave_model.ParticleCollection.boundary; shading=NoShading)
    sp1 = GeoMakie.surface!(ax_left, lons, lats, field; kargs...)
    sp2 = GeoMakie.surface!(ax_right, lons, lats, field; kargs...)

    for ax in [ax_left, ax_right]
        lines!(ax, GeoMakie.coastlines(), color=:black, linewidth=0.5)
    end

    cb = Colorbar(fig[:, 3], sp1, label=title, height=Relative(0.7))

    # Add title
    Label(fig[0, :], title)##, size=24, tellwidth=false)
    return sp1, sp2, cb
end

function OrthographicTwoMapsSeam(fig, lons, lats, field; title="", kargs=nothing)

    ax_left  = GeoAxis(fig[:, 1]; dest="+proj=ortho +lon_0=75 +lat_0=90",  xticks=-180:30:180, yticks=-90:30:60)
    ax_right = GeoAxis(fig[:, 2]; dest="+proj=ortho +lon_0=75 +lat_0=0" ,  xticks=-180:30:180, yticks=-90:30:60)

    if kargs == nothing
        kargs = (shading=NoShading, colormap=:dense)#, colorrange=range)
    end
    #sp = GeoMakie.surface!(ax, lons, lats, wave_model.ParticleCollection.boundary; shading=NoShading)
    sp1 = GeoMakie.surface!(ax_left, lons, lats, field; kargs...)
    sp2 = GeoMakie.surface!(ax_right, lons, lats, field; kargs...)

    for ax in [ax_left, ax_right]
        lines!(ax, GeoMakie.coastlines(), color=:black, linewidth=0.5)
    end

    try
        cb = Colorbar(fig[:, 3], sp1, label=title, height=Relative(0.7))
    catch
        cb = 0
    end

    # Add title
    Label(fig[0, :], title)##, size=24, tellwidth=false)
    return sp1, sp2
end


function sym_range(z)
    max_abs_value = maximum(abs, z)
    (-max_abs_value, max_abs_value)
end

function zero_range(z)
    max_abs_value = maximum(abs, z[.!isnan.(z)])
    (0, max_abs_value)
end


function PlotState_DoubleGlobe(model)

    lons, lats = model.grid.data.x, model.grid.data.y

    fig = Figure(size=(900, 1200), fontsize=22)

    #Hs = 4 * sqrt.(wave_model.State[:, :, 1]);
    Energy = model.State[:, :, 1]

    kargs = (shading=NoShading, colormap=:dense, colorrange=zero_range(Energy))
    ax_energy = OrthographicTwoMaps(fig[1, :], lons, lats, Energy; title="Energy", kargs=kargs)

    mom1 = model.State[:, :, 2]
    mom_range = sym_range(model.State[:, :, 2:3])
    kargs = (shading=NoShading, colormap=:balance, colorrange=mom_range)
    ax_mom1 = OrthographicTwoMaps(fig[2, :], lons, lats, mom1; title="x - Momentum", kargs=kargs)

    mom2 = model.State[:, :, 3]
    kargs = (shading=NoShading, colormap=:balance, colorrange=mom_range)
    ax_mom3 = OrthographicTwoMaps(fig[3, :], lons, lats, mom2; title="y - Momentum", kargs=kargs)

    return fig

end


function PlotState_DoubleGlobeSeam(model; scaled=true)

    lons, lats = model.grid.data.x, model.grid.data.y

    fig = Figure(size=(900, 1200), fontsize=22)

    #Hs = 4 * sqrt.(wave_model.State[:, :, 1]);
    Energy = model.State[:, :, 1]
    Energy[model.grid.data.mask.==0] .= NaN
    Energy[model.grid.data.mask.==2] .= NaN
    if scaled
        kargs = (shading=NoShading, colormap=:dense, colorrange=zero_range(Energy))
    else
        kargs = (shading=NoShading, colormap=:balance)#, colorrange=range)
    end
    ax_energy = OrthographicTwoMapsSeam(fig[1, :], lons, lats, Energy; title="Energy", kargs=kargs)


    mom1 = model.State[:, :, 2]
    mom1[model.grid.data.mask.==0] .= NaN
    mom1[model.grid.data.mask.==2] .= NaN

    mom_range = sym_range(model.State[:, :, 2:3])

    if scaled
        kargs = (shading=NoShading, colormap=:balance, colorrange=mom_range)
    else
        kargs = (shading=NoShading, colormap=:balance)#, colorrange=range)
    end

    ax_mom1 = OrthographicTwoMapsSeam(fig[2, :], lons, lats, mom1; title="zonal - Momentum", kargs=kargs)

    mom2 = model.State[:, :, 3]
    mom2[model.grid.data.mask.==0] .= NaN
    mom2[model.grid.data.mask.==2] .= NaN

    if scaled
        kargs = (shading=NoShading, colormap=:balance, colorrange=mom_range)
    else
        kargs = (shading=NoShading, colormap=:balance)#, colorrange=range)
    end

    ax_mom3 = OrthographicTwoMapsSeam(fig[3, :], lons, lats, mom2; title="meriod. - Momentum", kargs=kargs)


    return fig

end


# --------------------

function StandardSingleMap(fig, lons, lats, field; title="", kargs=nothing)

    ax_left = GeoAxis(fig;)
    # ax_right = GeoAxis(fig[:, 2];)

    if kargs == nothing
        kargs = (shading=NoShading, colormap=:dense)#, colorrange=range)
    end
    #sp = GeoMakie.surface!(ax, lons, lats, wave_model.ParticleCollection.boundary; shading=NoShading)
    sp1 = GeoMakie.surface!(ax_left, lons, lats, field; kargs...)
    # sp2 = GeoMakie.surface!(ax_right, lons, lats, field; kargs...)

    for ax in [ax_left]
        lines!(ax, GeoMakie.coastlines(), color=:black, linewidth=0.5)
    end

    cb = Colorbar(fig[3 ,:], sp1, label=title, vertical=false, flipaxis = false, width=Relative(2.5) )# orientation=:horizontal)

    # Add title
    Label(fig[0, :], title)##, size=24, tellwidth=false)
    return sp1, cb
end

# make a meshgrid
function make_global_meshgrid(; dx=0.5, dy=0.5)

    x = -180:dx:180
    y = -90:dy:90

    XX = transpose(reshape(repeat(x, inner=length(y)), length(y), length(x)))
    YY = transpose(reshape(repeat(y, outer=length(x)), length(y), length(x)))
    return XX, YY
end


function PlotState_SingleGlobe(model; winds=true)

    lons, lats = model.grid.data.x, model.grid.data.y
    fig = Figure(size=(1000, 800), fontsize=20)

    if winds
        XX, YY = make_global_meshgrid()

        u_speed = sqrt.(model.winds.u.(XX, YY, model.clock.time) .^ 2 + model.winds.v.(XX, YY, model.clock.time) .^ 2)
        kargs = (shading=NoShading, colormap=:dense, colorrange=zero_range(u_speed))
        ax_winds = StandardSingleMap(fig[1, 1:2], XX, YY, u_speed; title="Wind Speed", kargs=kargs)
    end

    #Hs = 4 * sqrt.(wave_model.State[:, :, 1]);
    Energy = model.State[:, :, 1]
    kargs = (shading=NoShading, colormap=:dense, colorrange=zero_range(Energy))
    ax_energy = StandardSingleMap(fig[1, 3:4], lons, lats, Energy; title="Energy", kargs=kargs)

    mom1 = model.State[:, :, 2]
    mom_range = sym_range(model.State[:, :, 2:3])
    kargs = (shading=NoShading, colormap=:balance, colorrange=mom_range)
    ax_mom1 = StandardSingleMap(fig[2, 1:2], lons, lats, mom1; title="x - Momentum", kargs=kargs)

    mom2 = model.State[:, :, 3]
    kargs = (shading=NoShading, colormap=:balance, colorrange=mom_range)
    ax_mom3 = StandardSingleMap(fig[2, 3:4], lons, lats, mom2; title="y - Momentum", kargs=kargs)

    return fig

end

