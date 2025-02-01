using Pkg
Pkg.activate("PiCLES/")

using CairoMakie, GeoMakie

fig = Figure()
ga = GeoAxis(
    fig[1, 1]; # any cell of the figure's layout
    dest="+proj=wintri", # the CRS in which you want to plot
)
lines!(ga, GeoMakie.coastlines()) # plot coastlines from Natural Earth as a reference
# You can plot your data the same way you would in Makie
scatter!(ga, -120:15:120, -60:7.5:60; color=-60:7.5:60, strokecolor=(:black, 0.2))
fig



# %%

lons = -180:180
lats = -90:90
field = [exp(cosd(l)) + 3(y / 90) for l in lons, y in lats]


function simple_geomap(lons, lats, field)
    fig = Figure()
    ax = GeoAxis(fig[1, 1])
    sp = surface!(ax, lons, lats, field; shading=NoShading)
    # lines!(ax, GeoMakie.coastlines())

    land = GeoMakie.land()

    poly!(ax, land, color=:black)

    fig
end

simple_geomap(lons, lats, field)

# %%

using GeoMakie, CairoMakie

xs = LinRange(311.9, 391.1, 30)
ys = LinRange(-23.6, 24.8, 20)

us = @. 1 * (2 * cos(2 * deg2rad(xs) + 3 * deg2rad(ys' + 30))^2)
vs = @. 2 * cos(6 * deg2rad(xs)) .+ ys' * 0 # that last part is just to establish the shape


pole_longitude = 177.5
pole_latitude = 37.5
arrow_crs = "+proj=ob_tran +o_proj=latlon +o_lon_p=0 +o_lat_p=$(pole_latitude) +lon_0=$(180+pole_longitude) +to_meter=$(deg2rad(1) * 6378137.0)"

f, a, p = arrows(
    xs, ys, us, vs;
    arrowsize=4,
    source=arrow_crs,
    axis=(; type=GeoAxis, dest="+proj=ortho +lon_0=-10 +lat_0=45")
)

# %%
ep = surface!(a,
    -180 .. 180, -90 .. 90,
    zeros(axes(rotr90(GeoMakie.earth())));
    shading=NoShading, color=rotr90(GeoMakie.earth())
)

translate!(ep, 0, 0, -1)
f