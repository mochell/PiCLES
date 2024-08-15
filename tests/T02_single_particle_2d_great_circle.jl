using Pkg
Pkg.activate("PiCLES/")

using DifferentialEquations
using Plots
using Setfield
using IfElse

using PiCLES.ParticleSystems: particle_waves_v5 as PW
using PiCLES.Utils: Init_Standard

import PiCLES: FetchRelations, ParticleTools
using PiCLES.Operators.core_2D: InitParticleInstance
using Oceananigans.Units

using PiCLES.Grids.TripolarGridMOM6: TripolarGridMOM6

using PiCLES.Grids: SphericalPropagationCorrection

plot_path_base = "plots/tests/T04_2D_single_particle/"
mkpath(plot_path_base)

# %%
load_path = "PiCLES/src/Grids/files/";
gridd = TripolarGridMOM6.MOM6GridMesh(load_path * "ocean_hgrid_221123.nc", 6; mask_radius=5);

using Makie
using PiCLES.Plotting: PlotState_DoubleGlobe, PlotState_SingleGlobe, PlotState_DoubleGlobeSeam, OrthographicTwoMaps, OrthographicTwoMapsSeam

data = SphericalPropagationCorrection.(gridd.data.y, -10.0)

fig = Figure(size=(900, 500), fontsize=22)
OrthographicTwoMapsSeam(fig, gridd.data.x[1:end, 1:end-10], gridd.data.y[1:end, 1:end-10], data[1:end, 1:end-10])
# heatmap(gridd.data.angle_dx)
fig

# %%
deg_x = -89.9:0.1:89.9
Gg = SphericalPropagationCorrection.(deg_x)
G_test = [Gg_i(10) for Gg_i in Gg]
Plots.plot(deg_x, G_test)



# define Particle Position
# %%
ij = (180, 86)
ij_mesh, gridstats = gridd.data[ij[1], ij[2]], gridd.stats
#ij_mesh.y = 0.0
x_lon, y_lat = ij_mesh.x, ij_mesh.y  
SGC = SphericalPropagationCorrection(ij_mesh, gridstats)


Revise.retry()
u(x::Number, y::Number, t::Number) = 0.0#(10.0 * cos(t / (3 * 60 * 60 * 2π)) + 0.1) + x * 0 + y * 0
v(x::Number, y::Number, t::Number) = 0.0#-(10.0 * sin(t / (3 * 60 * 60 * 2π)) + 0.1) + x * 0 + y * 0

DT = 30minutes

using PiCLES.Operators.core_2D: ParticleDefaults

ParticleState, default_ODE_parameters, WindSeamin, Const_ID = Init_Standard(15.0, 00.0, 4hour)
#FetchRelations.get_initial_windsea(15.0, 0.0, 4hour, particle_state=true)
ParticleState = ParticleDefaults(FetchRelations.get_initial_windsea(1.0, 10.0, 4hour, particle_state=true))

#particle_system = PW.particle_equations(u, v, γ=Const_ID.γ, q=Const_ID.q)

particle_system = PW.particle_equations(u, v, γ=Const_ID.γ, q=Const_ID.q,
    propagation=true,
    input=false,
    dissipation=false,
    peak_shift=false,
    direction=false,
);


used_ODE_params = (default_ODE_parameters..., x=x_lon, y=y_lat, PC=SphericalPropagationCorrection(ij_mesh, gridstats))

# test execution
ODE_settings = PW.ODESettings(
    Parameters=used_ODE_params,
    # define mininum energy threshold
    log_energy_minimum=log(WindSeamin["E"]),
    #maximum energy threshold
    log_energy_maximum=log(17),  # correcsponds to Hs about 16 m
    saving_step=10minutes,
    timestep=DT,
    total_time=T = 30days,
    save_everystep=false,
    maxiters=1e4,
    adaptive=true,
    dt=10,#60*10, 
    dtmin=1,#60*5, 
    force_dtmin=true,)

PI = InitParticleInstance(particle_system, ParticleState, ODE_settings, ij, (ij_mesh.x, ij_mesh.y), false, true)


R = 6.3710E+6
PI.position_xy
M_geocoords = @SArray [
    1/(R*cos(y_lat * pi / 180)) 0;
    0 1/R
]
PI.ODEIntegrator.u[4:5]
# propagation is in meters zonal and meridinal, there is no projection


function time_step_local!(PI, DT)
    "take 1 step over DT"

    M = PI.Parameters.M
    x_lon = PI.position_xy[1] 
    y_lat = PI.position_xy[2]

    used_ODE_params = (default_ODE_parameters..., x=x_lon, y=y_lat, PC=SphericalPropagationCorrection(ij_mesh, gridstats))

    # test execution
    ODE_settings = PW.ODESettings(
        Parameters=used_ODE_params,
        # define mininum energy threshold
        log_energy_minimum=log(WindSeamin["E"]),
        #maximum energy threshold
        log_energy_maximum=log(17),  # correcsponds to Hs about 16 m
        saving_step=10minutes,
        timestep=DT,
        total_time=T = 30days,
        save_everystep=false,
        maxiters=1e4,
        adaptive=true,
        dt=10,#60*10, 
        dtmin=1,#60*5, 
        force_dtmin=true,)

    PI = InitParticleInstance(particle_system, ParticleState, ODE_settings, ij, (x_lon, y_lat), false, true)

    #@info "proposed dt", get_proposed_dt(PI.ODEIntegrator) / 60
    step!(PI.ODEIntegrator, DT, true)

    # #@info "u:", PI.ODEIntegrator.u
    # #clock_time += DT
    # last_t = PI.ODEIntegrator.t

    # ## define here the particle state at time of resetting
    # #ui = [log(exp(PI.ODEIntegrator.u[1]) * 0.5), PI.ODEIntegrator.u[2] / 2, PI.ODEIntegrator.u[3] / 2, 0.0, 0.0]
    # #ui = [lne_local, cg_u_local, cg_v_local, 0.0, 0.0]
    # ui = PI.ODEIntegrator.u

    # #ui = [PI.ODEIntegrator.u[1], PI.ODEIntegrator.u[2], PI.ODEIntegrator.u[3], 0.0, 0.0]
    # WindSeamin = FetchRelations.get_initial_windsea(u(0.0, 0.0, last_t), v(0.0, 0.0, last_t), DT / 2)
    # ui = [log(WindSeamin["E"]), WindSeamin["cg_bar_x"], WindSeamin["cg_bar_y"], 0.0, 0.0]

    # set_u_and_t!(PI.ODEIntegrator, ui, last_t)
    # # #set_u!(PI.ODEIntegrator, ui)
    # #reinit!(PI.ODEIntegrator, ui, erase_sol=false, reset_dt=true, reinit_cache=true)
    # #reinit!(PI3.ODEIntegrator, ui, erase_sol=false, reset_dt=true, reinit_cache=true)

    # # #set_t!(PI.ODEIntegrator, last_t )
    # u_modified!(PI.ODEIntegrator, true)

    # # add_saveat!(PI.ODEIntegrator, PI.ODEIntegrator.t)
    # # savevalues!(PI.ODEIntegrator)

    return PI
end


for i in Base.Iterators.take(PI.ODEIntegrator, 600)
    PI = time_step_local!(PI, DT)
    #step!(PI.ODEIntegrator, DT, true)
end

#2.6 * 20 * DT /1e3

PID = ParticleTools.ParticleToDataframe(PI)
PI.ODEIntegrator.u

gr(display_type=:inline)
# plit each row in PID and a figure
tsub = range(start=1, stop=length(PID[:, 1]), step=10)
p3 = Plots.plot( x_lon  .+  PID[tsub, 5] / (110e3 * cos(y_lat)), y_lat .+  PID[tsub, 6] / 110e3, marker=3, title="position on Sphere ", ylabel="lat", xlabel= "lon", label="v4") #|> display

Plots.plot!(p3,  [x_lon]  , [y_lat] , marker=10, title="lon ", ylabel="lat", label="origin") #|> display

p3
# subtitle = "u$(U10)_v$(V10)_reset_to_windsea_dt$(DT)"
# subtitle = "u10_v10_reset_to_windsea_dt$(DT)"
# savefig(joinpath(plot_path_base, subtitle*"_continous_foreward.png"))



# %%

