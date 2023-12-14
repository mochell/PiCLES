

#using Plots
import Plots as plt
using Setfield, IfElse

using PiCLES.ParticleSystems: particle_waves_v5 as PW

import PiCLES: FetchRelations, ParticleTools
using PiCLES.Operators.core_2D: ParticleDefaults, InitParticleInstance, GetGroupVelocity
using PiCLES.Operators: TimeSteppers
using PiCLES.Simulations
using PiCLES.Operators.TimeSteppers: time_step!, movie_time_step!

using PiCLES.ParticleMesh: TwoDGrid, TwoDGridNotes, TwoDGridMesh
using PiCLES.Models.WaveGrowthModels2D

using Oceananigans.TimeSteppers: Clock, tick!
import Oceananigans: fields
using Oceananigans.Units
import Oceananigans.Utils: prettytime

using Architectures
using GLMakie

using PiCLES.Operators.core_2D: GetGroupVelocity, speed
using PiCLES.Plotting.movie: init_movie_2D_box_plot

using NCDatasets
using Interpolations
using Dates: Dates as Dates

# %%
using Distributions
save_path = "plots/tests/T04_2D_regtest_netCDF/"
mkpath(save_path)

save_path_data = "data/work/T04_2D_regtest_netCDF/"
mkpath(save_path_data)

##### basic parameters
# timestep
DT = 20minutes
# Characterstic wind velocities and std
U10, V10 = 10.0, 10.0

# Define basic ODE parameters
r_g0            = 0.85
Const_ID        = PW.get_I_D_constant()
@set Const_ID.γ = 0.88
Const_Scg       = PW.get_Scg_constants(C_alpha=-1.41, C_varphi=1.81e-5)


function interpolate_winds(ds)
    # define grid based on
    grid = TwoDGrid(ds["x"][end], Int(ceil(ds.attrib["Nx"] )),
        ds["y"][end], Int(ceil(ds.attrib["Ny"] )))
    #grid = TwoDGrid(ds["x"][end], 31, ds["y"][end], 31)
    grid_mesh = TwoDGridMesh(grid, skip=1)
    gn = TwoDGridNotes(grid)

    # define time
    time_rel = (ds["time"][:] - ds["time"][1]) ./ convert(Dates.Millisecond, Dates.Second(1))
    T = 1day#hours#time_rel[end]

    nodes = (ds["x"][:], ds["y"][:], time_rel)
    u_grid = LinearInterpolation(nodes, permutedims(ds["u10m"][:], [1, 2, 3]), extrapolation_bc=Periodic())
    v_grid = LinearInterpolation(nodes, permutedims(ds["v10m"][:] .+ 0.1, [1, 2, 3]), extrapolation_bc=Periodic())

    return grid, grid_mesh, gn, T, u_grid, v_grid
end


# define ODE system and parameters
#particle_system = PW.particle_equations(u, v, γ=0.88, q=Const_ID.q);
 

default_ODE_parameters = (r_g=r_g0, C_α=Const_Scg.C_alpha,
    C_φ=Const_ID.c_β, C_e=Const_ID.C_e, g=9.81)


Revise.retry()
# Default initial conditions based on timestep and chaeracteristic wind velocity
WindSeamin = FetchRelations.get_minimal_windsea(U10, V10, DT)
default_particle = ParticleDefaults(WindSeamin["lne"], WindSeamin["cg_bar_x"], WindSeamin["cg_bar_y"], 0.0, 0.0)


function make_reg_test(wave_model, save_path; plot_name="dummy", N=36)

    ### build Simulation
    wave_simulation = Simulation(wave_model, Δt=DT, stop_time=wave_model.ODEsettings.total_time)#1hours)
    initialize_simulation!(wave_simulation)

    # run simulation & save simulation
    init_state_store!(wave_simulation, save_path)
    #run!(wave_simulation, cash_store=true, debug=true)
    run!(wave_simulation, store=true, cash_store=false, debug=false)
    close_store!(wave_simulation)
    # # or, alternatively, make movie
    # fig, n = init_movie_2D_box_plot(wave_simulation, name_string="T01")

    # #wave_simulation.stop_time += 1hour
    # #N = 36
    # #plot_name = "dummy"
    # record(fig, save_path * plot_name * ".gif", 1:N, framerate=10) do i
    #     @info "Plotting frame $i of $N..."
    #     @info wave_simulation.model.clock
    #     movie_time_step!(wave_simulation.model, wave_simulation.Δt)

    #     n[] = 1
    # end

end


# %%
# loop over U10 and V10 range
#case_list = ["Test02_2D", "Test03_2D", "Test04_2D"]
#case_list = ["Test01_2D", "Test02_2D"]#, "Test03_2D", "Test04_2D", "Test05_2D", "Test06_2D", "Test07_2D"]
case_list = ["Test03_2D", "Test04_2D", "Test05_2D", "Test06_2D", "Test07_2D"]
#for I in CartesianIndices(gridmesh)
for case in case_list
    # load netCDF file
    ncfile = "data/work/wind_data/" * case * ".nc"
    ds = Dataset(ncfile, "r")

    grid, grid_mesh, gn, T, u_grid, v_grid = interpolate_winds(ds)

    u(x, y, t) = u_grid(x, y, t)
    v(x, y, t) = v_grid(x, y, t)
    winds = (u=u, v=v)

    #winds, u, v  =convert_wind_field_functions(u_func, v_func, x, y, t)
    particle_system = PW.particle_equations(u, v, γ=0.88, q=Const_ID.q)

    # ... and ODESettings
    ODE_settings =  PW.ODESettings(
        Parameters=default_ODE_parameters,
        # define mininum energy threshold
        log_energy_minimum=WindSeamin["lne"],
        #maximum energy threshold
        log_energy_maximum=log(27),#log(17),  # correcsponds to Hs about 16 m
        saving_step=DT,
        timestep=DT,
        total_time=T,
        adaptive=true,
        dt=1e-3, #60*10, 
        dtmin=1e-4, #60*5, 
        force_dtmin=true,
        callbacks=nothing,
        save_everystep=false)

    ## Define wave model
    wave_model = WaveGrowthModels2D.WaveGrowth2D(; grid=grid,
        winds=winds,
        ODEsys=particle_system,
        ODEsets=ODE_settings,  # ODE_settings
        ODEinit_type="wind_sea",  # default_ODE_parameters
        periodic_boundary=false,
        boundary_type="same",
        minimal_particle=FetchRelations.MinimalParticle(U10, V10, DT), #
        #minimal_state=FetchRelations.MinimalState(2, 2, DT) * 1,
        movie=true)


    NN = Int(floor(wave_model.ODEsettings.total_time / wave_model.ODEsettings.timestep))
    mkpath(save_path_data * case)
    make_reg_test(wave_model, save_path_data * case, plot_name=case, N=NN)
end
