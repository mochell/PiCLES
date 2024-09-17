
using Pkg

test_flag = true

Pkg.activate("../")

using Base.Threads
@info "Num. of threads", nthreads()

#using Plots
#import Plots as plt

using Setfield
using PiCLES.ParticleSystems: particle_waves_v5 as PW

import PiCLES: FetchRelations, ParticleTools
using PiCLES.Operators.core_2D: ParticleDefaults, InitParticleInstance, GetGroupVelocity
using PiCLES.Simulations
using PiCLES.Operators.TimeSteppers: time_step!, movie_time_step!

using PiCLES.ParticleMesh: TwoDGrid, TwoDGridNotes, TwoDGridMesh
using PiCLES.Models.WaveGrowthModels2D

using Oceananigans.TimeSteppers: Clock, tick!
using Oceananigans.Units
import Oceananigans.Utils: prettytime

using PiCLES.Architectures

using Interpolations
using Dates: Dates as Dates

using DifferentialEquations

using Revise
using Printf

# %%
if ~test_flag
    if length(ARGS) < 2
        println("Please provide a case [1] and a multiplyer [2] as arguments.")
        exit(1)
    end

else
    ARGS = ["32", "0"]
end

# wind data case that is used
#case = 32
case = parse(Int, ARGS[1])

# multiplyer for grid size (2^multiplyer)
#multiplyer = 0
#multiplyer = parse(Float64, ARGS[2])  # assuming multiplyer is a floating point number
multiplyer = parse(Int, ARGS[2])  # assuming multiplyer is an integer

#load_path = base_path * "data/work/wind_data_scaling/"
#load_path = "/glade/work/mhell/2022_particle_waves/wind_data_scaling/"

# save_path = "plots/tests/T04_2D_regtest_netCDF/"
# mkpath(save_path)

# save_path_data = "data/work/T04_2D_regtest_netCDF/"
# mkpath(save_path_data)

##### basic parameters
# %%
# timestep
DT = 10minutes
#DT = 20minutes # only use for strong scaling test comparison 

# Characterstic wind velocities and std
U10, V10 = 10.0, 10.0

u_func(x::Number, y::Number, t::Number) = 10.0 #* sin(t / (6 * 60 * 60 * 2π)) * sin(x / 50e3) * sin(y / 50e3)
v_func(x::Number, y::Number, t::Number) = 10.0 #* cos(t / (6 * 60 * 60 * 2π)) * sin(x / 50e3) * sin(y / 50e3)

# Define basic ODE parameters
r_g0 = 0.85
Const_ID = PW.get_I_D_constant()
@set Const_ID.γ = 0.88
Const_Scg = PW.get_Scg_constants(C_alpha=-1.41, C_varphi=1.81e-5)

# function interpolate_winds(ds, multiplyer=0 )

#     Nxx = Int(ceil(ds.attrib["Nx"]  * sqrt(2)^multiplyer ))
#     Nyy = Int(ceil(ds.attrib["Ny"]  * sqrt(2)^multiplyer ))

#     # define grid based on
#     grid = TwoDGrid(ds["x"][end], Nxx, ds["y"][end], Nyy)
#     grid_mesh = TwoDGridMesh(grid, skip=1)
#     gn = TwoDGridNotes(grid)

#     # define time
#     time_rel = (ds["time"][:] - ds["time"][1]) ./ convert(Dates.Millisecond, Dates.Second(1))
#     T = 1day#hours#time_rel[end]

#     nodes = (ds["x"][:], ds["y"][:], time_rel)
#     u_grid = LinearInterpolation(nodes, permutedims(ds["u10m"], [1, 2, 3]), extrapolation_bc=Flat());
#     v_grid = LinearInterpolation(nodes, permutedims(ds["v10m"], [1, 2, 3]), extrapolation_bc=Flat());

#     return grid, grid_mesh, gn, T, u_grid, v_grid
# end


# define ODE system and parameters
# particle_system = PW.particle_equations(u, v, γ=0.88, q=Const_ID.q);

default_ODE_parameters = (r_g=r_g0, C_α=Const_Scg.C_alpha,
    C_φ=Const_ID.c_β, C_e=Const_ID.C_e, g=9.81)

Revise.retry()
# Default initial conditions based on timestep and chaeracteristic wind velocity
WindSeamin = FetchRelations.MinimalWindsea(U10, V10, DT)
default_particle = ParticleDefaults(WindSeamin["lne"], WindSeamin["cg_bar_x"], WindSeamin["cg_bar_y"], 0.0, 0.0)

# load netCDF file
# ncfile = load_path * case * ".nc"


# println("Input file ", ncfile, "\n")
# println("Multiplyer ", multiplyer, "\n")
# println("Num. of threads ", nthreads(), "\n")

# ds = Dataset(ncfile, "r");

# grid, grid_mesh, gn, T, u_grid, v_grid = interpolate_winds(ds, multiplyer);

dx = 20e3 # basic grid size 20 km 

grid = TwoDGrid(dx * (case - 1), case, dx * (case - 1), case)
grid_mesh = TwoDGridMesh(grid, skip=1);
gn = TwoDGridNotes(grid);

T = 1day#hours#time_rel[end]

u(x, y, t) = u_func(x, y, t)
v(x, y, t) = v_func(x, y, t)
winds = (u=u, v=v)

#winds, u, v  =convert_wind_field_functions(u_func, v_func, x, y, t)
particle_system = PW.particle_equations(u, v, γ=0.88, q=Const_ID.q)

# ... and ODESettings
ODE_settings = PW.ODESettings(
    Parameters=default_ODE_parameters,
    # define mininum energy threshold
    log_energy_minimum=WindSeamin["lne"],
    #maximum energy threshold
    solver=DP5(),
    log_energy_maximum=log(27),#log(17),  # correcsponds to Hs about 16 m
    saving_step=6days,
    timestep=DT,
    total_time=T,
    adaptive=true,
    dt=10, #60*10, 
    dtmin=1, #60*5, 
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


# for saving data
# when saving data
# save_path_select = save_path_data
# mkpath(save_path_select * case)

### build Simulation
wave_simulation = Simulation(wave_model, Δt=DT, stop_time=20minutes)
initialize_simulation!(wave_simulation)

# run simulation
run!(wave_simulation, cash_store=false, debug=false);
reset_simulation!(wave_simulation);

# run actual test
wave_simulation.stop_time = 6hour
@time run!(wave_simulation, cash_store=false, debug=false);

println("problem size: N=", Int(grid.Nx * grid.Ny), "\n")
println("----------------------------------------------\n")


# save simulation
# statestore version
# init_state_store!(wave_simulation, save_path)
# run!(wave_simulation, store=true, cash_store=false, debug=false)
# close_store!(wave_simulation)

# chash store version
#run!(wave_simulation, cash_store=true, debug=true)

