ENV["JULIA_INCREMENTAL_COMPILE"]=true
#using ModelingToolkit, DifferentialEquations
using ModelingToolkit#: @register_symbolic
#using Plots
import Plots as plt
using Setfield

push!(LOAD_PATH, joinpath(pwd(), "code/"))
push!(LOAD_PATH, joinpath(pwd(), "code/Core"))
push!(LOAD_PATH, joinpath(pwd(), "code/Simulations"))

import particle_waves_v4#: particle_equations, ODESettings
PW4 = particle_waves_v4

import FetchRelations, ParticleTools
using core_2D: ParticleDefaults, InitParticleInstance, InitParticleState
using TimeSteppers

using ParticleMesh: TwoDGrid, TwoDGridNotes
using WaveGrowthModels2D

using Simulations

using Oceananigans.TimeSteppers: Clock, tick!
import Oceananigans: fields
using Oceananigans.Units
import Oceananigans.Utils: prettytime

using Architectures
using GLMakie

# %%
save_path = "data/2D_test/"
mkpath(save_path)

# % Parameters
U10,V10           = 5.0, 4.0
dt_ODE_save       = 60 * 30 # 3 min
DT                = Float64(60 * 60) * 12 # seconds
T                 = 6 * 24 * 60 * 60 # seconds

# version 3
r_g0              = 0.85

# function to define constants 
Const_ID = PW4.get_I_D_constant()
@set Const_ID.γ = 0.88
Const_Scg = PW4.get_Scg_constants(C_alpha=- 1.41, C_varphi=1.81e-5)
Const_ID
@register_symbolic u(x, y, t)
@register_symbolic v(x, y, t)
# u(x, y, t) = 0.01 + U10 * sin(t / (6 * 60 * 60 * 2π)) * sin(x / 50e3) * sin(y / 50e3)
# v(x, y, t) = 0.01 - V10 * cos(t / (6*60*60 * 2π) ) * sin(x / 50e3) * sin(y / 50e3)
u_std= 30e3
v_std= 30e3
u(x, y, t) = 0.01 + U10  * exp( - (x - 30e3)^2/ u_std^2) * exp( - (y - 30e3)^2/ v_std^2) #* sin(t / (2 * 60 * 60 * 2π))
v(x, y, t) = 0.01 + V10  * exp( - (x - 30e3)^2/ u_std^2) * exp( - (y - 30e3)^2/ v_std^2) #* sin(t / (2 * 60 * 60 * 2π) )   
winds = (u=u, v=v);

# define variables based on particle equation
t, x, y, c̄_x, c̄_y, lne, Δn, Δφ_p, r_g, C_α, C_φ, g, C_e = vars =  PW4.init_vars();

particle_equations = PW4.particle_equations(u, v, γ=0.88, q=Const_ID.q, input=true, dissipation=true);
@named particle_system = ODESystem(particle_equations);

# define V4 parameters absed on Const NamedTuple:
default_ODE_parameters = Dict(      r_g => r_g0, C_α => Const_Scg.C_alpha, 
                                    C_φ => Const_ID.c_β, C_e => Const_ID.C_e);


# define setting and standard initial conditions
cg_u_local                  = 0.1#FetchRelations.c_g_U_tau(U10, DT) / 1
cg_v_local                  = 0.1#FetchRelations.c_g_U_tau(V10, DT) / 1
lne_local                   = log(FetchRelations.Eⱼ(0.1*abs(U10), DT))

ODE_settings    = PW4.ODESettings(
    Parameters=default_ODE_parameters,
    # define mininum energy threshold
    log_energy_minimum=lne_local,#log(FetchRelations.Eⱼ(0.1, DT)),
    #maximum energy threshold
    log_energy_maximum=log(27),#log(17),  # correcsponds to Hs about 16 m
    saving_step=dt_ODE_save,
    timestep=DT,
    total_time=T,
    callbacks=nothing,
    save_everystep=false)

particle_defaults = ParticleDefaults(ODE_settings.log_energy_minimum, cg_u_local, cg_v_local , 0.0, 0.0)


# % Define grid
grid = TwoDGrid(150e3, 50, 150e3, 50)
periodic_boundary =false;
#bb = mark_boundary(grid)

function TwoDGridMesh(gridl; skip=1)
    gn = TwoDGridNotes(grid)
    gridmesh = [(i, j) for i in gn.x[1:skip:end], j in gn.y[1:skip:end]];
    gridmesh_x = [i for i in gn.x[1:skip:end], j in gn.y[1:skip:end]]
    gridmesh_y = [j for i in gn.x[1:skip:end], j in gn.y[1:skip:end]]
    return (tuples=gridmesh, x=gridmesh_x, y=gridmesh_y)
end

mesh = TwoDGridMesh(grid, skip = 1);
gn = TwoDGridNotes(grid);

# gridmesh = [(i, j) for i in gn.x, j in gn.y];
# mesh.x = [i for i in gn.x, j in gn.y];
# mesh.y = [j for i in gn.x, j in gn.y];
# # mesh.x = mesh.x[1:3:end, 1:3:end];
# # mesh.y = mesh.y[1:3:end, 1:3:end];


# plt.heatmap(gn.x / 1e3, gn.y / 1e3, transpose(u.(mesh.x, mesh.y, 0)))
# # %%

### build model
wave_model = WaveGrowthModels2D.WaveGrowth2D( ;grid = grid, 
                                winds             = winds,
                                ODEsys            = particle_system,
                                ODEvars           = vars,
                                layers            = 1, 
                                ODEsets           = ODE_settings,  # ODE_settings
                                ODEdefaults       = particle_defaults,  # default_ODE_parameters
                                periodic_boundary = periodic_boundary)

#Revise.retry()

### build Simulation
wave_simulation = Simulation(wave_model, Δt=20minutes, stop_time=3hours)
initialize_simulation!(wave_simulation, defaults=wave_model.ODEdefaults)

#init_state_store!(wave_simulation, save_path)


#reset_simulation!(wave_simulation)
## run simulation
#run!(wave_simulation, store=false)

#wave_simulation.stop_time += 1hour

# %%
#pcolormesh in Julia 
# gn = TwoDGridNotes(grid);


# gridmesh = [(i,j) for i in gn.x, j in gn.y];
# mesh.x = [i for i in gn.x, j in gn.y];
# mesh.y = [j for i in gn.x, j in gn.y];
# mesh.x = mesh.x[1:3:end, 1:3:end];
# mesh.y = mesh.y[1:3:end, 1:3:end];

# p1 = plt.heatmap(gn.x / 1e3, gn.y / 1e3, wave_simulation.model.State[:, :, 1])
# plt.quiver!(p1, mesh.x / 1e3, mesh.y / 1e3, quiver=(wave_simulation.model.State[1:3:end, 1:3:end, 2] * 20, wave_simulation.model.State[1:3:end, 1:3:end, 3] * 20), color=:red, scale_unit=:data, label="wind")
# # scale_unit options: :width, :height, :x, :y, :xy, :data

# #plt.heatmap(p1, gn.x / 1e3, gn.y / 1e3, u.(mesh.x, mesh.y, 0))


# # heatmap(wave_simulation.model.State[:, :, 2])
# # heatmap(wave_simulation.model.State[:, :, 3])

# # %%
# # use GLMakie to plot model.State each timestep 


# %%
reset_simulation!(wave_simulation)

n = Observable(1) # for visualization
# Ocean vorticity

model_time = @lift ($n; wave_simulation.model.clock.time)
uo, vo = wave_simulation.model.winds
ocean_wind_u = @lift(uo.(mesh.x, mesh.y, $model_time))
ocean_wind_v = @lift(vo.(mesh.x, mesh.y, $model_time))
#ocean_wind = @lift(  sqrt( vo.(mesh.x, mesh.y, $model_time)^2 + uo.(mesh.x, mesh.y, $model_time)^2) )
#ocean_wind = @lift(sqrt(vo.(mesh.x, mesh.y)^2 + uo.(mesh.x, mesh.y)^2)($model_time) )


wave_energy = @lift ($n; 4 * sqrt.(wave_simulation.model.State[:,:, 1]))
wave_momentum_x = @lift ($n; wave_simulation.model.State[:,:, 2])
wave_momentum_y = @lift ($n; wave_simulation.model.State[:,:, 3])

# for testing
#we = wave_simulation.model.State[:, :, 1]


#compute!(wave_simulation.model.State)
# Ice speed
# ui, vi = ice_model.velocities
# ice_speed = Field(sqrt(ui^2 + vi^2))
# ice_speed_n = @lift ($n; interior(compute!(ice_speed))[:, :, 1])

# Make figure
#x, y, z = nodes((Center, Center, Center), grid)
fig = Figure(resolution=(800, 900))

ax_wind = Axis(fig[1, 1], aspect=1, xlabel="x (km)", ylabel="y (km)", title="Winds")
ax_o = Axis(fig[1, 2], aspect=1, xlabel="x (km)", ylabel="y (km)", title="Hs")
ax_mx = Axis(fig[2, 1], aspect=1, xlabel="x (km)", ylabel="y (km)", title="x momentum")
ax_my = Axis(fig[2, 2], aspect=1, xlabel="x (km)", ylabel="y (km)", title="y momentum")

#hm_i = heatmap!(ax_wind, 1e-3 * x, 1e-3 * y, ice_speed_n)
hm_wind = heatmap!(ax_wind, 1e-3 * gn.x[1:3:end], 1e-3 * gn.y[1:3:end], ocean_wind_u, colormap=:blues)
# colormap options for heatmap 

#quiver!(ax_wind, mesh.x / 1e3, mesh.y / 1e3, quiver=(ocean_wind_u, ocean_wind_v))#, color=:red, scale_unit=:data, label="wind")
#scatter!(ax_wind, vec(gridmesh.* 1e-3), rotations=0, markersize=20, marker='↑')

hm_o = heatmap!(ax_o, 1e-3 * gn.x, 1e-3 * gn.y, wave_energy, colormap=:dense, colorrange=(0.1, 15))
hm_x = heatmap!(ax_mx, 1e-3 * gn.x, 1e-3 * gn.y, wave_momentum_x, colormap=:dense, colorrange=(0, 1))
hm_y = heatmap!(ax_my, 1e-3 * gn.x, 1e-3 * gn.y, wave_momentum_y, colormap=:dense, colorrange=(0, 1))
#colormaps

#colorbar(ax_my)
Colorbar(fig[1, 3], hm_o, label = "Wave energy [m^2]")
Colorbar(fig[2, 3], hm_x, label = "Wave momentum x []")
# sc

#hm_o = heatmap(ax_o, 1e-3 * mesh.x, 1e-3 * mesh.y, ocean_wind_v, colormap=:redblue)

title = @lift ($n; "Coupled ice-ocean simulation after " * prettytime(wave_simulation.model.clock.time))
Label(fig[0, :], title)
display(fig)


#wave_simulation.stop_time += 1hour

N = 30
record(fig, save_path*"viscoelastic_ice.gif", 1:N, framerate=12) do i
    @info "Plotting frame $i of $N..."
    #@time for _ = 1:10
    #run!(wave_simulation, store=false)
    @info wave_simulation.model.clock
    time_step!(wave_simulation.model, wave_simulation.Δt)
    #coupled_time_step!(ocean_simulation, ice_simulation)
    #end
    n[] = 1
end
