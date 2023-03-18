
using DataStructures
using Oceananigans: AbstractDiagnostic, AbstractOutputWriter, fields

import Oceananigans: fields
using Oceananigans.Units

#using storing: StateOrNothing, StateStore, AbstractStore, EmptyStore

mutable struct Simulation{ML, TS, DT, ST, DI, OW, CB}
        model            :: ML
        timestepper      :: TS
        Δt               :: DT
        stop_iteration   :: Float64
        stop_time        :: ST
        wall_time_limit  :: Float64
        diagnostics      :: DI
        output_writers   :: OW
        callbacks        :: CB
        run_wall_time    :: Float64
        running          :: Bool
        initialized      :: Bool
        verbose          :: Bool
        store            :: AbstractStore
        store_itereation :: Int
end

"""
Simulation(model; Δt,
     verbose = true,
     stop_iteration = Inf,
     stop_time = Inf,
     wall_time_limit = Inf)

Construct a `Simulation` for a `model` with time step `Δt`.

Keyword arguments
=================

- `Δt`: Required keyword argument specifying the simulation time step. Can be a `Number`
for constant time steps or a `TimeStepWizard` for adaptive time-stepping.

- `stop_iteration`: Stop the simulation after this many iterations.

- `stop_time`: Stop the simulation once this much model clock time has passed.

- `wall_time_limit`: Stop the simulation if it's been running for longer than this many
           seconds of wall clock time.
"""
function Simulation(model; Δt,
          verbose = true,
          stop_iteration = Inf,
          stop_time = Inf,
          wall_time_limit = Inf)

    if stop_iteration == Inf && stop_time == Inf && wall_time_limit == Inf
    @warn "This simulation will run forever as stop iteration = stop time " *
    "= wall time limit = Inf."
    end

    diagnostics = OrderedDict{Symbol, AbstractDiagnostic}()
    output_writers = OrderedDict{Symbol, AbstractOutputWriter}()
    #callbacks = OrderedDict{Symbol, Callback}()

    # callbacks[:stop_time_exceeded] = Callback(stop_time_exceeded)
    # callbacks[:stop_iteration_exceeded] = Callback(stop_iteration_exceeded)
    # callbacks[:wall_time_limit_exceeded] = Callback(wall_time_limit_exceeded)

    # Check for NaNs in the model's first prognostic field every 100 iterations.
    model_fields = fields(model)
    field_to_check_nans = NamedTuple{keys(model_fields) |> first |> tuple}(first(model_fields) |> tuple)
    #nan_checker = NaNChecker(field_to_check_nans)
    #callbacks[:nan_checker] = Callback(nan_checker, IterationInterval(100))

    # Convert numbers to floating point; otherwise preserve type (eg for DateTime types)
    FT = Float64 #eltype(model.grid)
    Δt = Δt isa Number ? FT(Δt) : Δt
    stop_time = stop_time isa Number ? FT(stop_time) : stop_time

    store = EmptyStore(1)

    return Simulation(model,
            model.timestepper,
            Δt,
            Float64(stop_iteration),
            stop_time,
            Float64(wall_time_limit),
            diagnostics,
            output_writers,
            nothing, #callbacks,
            0.0,
            false,
            false,
            verbose,
            store,
            0)
end


#fields(model::WaveGrowth1D) = (State = model.State, )
