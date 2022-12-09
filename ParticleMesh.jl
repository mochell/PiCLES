module ParticleMesh

abstract type AbstractGrid end

export OneDGrid, OneDGridNotes, TwoDGrid, TwoDGridNotes


"""
    TwoDGrid( dimx, nx, dimy, ny)

Generate a cartesians mesh on rectangle `dimx`x `dimy` with `nx` x `ny` points

- `nx` : indices are in [1:nx]
- `ny` : indices are in [1:ny]
- `dimx = xmax - xmin`
- `dimy = ymax - ymin`
- `x, y` : node positions
- `dx, dy` : step size
"""
struct TwoDGrid <: AbstractGrid
    Nx::Int
    Ny::Int
    Ndx::Int
    Ndy::Int

    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64

    dimx::Float64
    dimy::Float64

    dx::Float64
    dy::Float64

    function TwoDGrid(xmin, xmax, Nx::Int, ymin, ymax, Ny::Int)
        dimx = xmax - xmin
        dimy = ymax - ymin

        Ndx = Nx -1
        Ndy = Ny -1

        dx = dimx / Ndx
        dy = dimy / Ndy

        return new(Nx, Ny, Ndx, Ndy, xmin, xmax, ymin, ymax, dimx, dimy, dx, dy)
    end
end

"""
    TwoDGrid( xmin, xmax, nx, ymin, ymax, ny )

Simple structure to store mesh data from 2 dimensions
"""
TwoDGrid(dimx, nx::Int, dimy, ny::Int) = TwoDGrid(0.0, dimx, nx, 0.0, dimy, ny)



struct TwoDGridNotes <: AbstractGrid
    Nx::Int
    Ny::Int
    Ndx::Int
    Ndy::Int

    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64

    dimx::Float64
    dimy::Float64
    dx::Float64
    dy::Float64

    x::Vector{Float64}
    y::Vector{Float64}

    function TwoDGridNotes( gridi::TwoDGrid )
        # x = collect(LinRange(0, gridi.dimx, gridi.Nx))
        # y = collect(LinRange(0, gridi.dimy, gridi.Ny))
        x = collect(range(gridi.xmin, stop = gridi.xmax, step = gridi.dx))
        y = collect(range(gridi.ymin, stop = gridi.ymax, step = gridi.dy))

        return new(gridi.Nx, gridi.Ny, gridi.Ndx, gridi.Ndy,
                    gridi.xmin, gridi.xmax, gridi.ymin, gridi.ymax,
                    gridi.dimx, gridi.dimy, gridi.dx, gridi.dy,x, y)
    end
end

struct OneDGrid <: AbstractGrid
    Nx::Int # number grid points
    Ndx::Int # number cells
    xmin::Float64
    xmax::Float64
    dimx::Float64
    dx::Float64

    function OneDGrid(xmin, xmax, Nx::Int)
        dimx = xmax - xmin
        Ndx = Nx -1
        dx = dimx / Ndx

        return new(Nx, Ndx, xmin, xmax, dimx, dx)#, x)
    end
end

struct OneDGridNotes <: AbstractGrid
    Nx::Int # number grid points
    Ndx::Int # number of cells
    xmin::Float64
    xmax::Float64
    dimx::Float64
    dx::Float64
    x::Vector{Float64}

    function OneDGridNotes(grid::OneDGrid )
        x = collect(LinRange(0, grid.dimx, grid.Nx))
        return new(grid.Nx, grid.Ndx, grid.xmin, grid.xmax, grid.dimx, grid.dx, x)
    end
end

export get_x

"""
    get_x( mesh, i )

Get position
"""
function get_x(m::OneDGrid, i)
    return m.xmin + (i - 1) * m.dx
end

"""
    get_cell_and_offset( mesh, x )

Get cell and offset

We compute the cell indices where the particle is and its relative
normalized position inside the cell

"""
function get_cell_and_offset(m::OneDGrid, x)
    cell = floor(Int64, ((x - m.xmin) / m.Lx) * m.nx) + 1
    offset = (x - get_x(m, cell)) / dx

    return cell, offset
end

end
