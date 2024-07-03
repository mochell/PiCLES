module CartesianGrid

using ...Architectures: AbstractGrid, AbstractGridStatistics, CartesianGrid1D, CartesianGrid2D

#using LinearAlgebra
using StructArrays
using StaticArrays


"""
    CartesianGrid2D( dimx, nx, dimy, ny)

    Generate a cartesians mesh on rectangle `dimx`x `dimy` with `nx` x `ny` points

    - `nx` : indices are in [1:nx]
    - `ny` : indices are in [1:ny]
    - `dimx = xmax - xmin`
    - `dimy = ymax - ymin`
    - `x, y` : node positions
    - `dx, dy` : step size
"""
struct TwoDCartesianGridStatistics <: AbstractGridStatistics
    
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

    area::Float64
    angle_dx::Float64

    function TwoDCartesianGridStatistics(xmin, xmax, Nx::Int, ymin, ymax, Ny::Int ; angle = 0.0)
        dimx = xmax - xmin
        dimy = ymax - ymin

        Ndx = Nx - 1
        Ndy = Ny - 1

        dx = dimx / Ndx
        dy = dimy / Ndy

        area = dx * dy

        return new(Nx, Ny, Ndx, Ndy, xmin, xmax, ymin, ymax, dimx, dimy, dx, dy, area, angle)
    end
end


struct TwoDCartesianGridMesh <: CartesianGrid2D
    data::StructArray{<:Any}
    stats::TwoDCartesianGridStatistics
end
    
function TwoDCartesianGridMesh(grid::TwoDCartesianGridStatistics; mask = nothing)

    x = collect(range(grid.xmin, stop=grid.xmax, step=grid.dx))
    y = collect(range(grid.ymin, stop=grid.ymax, step=grid.dy))

    XX = reshape(repeat(x, inner=length(y)), length(y), length(x))
    YY = reshape(repeat(y, outer=length(x)), length(y), length(x))

    if mask == nothing
        mask = fill(1, size(XX))
    else
        mask = mask
    end

    return StructArray(
        x=XX,
        y=YY,
        mask=mask
        )

end

# initalization
function TwoDCartesianGridMesh(xmin, xmax, Nx::Int, ymin, ymax, Ny::Int, mask=nothing; angle = 0.0)
    GS = TwoDCartesianGridStatistics(xmin, xmax, Nx, ymin, ymax, Ny, angle = angle)
    GMesh = TwoDCartesianGridMesh(GS, mask= mask)
    return TwoDCartesianGridMesh(GMesh, GS)
end

# alternative short hand
TwoDCartesianGridMesh(dimx, nx::Int, dimy, ny::Int; angle =0.0) = TwoDCartesianGridMesh(0.0, dimx, nx, 0.0, dimy, ny; angle=angle)


function ProjetionKernel(G::TwoDCartesianGridMesh)
    if G.stats.angle_dx == 0.0
        M = [
            1/G.stats.dx 0;
            0 1/G.stats.dy
        ]
    else
        cosa = cos(G.stats.angle_dx * pi / 180)
        sina = sin(G.stats.angle_dx * pi / 180)

        M = @SArray [
            cosa/G.stats.dx sina/G.stats.dy;
            sina/G.stats.dx cosa/G.stats.dy
        ]
    end
    return M
end

# %% ADD 1D version here



end # module