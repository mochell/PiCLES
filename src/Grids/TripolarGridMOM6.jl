module TripolarGridMOM6

using NCDatasets

using ...Architectures: AbstractGrid, AbstractGridStatistics, TripolarGrid

using StructArrays
using StaticArrays


### ------ basic functions to create the grid -------

"""
    extract_grid_points(x, y, angle_dx, k) -> NamedTuple

    Constructs C-grid Points for tripolar grids, returning the locations of T, U, V, and Q points along with the angle at T points.

    # Inputs:
    - `x`: Array of longitude positions (px, py), where px and py are the dimensions of the grid.
    - `y`: Array of latitude positions (px, py), must have the same dimensions as `x`.
    - `angle_dx`: Array of angles at dx positions (px, py), must have the same dimensions as `x`.
    - `k`: Integer, the stride used to sample the grid points. Determines the spacing between the points in the output. Possible options seam to be k = 2,4,6,8
    - `mask` (optional) mask array given for k =2 on the same grid as x and y. land=1, ocean==0

    # Outputs:
    A named tuple containing the following fields:
    - `Tpoint`: Named tuple with fields `lon` and `lat`, representing the longitude and latitude of T (center) points.
    - `Upoint`: Named tuple with fields `lon` and `lat`, representing the longitude and latitude of U points.
    - `Vpoint`: Named tuple with fields `lon` and `lat`, representing the longitude and latitude of V points.
    - `Qpoint`: Named tuple with fields `lon` and `lat`, representing the longitude and latitude of corner (Q) points.
    - `angle`: Array of angles at T points.
    - `center_index`: Named tuple with fields `x` and `y`, representing the indices of the T points.
    - `corner_index`: Named tuple with fields `x` and `y`, representing the indices of the Q points.
    - `k`: Integer, the stride used to sample the grid points.
    - `khalf`: Integer, half of `k`. Position where to start center points, relative to the corner points.

"""
function extract_grid_points(x, y, angle_dx, k; mask=nothing)

    all_shape = size(x)
    # test if all_shape is equal to size(y)
    if all_shape != size(y)
        error("x and y have different shapes")
    end

    khalf = Int(k / 2) # this is not python, intex start with 1

    center_points_x = khalf+1:k:all_shape[1]
    center_points_y = khalf+1:k:all_shape[2]

    corner_points_x = 1:k:all_shape[1]
    corner_points_y = 1:k:all_shape[2]

    # T point locations
    tlon = x[center_points_x, center_points_y]
    tlat = y[center_points_x, center_points_y]

    # U point locations
    ulon = x[corner_points_x, center_points_y]
    ulat = y[corner_points_x, center_points_y]

    # V point locations
    vlon = x[center_points_x, corner_points_y]
    vlat = y[center_points_x, corner_points_y]

    # Corner point locations
    qlon = x[corner_points_x, corner_points_y]
    qlat = y[corner_points_x, corner_points_y]

    # t-point angle
    angle = angle_dx[center_points_x, center_points_y]

    if mask != nothing
        if k == 2
            mask = mask .== 1 # land is 1
        elseif k in [4, 6, 8]
            mask = mask[1:Int(k / 2):end, 1:Int(k / 2):end] .== 1 # land is 1
        else
            error("k must be 2, 4, 6, or 8")
        end
    else
        mask = nothing
    end

    # Nested tuple output with named nested tuples
    output = (
        Tpoint=(lon=tlon, lat=tlat),
        u=(lon=ulon, lat=ulat),
        v=(lon=vlon, lat=vlat),
        Qpoint=(lon=qlon, lat=qlat),
        angle=angle,
        center_index=(x=center_points_x, y=center_points_y),
        corner_index=(x=corner_points_x, y=corner_points_y),
        k=k,
        khalf=khalf,
        mask=mask
    )
    return output
end


"""
    calculate_distances_k2(area, dx, dy) -> (tarea, dxt, dyt, dxCv, dyCu, dxr, dyr)

    Calculates distances and areas for a grid with a specific stride of k=2. This function only works for the downsampling factor is exactly 2. Its a reproduction of mosaic2mom6.py function.

    # Inputs:
    - `area`: 2D array representing the area of each cell in the original grid.
    - `dx`: 2D array representing the distance between cells in the x-direction in the original grid.
    - `dy`: 2D array representing the distance between cells in the y-direction in the original grid.

    # Outputs:
    - `Tarea`: 2D array of the T cell area, centered at T, aggregated over 2x2 cell blocks.
    - `UVdist`: Named tuple containing:
    - `dxt`: 2D array of the x-distance between U points, centered at T points.
    - `dyt`: 2D array of the y-distance between V oints, centered at T points.
    - `Qdist`: Named tuple containing:
    - `dxCv`: 2D array of the x-distance between Q (corner) points, centered at V points.
    - `dyCu`: 2D array of the y-distance between Q (corner) points, centered at U points.
    - `Tdist`: Named tuple containing:
    - `dxCu`: 2D array of the x-distance between T points, centered at U points, adjusted for the stride.
    - `dyCv`: 2D array of the y-distance between T points, centered at V points, adjusted for the stride and including adjustments for north seam periodicity.

    This function is specifically designed for use with the MOM6 ocean model grid or similar grid structures.
        May also work for other tripolar grids.
"""
function calculate_distances_k2(area, dx, dy)

    k, khalf = 2, 1
    # T cell area (sum of four supergrid cells)
    tarea = area[1:k:end, 1:k:end] + area[1:k:end, 2:k:end] +
            area[2:k:end, 1:k:end] + area[2:k:end, 2:k:end]

    # x-distance between u-points, centered at t
    dxt = dx[1:k:end, khalf+1:k:end] + dx[2:k:end, khalf+1:k:end]

    # y-distance between v-points, centered at t
    dyt = dy[khalf+1:k:end, 1:k:end] + dy[khalf+1:k:end, 2:k:end]

    # x-distance between q-points, centered at v
    dxCv = dx[1:k:end, k+1:k:end] + dx[2:k:end, k+1:k:end]

    # y-distance between q-points, centered at u
    dyCu = dy[k+1:k:end, 1:k:end] + dy[k+1:k:end, 2:k:end]

    ##  ----------------
    # x-distance between t-points, centered at u
    dxr = circshift(dx, (-khalf, 0))
    dxCu = dxr[1:k:end, khalf+1:k:end] + dxr[2:k:end, khalf+1:k:end]

    # y-distance between t-points, centered at v
    dyr = circshift(dy, (0, -khalf))
    # North seam periodicity
    dyr[:, end] = dyr[end:-1:1, end-3]
    dyr[:, end-1] = dyr[end:-1:1, end-2]
    dyCv = dyr[khalf+1:k:end, 1:k:end] + dyr[khalf+1:k:end, 2:k:end]

    # Nested tuple output with named nested tuples
    return (
        Tarea=tarea,
        # x-/y-distance between u-/v-points, centered at t
        UVdist=(dxt=dxt, dyt=dyt),
        # x-/y-distance between q-points, centered at v/u
        Qdist=(dxCv=dxCv, dyCu=dyCu),
        # x-/y-distance between t-points, centered at u/v
        Tdist=(dxCu=dxCu, dyCv=dyCv)
    )
end


# Ndim version:
"""
    calculate_distances(area, dx, dy, k, kp2) -> (tarea, UVdist, Qdist, Tdist)

    Calculates distances and areas for a grid with a specified stride `k`. This function is for aggregating finer grid data into coarser representations

    # Inputs:
    - `area`: 2D array representing the area of each cell in the original grid.
    - `dx`: 2D array representing the distance between cells in the x-direction in the original grid.
    - `dy`: 2D array representing the distance between cells in the y-direction in the original grid.
    - `k`: Integer specifying the stride between grid points, used to determine the output grid's resolution.
    - `kp2`: Integer, typically `k/2`, Position where to start center points, relative to the corner points.

    # Outputs:
    A tuple containing:
    - `Tarea`: 2D array of the T cell area, centered at T, aggregated over 2x2 cell blocks.
    - `UVdist`: Named tuple containing:
    - `dxt`: 2D array of the x-distance between U points, centered at T points.
    - `dyt`: 2D array of the y-distance between V oints, centered at T points.
    - `Qdist`: Named tuple containing:
    - `dxCv`: 2D array of the x-distance between Q (corner) points, centered at V points.
    - `dyCu`: 2D array of the y-distance between Q (corner) points, centered at U points.
    - `Tdist`: Named tuple containing:
    - `dxCu`: 2D array of the x-distance between T points, centered at U points, adjusted for the stride.
    - `dyCv`: 2D array of the y-distance between T points, centered at V points, adjusted for the stride and including adjustments for north seam periodicity.

    This should work for tripolar MOM6 grid with k =2,4,6,8. Others are not tested
"""
function calculate_distances(area, dx, dy, k, khalf)
    # Initialize matrices to store results
    tarea = zeros(size(area)[1] ÷ k, size(area)[2] ÷ k)
    dxt   = zeros(size(dx)[1] ÷ k, size(dx)[2] ÷ k)
    dyt   = zeros(size(dy)[1] ÷ k, size(dy)[2] ÷ k)
    dxCv  = zeros(size(dx)[1] ÷ k, size(dx)[2] ÷ k)
    dyCu  = zeros(size(dy)[1] ÷ k, size(dy)[2] ÷ k)
    dxCu  = zeros(size(dx)[1] ÷ k, size(dx)[2] ÷ k)
    dyCv  = zeros(size(dy)[1] ÷ k, size(dy)[2] ÷ k)
    dxr   = copy(dx)
    dyr   = copy(dy)
    # area_copy  = copy(area)

    # T cell area (sum of four supergrid cells)
    for i in 1:k
        for j in 1:k
            # @info i, j, size(area[i:k:end, j:k:end])
            tarea .+= area[i:k:end, j:k:end]
            # area_copy[i:k:end, j:k:end] .= 0
        end
    end
    # x-distance between u-points, centered at t
    dxt .= sum(dx[j:k:end, khalf+1:k:end] for j in 1:k)

    # y-distance between v-points, centered at t
    dyt .= sum(dy[khalf+1:k:end, i:k:end] for i in 1:k)

    # x-distance between q-points, centered at v
    dxCv .= sum(dx[j:k:end, k+1:k:end] for j in 1:k)

    # y-distance between q-points, centered at u
    dyCu .= sum(dy[k+1:k:end, i:k:end] for i in 1:k)

    # x-distance between t-points, centered at u
    dxr = circshift(dxr, (-khalf, 0))
    for j in 1:k
        dxCu .+= dxr[j:k:end, khalf+1:k:end]
    end

    ##  -------------- upto here ------------
    # y-distance between t-points, centered at v
    dyr = circshift(dyr, (0, -khalf))
    # North seam periodicity
    dyr[:, end] = dyr[end:-1:1, end-3]
    dyr[:, end-1] = dyr[end:-1:1, end-2]
    for i in 1:k
        dyCv .+= dyr[khalf+1:k:end, i:k:end]
    end

    #return tarea, dxt, dyt, dxCv, dyCu, dxr, dyr
    # Nested tuple output with named nested tuples
    return (
        Tarea=tarea,
        # x-/y-distance between u-/v-points, centered at t
        UVdist=(dxt=dxt, dyt=dyt),
        # x-/y-distance between q-points, centered at v/u
        Qdist=(dxCv=dxCv, dyCu=dyCu),
        # x-/y-distance between t-points, centered at u/v
        Tdist=(dxCu=dxCu, dyCv=dyCv),
        # area_copy=area_copy
    )
end


### ------- Structure Defintions ----------
"""
    struct MOM6GridStatistic <: AbstractGrid

    The `MOM6GridStatistic` struct represents the parameters of a MOM6 grid. It is a subtype of the `AbstractGrid` abstract type.

    Fields:
    - `file::String`: Name of the file containing the grid data.
    - `Nx::Int`: Number of grid points in the x-direction.
    - `Ny::Int`: Number of grid points in the y-direction.
    - `Ndx::Int`: Number of spaces in the x-direction.
    - `Ndy::Int`: Number of spaces in the y-direction.
    - `xmin::Float64`: Minimum longitude value.
    - `xmax::Float64`: Maximum longitude value.
    - `ymin::Float64`: Minimum latitude value.
    - `ymax::Float64`: Maximum latitude value.
    - `dimx::Float64`: Longitude range.
    - `dimy::Float64`: Latitude range.
    - `mask_value::Int`: Value used to mask grid points.

"""
struct MOM6GridStatistic <: AbstractGridStatistics

    file::String

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

    mask_value::Int

    function MOM6GridStatistic(Grid; mask_value=1, file="unknown")
        Nx = size(Grid.Tpoint.lon)[1]
        Ny = size(Grid.Tpoint.lon)[2]
        Ndx = Nx - 1
        Ndy = Ny - 1

        xmin = minimum(Grid.Tpoint.lon)
        xmax = maximum(Grid.Tpoint.lon)
        ymin = minimum(Grid.Tpoint.lat)
        ymax = maximum(Grid.Tpoint.lat)

        dimx = xmax - xmin
        dimy = ymax - ymin

        return new(file, Nx, Ny, Ndx, Ndy, xmin, xmax, ymin, ymax, dimx, dimy, mask_value)
    end

end


struct MOM6GridMesh <: TripolarGrid

    data::StructArray{<:Any}  # Adjust the type as necessary
    stats::MOM6GridStatistic

    # initialize with Grid and GridAreaGen objects    
    function MOM6GridMesh(G, GA; mask=nothing, file="unknown")
        # Extracting fields from Grid and GridAreaGen objects

        if mask == nothing
            try
                if G.mask == nothing
                    @info "mask is not provided, using default"
                    # mask = Array{Union{Nothing,Float64}}(1, size(G.Tpoint.lon))
                    mask = trues(size(G.Tpoint.lon))
                else
                    @info "Using provided mask"
                    mask = G.mask
                end
            catch
                # mask = Array{Union{Nothing,Float64}}(1, size(G.Tpoint.lon))
                mask = trues(size(G.Tpoint.lon))
            end
        else
            @info "Using provided mask"
            if size(mask) != size(G.Tpoint.lon)
                error("Mask size must be the same as the grid size")
            end
            mask = mask
        end

        # test mask size

        data = StructArray(
            x=G.Tpoint.lon,
            y=G.Tpoint.lat,
            angle_dx=G.angle,
            dx=GA.Tdist.dxCu,
            dy=GA.Tdist.dyCv,
            area=GA.Tarea,
            mask=mask
        )

        info = MOM6GridStatistic(G; file=file)

        # Create and return the extended struct with Grid and GridAreaGen objects included
        return new(data, info)
    end

    # initialize with Grid and GridAreaGen files
    function MOM6GridMesh(GridFile::String, k::Int; MaskFile::Union{String,Nothing}=nothing)
        # Extracting fields from Grid and GridAreaGen objects

        hgrd     = NCDataset(GridFile)
        x        = hgrd["x"]
        y        = hgrd["y"]
        dx       = hgrd["dx"]
        dy       = hgrd["dy"]
        area     = hgrd["area"]
        angle_dx = hgrd["angle_dx"]


        if MaskFile != nothing
            topo = NCDataset(MaskFile)
            mask = topo["mask"]
        else
            mask = nothing
        end

        Grid     = extract_grid_points(x, y, angle_dx, k; mask=mask)
        GridArea = calculate_distances(area, dx, dy, Grid.k, Grid.khalf)

        close(hgrd)
        if MaskFile != nothing
            close(topo)
        end

        # @info Grid.mask

        return MOM6GridMesh(Grid, GridArea; file=GridFile)#; mask=Grid.mask)
    end

end


## projection Kernel for this grid:

function ProjetionKernel(G::TripolarGrid)
    cosa = cos.(G.data.angle_dx * pi / 180)
    sina = sin.(G.data.angle_dx * pi / 180)
    #@SArray
    proj(cosa, sina, dx, dy, mask) = mask == 1 ? SMatrix{2,2}([
            cosa/dx sina/dy;
            sina/dx cosa/dy]) : nothing

    M = map(proj, cosa, sina, G.data.dx, G.data.dy, G.data.mask)
    return M
end

end 