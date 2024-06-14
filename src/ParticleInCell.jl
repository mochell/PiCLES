module ParticleInCell

using Statistics

using SharedArrays


# %%
using ..ParticleMesh
using ...Architectures: AbstractParticleInstance
# Tolerance for comparison of real numbers: set it here!

# Set parameters
# eta_min = 0.0
# eta_max = 20.0
#
# num_cells = 22
# delta_t = 0.1
#
# grid2d = TwoDGrid(eta_max,num_cells, eta_max ,num_cells)
# gridnotes = TwoDGridNotes(grid2d)

export push_to_grid!
# %%

function get_i_and_w(zp_normed::Float64)

    # floor points x
    ip_base = floor(zp_normed) # ceil position
    #ip_floor = Int(ip_base)
    # if sign(ip_base) == -1.0
    #     ip_floor = Int(ip_base)
    # else
    #     ip_floor = Int(ip_base+ 1)
    # end
    ip_floor = Int(ip_base+ 1)

    dxp_ceil = zp_normed - ip_base   # floor weight
    dxp_floor = 1.0 - dxp_ceil #ip_base+1 - xp_normed # ceil weight
    return [ip_floor , ip_floor+1] , [dxp_floor, dxp_ceil]
end

"""
compute_weights_and_index(g_pars::TwoDGrid, xp::Float64, yp:: Float64 )
returns indexes and weights for in 2D for single x,y point
"""
function compute_weights_and_index(g_pars::TwoDGrid, xp::Float64, yp:: Float64 )
    """
    2d wrapper for 1d function
    """

    xp_normed = (xp - g_pars.xmin) / g_pars.dx # multiples of grid spacing
    yp_normed = (yp - g_pars.ymin) / g_pars.dy # multiples of grid spacing

    xi, xw = get_i_and_w(xp_normed)
    yi, yw = get_i_and_w(yp_normed)

    idx  = [ (xi[1], yi[1]) , (xi[2], yi[1]) , (xi[1], yi[2]), (xi[2], yi[2]) ]
    wtx  = [ (xw[1], yw[1]) , (xw[2], yw[1]) , (xw[1], yw[2]), (xw[2], yw[2]) ]

    return idx, wtx
end

"""
compute_weights_and_index(g_pars::OneDGrid, xp::Float64 )
returns indexes and weights for in 2D for single x point
"""
function compute_weights_and_index(g_pars::OneDGrid, xp::Float64 )

    xp_normed = (xp - g_pars.xmin) / g_pars.dx # multiples of grid spacing
    xi, xw = get_i_and_w(xp_normed)

    idx  = [ (xi[1]) , (xi[2]) ]
    wtx  = [ (xw[1]) , (xw[2]) ]

    return idx, wtx
end

#indexes, weights = compute_weights_and_index(grid2d, 10.0, 9.9 )

"""
compute_weights_and_index(g_pars::TwoDGrid, xp::Float64, yp:: Float64 )
returns indexes and weights for in 2D for vectors
"""
function compute_weights_and_index(grid::TwoDGrid, xp::Vector{Float64}, yp:: Vector{Float64} )
    """
    returns:
    weight list        2 N X 1 vector of tuples with indices and weights
    """
    #return reduce(vcat,[compute_weights_and_index_1d(g_pars, xi) for xi in xp])
    index_list, weight_list  = [], []
    for (xi,yi) in zip(xp, yp)
        dd, ff = compute_weights_and_index(grid, xi, yi)
        push!(index_list, dd)
        push!(weight_list,ff)
    end

    #return reduce(vcat, index_list ), reduce(vcat, weight_list)
    return index_list, weight_list
end




"""
compute_weights_and_index(g_pars::TwoDGrid, xp::Float64 )
returns indexes and weights for in 1D for vectors
"""
function compute_weights_and_index(grid::OneDGrid, xp::Vector{Float64})
    """
    returns:
    weight list        2 N X 1 vector of tuples with indices and weights
    """
    #return reduce(vcat,[compute_weights_and_index_1d(g_pars, xi) for xi in xp])
    index_list, weight_list  = [], []
    for xi in xp
        dd, ff = compute_weights_and_index(grid, xi)
        push!(index_list, dd)
        push!(weight_list,ff)
    end

    #return reduce(vcat, index_list ), reduce(vcat, weight_list)
    return index_list, weight_list
end


#index_positions, weights = compute_weights_and_index(grid2d, [0.1, 11.1, 0.1, 5.234, -4.3], [0.1, 0.5, 2.9, 9.99, -3.2])

# function wrap_index!(pos::Int, N::Int)
#     if pos < 0
#         pos = pos  % Int(N)
#         #pos =  N
#     elseif pos > N
#         pos = pos % Int(N)
#         #pos = N
#     elseif pos == 0
#         pos = pos + Int(N)
#     end
#     return pos
# end

function wrap_index!(pos::Int, N::Int)
        pos = pos % Int(N)

        if pos < 0
                pos = pos + Int(N)
                #pos =  N
        elseif pos == 0
                pos = pos + Int(N)
        end
    return pos
end


#### merging rules: ####

function merge!(grid_point::Vector{Float64}, charge::Vector{Float64}; verbose=false)

    verbose ? (@info "grid_point = $grid_point, charge = $charge") : nothing

    m_grid = grid_point[2]
    m_charge = charge[2]
    if (sign(m_grid) == sign(m_charge)) | (m_grid == 0)
        verbose ? (@info "add, same sign, or grid is zero") : nothing
        grid_point += charge
    elseif (sign(m_grid) != sign(m_charge)) & (sign(m_grid) * m_grid > sign(m_charge) * m_charge)
        verbose ? (@info "forget charge, State is larger") : nothing

    elseif (sign(m_grid) != sign(m_charge)) & (sign(m_grid) * m_grid <= sign(m_charge) * m_charge)
        verbose ? (@info "overwrite, charge is >=") : nothing
        grid_point = charge
    end
    return grid_point
end

function merge!(grid_point::Float64, charge::Float64; verbose=false)

    verbose ? (@info "grid_point = $grid_point, charge = $charge") : nothing

    m_grid = grid_point
    m_charge = charge
    if (sign(m_grid) == sign(m_charge)) | (m_grid == 0)
        verbose ? (@info "add, same sign, or grid is zero") : nothing
        m_grid += m_charge
    elseif (sign(m_grid) != sign(m_charge)) & (sign(m_grid) * m_grid > sign(m_charge) * m_charge)
        verbose ? (@info "forget charge, State is larger") : nothing

    elseif (sign(m_grid) != sign(m_charge)) & (sign(m_grid) * m_grid <= sign(m_charge) * m_charge)
        verbose ? (@info "overwrite, charge is >=") : nothing
        m_grid = m_charge
    end
    return m_grid
end



#define merging operator
(⊓)(g::Vector{Float64}, c::Vector{Float64}) = merge!(g, c, verbose=false)
(⊓)(g::Float64, c::Float64) = merge!(g, c, verbose=true)



#wrap_index!.(collect(range(0, gridnotes.Nx+1, step=1)), gridnotes.Nx)
#wrap_index!.(collect(range(-2, gridnotes.Ny*2, step=1)), gridnotes.Ny)

#charges_grid = zeros(grid2d.Nx, grid2d.Ny)


## 2D versions
function push_to_grid!(grid::Matrix{Float64},
                            particlesAtNode::Array{Array{Array{Any,1},1},1},
                            PI::AbstractParticleInstance,
                            charge::Float64,
                            index_pos::Tuple{Int, Int},
                            weights::Tuple{Float64, Float64},
                            Nx::Int,  Ny::Int )

        grid[ wrap_index!(index_pos[1], Nx) , wrap_index!(index_pos[2], Ny) ] += weights[1] * weights[2] * charge
        #grid[ index_pos[1] , index_pos[2] ] += weights[1] * weights[2] * charge

end

function push_to_grid!(grid::SharedArray{Float64, 3},
                            particlesAtNode::Array{Array{Array{Any,1},1},1},
                            PI::AbstractParticleInstance,
                            charge::Vector{Float64},
                            index_pos::Tuple{Int, Int},
                            weights::Tuple{Float64, Float64},
                            Nx::Int,  Ny::Int,
                            periodic::Bool = true)
    if periodic
        grid[ wrap_index!(index_pos[1], Nx) , wrap_index!(index_pos[2], Ny), : ] += weights[1] * weights[2] * charge
    else
        if sum(test_domain(index_pos, Nx, Ny)) != 2
            return
        else
            push!(particlesAtNode[index_pos[1]][index_pos[2]], [weights[1] * weights[2], charge[4], PI])
            grid[ index_pos[1] , index_pos[2], 1:3 ] += weights[1] * weights[2] * charge[1:3]
        end
    end

end


"""
test_domain(index_pos, Nx, Ny)
test if the index position is within the domain
"""
function test_domain(index_pos, Nx, Ny)
    return (index_pos[1] > 0 && index_pos[1] <= Nx), (index_pos[2] > 0 && index_pos[2] <= Ny)
end


#push_to_grid!(charges_grid,1.0 , index_positions[1][1], weights[1][1], grid2d.Nx , grid2d.Ny )


# wrapping over vectors of charges, index positions and weights
function push_to_grid!(grid::SharedArray{Float64,3},
                            particlesAtNode::Array{Array{Array{Any,1},1},1},
                            PI::AbstractParticleInstance,
                            charge::Vector{Float64},
                            index_pos::Vector{Tuple{Int, Int}},
                            weights::Vector{Tuple{Float64, Float64}},
                            Nx::Int, Ny::Int,
                            periodic::Bool=true)
    #@info "this is version D"
    for (i, w) in zip(index_pos, weights)
        push_to_grid!(grid, particlesAtNode, PI, charge , i, w , Nx, Ny, periodic)
    end
end

#push_to_grid!(charges_grid, 1.0 , index_positions[3], weights[3] , grid2d.Nx , grid2d.Ny )




###### 1D versions ####

# wrapper over 1D Vecors of chanegs and (nested) index positions and weights
function push_to_grid!(grid::SharedMatrix{Float64},
    particlesAtNode::Array{Array{Array{Any,1},1},1},
    PI::AbstractParticleInstance,
    charge::Vector{Float64},
    index_pos::Vector{Any},
    weights::Vector{Any},
    Nx::Int,
    periodic::Bool=true)
    #@info "this is version C"
    for (im, wm, c) in zip(index_pos, weights, charge)
        for (i, w) in zip(im, wm)
            push_to_grid!(grid, particlesAtNode, PI, c, i, w, Nx, periodic)
        end
    end
end



# multiple index positions and 1 charge
function push_to_grid!(grid::MM,
                            particlesAtNode::Array{Array{Array{Any,1},1},1},
                            PI::AbstractParticleInstance,
                            charge::Vector{Float64},
                            index_pos::Vector{Int},
                            weights::Vector{Float64},
                            Nx::Int,
                            periodic::Bool = true) where MM <: Union{SharedMatrix{Float64}, Matrix{Float64}}
    #@info "this is version B"
    if periodic
        for (im, wm) in zip(index_pos, weights)
            #grid[ wrap_index!(im, Nx), :] ⊓ wm * charge
            # old version
            #grid[wrap_index!(im, Nx), :] += wm * charge
            # merge rule version
            grid[wrap_index!(im, Nx), :] = merge!(grid[wrap_index!(im, Nx), :], wm * charge)
        end
    else
        for (im, wm) in zip(index_pos, weights)
            if (im <= Nx) & (im > 0)
                # old version
                #grid[im, : ] += wm * charge

                # merge rule version
                grid[im, :] = merge!(grid[im, :], wm * charge)
                #grid[im, :] = grid[im, :] ⊓ wm * charge
            end
        end
    end
end


# only 1 index position and 1 charge
function push_to_grid!(grid::MM,
                            particlesAtNode::Array{Array{Array{Any,1},1},1},
                            PI::AbstractParticleInstance,
                            charge::Float64,
                            index_pos::Int,
                            weights::Float64,
                            Nx::Int,
                            periodic::Bool = true) where MM <: Union{SharedMatrix{Float64}, Matrix{Float64}}
        #@info "this is version A"
        if periodic
            #grid[ wrap_index!(index_pos[1], Nx) ] += weights * charge
            #grid[ index_pos[1] , index_pos[2] ] += weights[1] * weights[2] * charge

            #grid[wrap_index!(index_pos[1], Nx)] = grid[wrap_index!(index_pos[1], Nx)] ⊓ weights * charge
            grid[wrap_index!(index_pos[1], Nx)] = merge!(grid[wrap_index!(index_pos[1], Nx)], weights * charge)
        else
            if (index_pos <= Nx) & (index_pos > 0)
                #grid[index_pos] ⊓ weights * charge
                #grid[index_pos] += weights * charge
                grid[index_pos] = merge!(grid[index_pos], weights * charge)
            end
        end
end


# function push_to_grid!(grid::MM,
#                             charge::Float64,
#                             index_pos::Int,
#                             weights::Float64,
#                             Nx::Int,
#                             periodic::Bool = true) where MM <: Union{SharedMatrix{Float64}, Matrix{Float64}}
#         if periodic
#             grid[ wrap_index!(index_pos[1], Nx) ] += weights * charge
#             #grid[ index_pos[1] , index_pos[2] ] += weights[1] * weights[2] * charge
#         else
#             if (im <= Nx) & (im > 0)
#                 grid[im, : ] += wm * charge
#             end
#         end
# end




# end of module
end
# %%
# charges_grid = zeros(grid2d.Nx, grid2d.Ny)
# x_list = [5.0, 1.234, 0.5]#[0.1, 0.1, 0.5]#, 5.234, -4.3]
# y_list = [6.0, 2.34, 0.1]#[0.5, 0.5, 0.5]#, 9.99, -3.2]
# gridnotes.dx
#
#
# x_list = gridnotes.x .+ gridnotes.dx * 0.5
# y_list = gridnotes.y .+ gridnotes.dy* 0.9999999 ##[0.1, 0.1, 0.5]#, 5.234, -4.3]
# #y_list = [6.0, 2.34, 0.1]#[0.5, 0.5, 0.5]#, 9.99, -3.2]
# charges = gridnotes.y *0 .+1.0#, 1.0]#, 1.0, 1.0, 1.0]
#
# index_positions, weights = compute_weights_and_index(grid2d,x_list, y_list )
# push_to_grid!(charges_grid, charges , index_positions, weights  , grid2d.Nx , grid2d.Ny )

# %%
# n_particles = 100
# grid2d = mesh.TwoDGrid(eta_max,n_particles, eta_max ,50)
# gridnotes = mesh.TwoDGridNotes(grid2d)
#
# xmesh = [ x for x in gridnotes.x, y in gridnotes.y] # .+ 0.3
# ymesh = [ y for x in gridnotes.x, y in gridnotes.y] # .+0.5
#
# xp = xmesh .+ grid2d.dx* 0.1
# yp = ymesh .+ grid2d.dy* 0.1
#
# charge = xmesh .* 0
# charge[48:75, 5:30] .= sin.(gridnotes.x[48:75])
# #charge[20:40, 5:20] .= -1.0
#
# heatmap(gridnotes.x, gridnotes.y, charge')

# %%

# charges_grid = xmesh .* 0
#
# x_list = reshape(xp , (1, grid2d.Nx * grid2d.Ny) )
# y_list = reshape(yp , (1, grid2d.Nx * grid2d.Ny) )
# charge_list = reshape(charge , (1, grid2d.Nx * grid2d.Ny) )
#
# x_list      = x_list[charge_list .!= 0]
# y_list      = y_list[charge_list .!= 0]
# charge_list = charge_list[charge_list .!= 0]
#
# index_positions, weights = compute_weights_and_index(grid2d, x_list, y_list)
# #index_positions, weights = compute_weights_and_index(grid2d,dropdims(x_list, dims= 1),dropdims(y_list, dims= 1))
#
# push_to_grid!(charges_grid, charge_list , index_positions, weights  , grid2d.Nx , grid2d.Ny )
#
# #charge = charges_grid
#
# #contour(gridnotes.x, gridnotes.y, charges_grid)
# heatmap(gridnotes.x, gridnotes.y, charges_grid')



# %%
#dropdims(reshape(xp , (1, grid2d.Nx * grid2d.Ny) ), dims=1 )
# %%

# p_collect=[]
#
#
# cg0 = 0.12332 *1
#
#
# i_charge =copy(charge)
# charge_sum = sum(charge)
#
# i_charge
# for ti in 1:1:200
#
#     # xp = xmesh .+ grid2d.dx* cg0 .+ randn(100, 50) *0.01 + i_charge
#     # yp = ymesh .- grid2d.dy* cg0 .+ randn(100, 50) *0.01
#     xp = xmesh .+ i_charge*0.2
#     yp = ymesh .+ i_charge*0.05
#
#     #xp                       = grid1d.x .+ cg0 .*ip_charges  #.* (0.9 .+ 0.05 *(1 .+  rand(grid1d.Nx)))
#     #xp                       = grid1d.x .+ cg0 .* (0.9 .+ 0.05 *(1 .+  rand(grid1d.Nx)))
#     #xp                       = grid1d.x .+  0.1 * (.-0.5 .+  rand(grid1d.Nx))
#     #xp                       = xp .+  0.1 * (.-0.2 .+  rand(grid1d.Nx))
#
#     x_list = dropdims(reshape(xp , (1, grid2d.Nx * grid2d.Ny) ), dims=1 )
#     y_list = dropdims(reshape(yp , (1, grid2d.Nx * grid2d.Ny) ) , dims= 1)
#     charge_list = dropdims(reshape(i_charge , (1, grid2d.Nx * grid2d.Ny) ), dims = 1)
#
#     x_list      = x_list[charge_list .!= 0]
#     y_list      = y_list[charge_list .!= 0]
#     charge_list = charge_list[charge_list .!= 0]
#
#     index_positions, weights = compute_weights_and_index(grid2d, x_list, y_list)
#
#     charges_grid = xmesh .* 0
#     push_to_grid!(charges_grid, charge_list , index_positions, weights  , grid2d.Nx , grid2d.Ny )
#
#     @show charge_sum - Float64( sum(charges_grid) )
#
#     # propagate
#     i_charge = charges_grid
#     push!(p_collect, charges_grid)
# end

# %%
# anim = @animate for i in 1:1:length(p_collect)
#
#     #scatter(x_collect[i], p_collect[i], label ="new" )
#     heatmap(gridnotes.x, gridnotes.y, p_collect[i]', clim=(-1.0, 1.0))
#     # if i > 1
#     #     scatter!(x_collect[i-1], p_collect[i-1], color= "gray", alpha = 0.4, label="old")
#     # end
#     title!( "sum = $( round(sum(p_collect[i]) )) ")
#
#     #ylims!( 0, 1.1)
#     #xlims!(eta_min, eta_max,)
# end
# #round(sum(p_collect[1]) )
#
#
# gif(anim, "PIC_2d_test_fps30.gif", fps = 10)
