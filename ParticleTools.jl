module ParticleTools

using DataFrames

using core_1D: ParticleInstance, MarkedParticleInstance

using Plots

function CreateIterationMask(time)
        mask = BitArray([z < 0 for z in diff(time)])
        #breaks = findall(x->x==1, mask)
        lower_bound = insert!(findall(x->x==1, mask), 1, 0)
        lower_bound .+= 1
        upper_bound = insert!(findall(x->x==1, mask), length(findall(x->x==1, mask))+1 , length(time))
        time_mask = time * 0
        for (i1, i2, timem) in zip(lower_bound, upper_bound, range(1,length(upper_bound), step=1))
                #@show i1, i2, time[i1:i2], ParticleData[i1:i2, :]
                time_mask[i1:i2] .= timem
        end
        return time_mask
end

function ParticleToDataframe(a_particle::ParticleInstance)
        return DataFrame(Tables.columntable(a_particle.ODEIntegrator.sol))
end

function ParticleToDataframe(a_particle::MarkedParticleInstance)
        return DataFrame(Tables.columntable(a_particle.Particle.ODEIntegrator.sol))
end

function ParticleToDataframe(Collection::Vector{ParticleInstance})
        DD = []
        for a_particle in Collection
                D         = ParticleToDataframe(a_particle)
                D.mask    = CreateIterationMask(D[!, "timestamp"])
                push!(DD, D[!, [1, 2, 3, 4]] ) # time, x(t), c_x(t), lne(t)
        end
        return DD
end

function ParticleToDataframe(Collection::Vector{MarkedParticleInstance})
        DD = []
        for a_particle in Collection
                D         = ParticleToDataframe(a_particle.Particle)
                D.mask    = CreateIterationMask(D[!, "timestamp"])
                push!(DD, D[!, [1, 2, 3, 4]] ) # time, x(t), c_x(t), lne(t)
        end
        return DD
end

function PlotFailedParticles(Collection, ID, DT, dx)
        energy_list = []
        cg_list = []
        time_list = []
        x_list = []

        time_fail_list =[]
        x_fail_list =[]
        for Fdi in Collection

                di = ParticleTools.ParticleToDataframe(Fdi)

                time_pos = di[:, 1]/DT
                x_pos = di[:, 2]/dx
                lne_pos = di[:, 4]
                #Energy = @. exp.(di[:, 4])
                cg = @. di[:, 3] + 0* time_pos
                Energy = @. di[:, 4] + 0* time_pos

                push!(cg_list, cg)
                push!(energy_list, Energy)
                push!(time_list, time_pos)
                push!(x_list, x_pos)

                push!(time_fail_list, Fdi.time/DT)
                push!(x_fail_list, x_pos[1])

                #scatter!(x, lne, marker=2) |> display
                #scatter!(time, Energy, marker=2) #|> display
                #scat = scatter(time, [ Energy,  cg], layout= (2, 1), marker=2) #|> display
                #plot(alpha_range,  Δ_β.(alpha_range)  , label="Delta_beta")

        end

        p1 = plot(x_list, time_list, title  = "Time") #|> display
        p2 = plot(x_list, energy_list, marker=2, title  = "Energy") #|> display
        p3 = plot(x_list, cg_list, marker=3, title  = "c_g") #|> display

        plot(p1, p2, p3,  layout= (3, 1) , legend= false, size= (600,1200) , ylabel= "time/DT" )
        scatter!( x_fail_list ,time_fail_list, color = "red",marker = 4) #

        plot!(xlabel ="position", plot_title=ID) |> display

end


function plot_cg_version1(Collection)
        plot()

        for Fdi in Collection
                di = ParticleTools.ParticleToDataframe(Fdi)
                #plot!( di[:,1]/DT, di[:,2]/dx, label = "particle trajectory")
                plot!(di[:,3])
                #scatter!([Fdi.time/DT],  [0] ,  marker =2,  color   ="red")
        end
        plot!(legend = :none, title="c_g of initally failed particles", ylabel = "cg m/s", xlabel="position in Particle array") |> display

end


function plot_cg_version2(State_collect::Vector{Any})
        plot()

        for i in 1:length(State_collect)

                statei = State_collect[i]
                energy  =  statei[:,1] #GP.sel(state= 0 ).T # energy
                m_x     =  statei[:,2] #GP.sel(state= 1 ).T # energy
                cg    = energy ./ m_x ./ 2 # c_g
                plot!(cg, color ="black")
        end

        plot!(legend = :none, title="c_g", ylabel = "cg m/s", xlabel="position index x") |> display
end


function plot_cg_version2(State_collect::Array{Float64, 3})
        plot()

        for i in 1:size(State_collect[:,:,:])[1]

                statei = State_collect[i,:,:]
                energy  =  statei[:,1] #GP.sel(state= 0 ).T # energy
                m_x     =  statei[:,2] #GP.sel(state= 1 ).T # energy
                cg    = energy ./ m_x ./ 2 # c_g
                plot!(cg, color ="black")
        end

        plot!(legend = :none, title="c_g", ylabel = "cg m/s", xlabel="position index x") |> display
end



end
