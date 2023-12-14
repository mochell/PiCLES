using Printf


# using Distributed
# addprocs(2)
# advance_workers = WorkerPool([2, 3])

push!(LOAD_PATH,   joinpath(pwd(), "code/")       )
# %%

function make_str(name::String, value::Number)
    return name*@sprintf("%i",value)
end


function make_str(name::String, value::Number, digits::Number)
    return name*@sprintf("%0.2f",round(value, digits= digits) )
end


function concat_str(alist)
    str_list =""
    for ai in alist
        str_list *= ai *  "_"
    end
    return str_list[1:end-1]
end


function concat_str_space(alist)
    str_list =""
    for ai in alist
        str_list *= ai *  " "
    end
    return str_list[1:end-1]
end

make_str

# %%

base_parset="1D_gaussian"
parset_name = "Nx-DT"
var1_list  = [20,50, 100, 150, 200, 300]
var2_list   = [5, 10, 15, 20, 30, 40, 50, 60]
periodic = false


# parset_name = "Nx-DT-periodic"
# var1_list  = [10]
# var2_list   = [3, 5, 10, 15, 20, 30, 60]
# var3_list   = [5, 8, 10, 20, 30, 40, 50, 80]
# periodic = true


str_list =[]
arg_list = []
for var1 in var1_list
    for var2 in var2_list
        push!(str_list ,  )

        ID = concat_str([ make_str("Nx", var1), make_str("DT", var2)  ])  # "U10-DT"
        #ID = concat_str([ make_str("rg", var2, 2), make_str("gamma", var1, 2),make_str("Nx", var3)  ])

        #ID = concat_str([ make_str("gamma", var1, 2), make_str("cbeta", var2, 1),make_str("Nx", var3)  ])  # "gamma-c_beta"

        arg_case = [  "--ID", ID,
                    "--Nx",  @sprintf("%i",var1),
                    #"--U10",  @sprintf("%i",var1),
                    #"--c_beta",  @sprintf("%0.2f",var2),
                    #"--gamma",  @sprintf("%0.2f",var1),
                    #"--rg",  @sprintf("%0.2f",var2),
                    "--DT",   @sprintf("%i",var2),
                    "--parset", parset_name]
        if periodic
            push!(arg_case, "--periodic")
        end
        push!(arg_list , arg_case)
        @show arg_case
    end
end



#localARGS= arg_list[5]

#
#aa = ["--ID", "Lx30_DT30_U10", "--Lx", 30, "--DT", 30, "--U", 10]
#push!(aa, "--periodic")
# ARGS = [    "--ID", "test1",
#                 "--Lx", "20",
#                 "--T",  "24",
#                 "--DT", "30"]

#include("./run_fake.jl")#, ARGS=["--foo", "123"])
mkdir(pwd()* "/data/"*base_parset*"/" * parset_name * "/")

# %%
localARGS = arg_list[3]

completed_runs = readdir( pwd()* "/data/"*base_parset*"/" * parset_name * "/" )


@show completed_runs
#@show arg_list
# include( joinpath(pwd(), "parameter_test/") * "/run_fake.jl")
# localARGS

for args in arg_list
    localARGS = args
    #push!(empty!(localARGS), args)
    @show localARGS

    if args[2] in completed_runs
        @printf "runner : case exist, obmitted "
    else
        @printf "runner : run code: "
        #include( joinpath( pwd(), "parameter_test/") * "/run_fake.jl")
        include( joinpath( pwd(), "code/") * "/run_1D_time_varing_forcing.jl")
    end
    @printf "\n runner : next\n"

end

#run(pipeline( `julia -e 'print("hello")' ` , `grep el`))
#localARGS
#run( `julia` )

# %%
# function do_run(args)
#         localARGS = args
#         #@show "before pass", localARGS
#         include("./run_fake.jl")
# end

#@everywhere carrier(argi)    = x -> arg_list[x]

# @everywhere carrier(f)  = x -> f(x)
#
#
# carrier(do_run)(arg_list[1])
#
# carrier( remesh!, State, t)
#
# ParticleCollection = fetch(pmap(carrier( remesh!, State, t) ,advance_workers, arg_list ));
