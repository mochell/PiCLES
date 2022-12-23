using Printf


# using Distributed
# addprocs(2)
# advance_workers = WorkerPool([2, 3])

push!(LOAD_PATH,   joinpath(pwd(), "code/")       )
# %%

function make_str(name::String, value::Number)
    return name*@sprintf("%i",value)
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

# %%


# parset_name = "U10-DT"
# #U10_list  = [1, 3, 5, 8, 10, 13, 15, 20]
# U10_list  = [20]
# DT_list   = [3, 5, 8, 10, 15, 20, 30, 40, 60]
# Nx_list   = [50]#10, 20, 20]
# periodic = false

# parset_name = "Nx-DT"
# U10_list  = [10]
# DT_list   = [3, 5, 10, 15, 20, 30, 60]
# Nx_list   = [5, 8, 10, 20, 30, 40, 50, 80]
# periodic = false


# parset_name = "Nx-DT-periodic"
# U10_list  = [10]
# DT_list   = [3, 5, 10, 15, 20, 30, 60]
# Nx_list   = [5, 8, 10, 20, 30, 40, 50, 80]
# periodic = true


parset_name = "U10-c_beta"
U10_list  = [1, 3, 5, 8, 10, 13, 15, 20]
DT_list   = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]
Nx_list   = [40]#10, 20, 20]
periodic = false

str_list =[]
arg_list = []
for u10 in U10_list

    for dt in DT_list
        for Nx in Nx_list
            push!(str_list ,  )

            ID = concat_str([ make_str("Nx", Nx), make_str("DT", dt), make_str("U", u10)  ])
            arg_case = [  "--ID", ID,
                        "--Nx",  @sprintf("%i",Nx),
                        #"--DT",   @sprintf("%i",dt),
                        "--c_beta",  @sprintf("%f",dt),
                        "--U10",  @sprintf("%i",u10),
                        "--parset", parset_name]
            if periodic
                push!(arg_case, "--periodic")
            end
            push!(arg_list , arg_case)
            @show arg_case
        end
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
# %%
localARGS = arg_list[3]

@show arg_list
# include( joinpath(pwd(), "parameter_test/") * "/run_fake.jl")
# localARGS

for args in arg_list
    localARGS = args
    #push!(empty!(localARGS), args)
    @show localARGS
    #include( joinpath( pwd(), "parameter_test/") * "/run_fake.jl")
    include( joinpath( pwd(), "code/") * "/run_1D.jl")
end

#run(pipeline( `julia -e 'print("hello")' ` , `grep el`))
localARGS
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
