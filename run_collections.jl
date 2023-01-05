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

# parset_name = "U10-DT"
# #var1_list  = [1, 3, 5, 8, 10, 13, 15, 20]
# var1_list  = [20]
# var2_list   = [3, 5, 8, 10, 15, 20, 30, 40, 60]
# var3_list   = [50]#10, 20, 20]
# periodic = false

# parset_name = "Nx-DT"
# var1_list  = [10]
# var2_list   = [3, 5, 10, 15, 20, 30, 60]
# var3_list   = [5, 8, 10, 20, 30, 40, 50, 80]
# periodic = false


# parset_name = "Nx-DT-periodic"
# var1_list  = [10]
# var2_list   = [3, 5, 10, 15, 20, 30, 60]
# var3_list   = [5, 8, 10, 20, 30, 40, 50, 80]
# periodic = true



parset_name = "gamma-c_beta"
var1_list  = round.(range(0.88* 0.2 , 0.88* 1.5; length= 10), digits= 2) # gamma
var2_list  = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6] # c_beta_parameter
var3_list  = [40]#10, 20, 20]
periodic   = false


parset_name = "gamma-c_beta-rg0.85"
var1_list  = round.(range(0.88* 0.2 , 0.88* 1.5; length= 10), digits= 2) # gamma
var2_list  = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6] # c_beta_parameter
var3_list  = [0.85]#10, 20, 20]
periodic   = false


parset_name = "gamma-rg0"
var1_list  = round.(range(0.88* 0.2 , 0.88* 1.5; length= 10), digits= 2) # gamma\ # gamma

sort!(push!(var1_list, 0.88))

var2_list  = round.(range(0.80 , 1/0.80; length= 10), digits= 2)
#[2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6] # c_beta_parameter
var3_list  = [40]#10, 20, 20]
periodic   = false

str_list =[]
arg_list = []
for var1 in var1_list
    for var2 in var2_list
        for var3 in var3_list
            push!(str_list ,  )

            ID = concat_str([ make_str("rg", var2, 2), make_str("gamma", var1, 2),make_str("Nx", var3)  ])
            arg_case = [  "--ID", ID,
                        #"--U10",  @sprintf("%i",var1),
                        #"--c_beta",  @sprintf("%0.2f",var2),
                        "--Nx",  @sprintf("%i",var3),
                        "--gamma",  @sprintf("%0.2f",var1),
                        "--rg",  @sprintf("%0.2f",var2),
                        #"--DT",   @sprintf("%i",var2),
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

completed_runs = readdir( pwd()* "/data/" * parset * "/" )


@show completed_runs
#@show arg_list
# include( joinpath(pwd(), "parameter_test/") * "/run_fake.jl")
# localARGS

for args in arg_list
    localARGS = args
    #push!(empty!(localARGS), args)
    @show localARGS

    if args[2] in completed_runs
        @printf "case exist, obmitted "
    else
        @printf "case not exist, run: "
        #include( joinpath( pwd(), "parameter_test/") * "/run_fake.jl")
        include( joinpath( pwd(), "code/") * "/run_1D.jl")
    end
    @printf "\n next"

end

#run(pipeline( `julia -e 'print("hello")' ` , `grep el`))
#localARGS
#run( `julia` )


# %%

if localARGS[2] in completed_runs
    @printf "case exist, obmitted "
else
    @printf "case not exist, run: "
    #include( joinpath( pwd(), "parameter_test/") * "/run_fake.jl")
    #include( joinpath( pwd(), "code/") * "/run_1D.jl")
end

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
""
