module InputOutput

using ArgParse

export parse_args
export Argsettings

#https://argparsejl.readthedocs.io/en/latest/argparse.html
Argsettings = ArgParseSettings()
@add_arg_table! Argsettings begin
    "--ID"
        help = "ID (or folder) of make -f run_model.make make_collection the model output"
        arg_type = String
    "--T"
        help = "run time in hours"
        arg_type = Float16
    "--DT"
        help = "re-meshing time step in minutes"
        arg_type = Float16

    "--Lx"
        help = "domain length in km"
        arg_type = Float16
    "--Nx"
        help = "# of Nodes"
        arg_type = Int

    "--U10"
        help = "10-meter windspeed amplitude"
        arg_type = Float16

    "--c_beta"
        help = "grow parameters in 1e-2"
        arg_type = Float16
        default = Float16(4.0)
    "--gamma"
        help = "input dissipation coefficient"
        arg_type = Float16
    "--rg"
        help = "input dissipation coefficient"
        arg_type = Float16
    "--periodic"
        help = "flag for periodic boundary condition"
        arg_type = Bool
        action = :store_true
    "--parset"
        help = "set, or group of experiments"
        arg_type = String
    # "Name"
    #     help = "Name of the experiment. if not defined the name is generated"
    #     arg_type = String
end

end
