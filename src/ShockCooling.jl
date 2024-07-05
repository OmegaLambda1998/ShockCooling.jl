module ShockCooling

# External packages
using TOML
using BetterInputFiles
using ArgParse
using OrderedCollections

# Internal Packages
include("RunModule.jl")
using .RunModule

# Exports
export main

"""
    get_args()

Helper function to get the ARGS passed to julia.
"""
function get_args()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--verbose", "-v"
        help = "Increase level of logging verbosity"
        action = :store_true
        "input"
        help = "Path to .toml file"
        required = true
    end
    return parse_args(s)
end

"""
    main()

Read the args, prepare the input TOML and run the actual package functionality.
"""
function main()
    args = get_args()
    toml_path = args["input"]
    verbose = args["verbose"]
    main(toml_path, verbose)
end

function main(toml_path::AbstractString, verbose::Bool)
    paths = OrderedDict(
        "data_path" => ("base_path", "Data"),
        "filter_path" => ("base_path", "Filter"),
    )
    toml = setup_input(toml_path, verbose; paths = paths)
    run_ShockCooling(toml)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
