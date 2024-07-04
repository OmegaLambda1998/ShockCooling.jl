module ShockCooling

# External packages
using TOML
using BetterInputFiles
using ArgParse
using StatProfilerHTML
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
        "--profile", "-p"
        help = "Run profiler"
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
    profile = args["profile"]
    main(toml_path, verbose, profile)
end

function main(toml_path::AbstractString, verbose::Bool, profile::Bool)
    paths = OrderedDict(
        "data_path" => ("base_path", "Data"),
        "filter_path" => ("base_path", "Filter"),
    )
    toml = setup_input(toml_path, verbose; paths = paths)
    if profile
        @warn "Running everything once to precompile before profiling"
        run_ShockCooling(toml)
        @profilehtml run_ShockCooling(toml)
    else
        run_ShockCooling(toml)
    end

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
