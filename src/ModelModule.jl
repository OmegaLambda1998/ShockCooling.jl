# Model Module
module ModelModule

# Internal Packages 

# External Packages 
using Supernovae
using Unitful, UnitfulAstro
UNITS = [Unitful, UnitfulAstro]
using Distributions
using LaTeXStrings

# Exports
export Model
export run_model
export setup_model

"""
    abstract type Model

Abstract type for models. All models must be subtypes of this type.
"""
abstract type Model end

"""
    function run_model(model::Model, param::Dict{String,Unitful.Quantity}, supernova::Supernova)

Run a model on a supernova with a given set of parameters. Each model should have its own implementation of this function. If a model is missing this implementation, this will through an error

# Arguments
- `model::Model`: The [`Model`](model) to run.
- `param::Dict{String,Unitful.Quantity}`: A dictionary of parameters to run the model with.
- `supernova::Supernova`: The [`Supernova`](supernova) to run the model on.
"""
function run_model(model::Model, param::Dict{String,Unitful.Quantity}, supernova::Supernova)
    error("generic run_model is being used, no method specified for model $(typeof(model))")
end

"""
    Include all [`Model`](model) files in models directory
"""
const model_path = joinpath(@__DIR__, "models")
for path in readdir(model_path; join = true)
    if isfile(path)
        include(path)
    end
end

"""
    function get_model(model_name::String)

Get a [`Model`](model) by name. This function will return the [`Model`](model) with the given name. If no model with the given name exists, this will throw an error.

# Arguments
- `model_name::String`: The name of the [`Model`](model) to get.
"""
function get_model(model_name::String)
    return getfield(ModelModule, Symbol(model_name))
end

function get_constraints(constraints_dict::Dict{String,Any})
    constraints = Dict{String,Tuple{Distribution,Unitful.FreeUnits}}()
    for param in keys(constraints_dict)
        param_dict = constraints_dict[param]
        param_unit = uparse(param_dict["UNIT"], unit_context = UNITS)
        param_prior = getfield(Distributions, Symbol(param_dict["PRIOR"]))
        param_values = param_dict["VALUES"]
        param_min = get(param_dict, "MIN", -Inf)
        param_max = get(param_dict, "MAX", Inf)
        if param_prior == Normal
            prior = TruncatedNormal(param_values..., param_min, param_max)
        else
            prior = Truncated(param_prior(param_values...), param_min, param_max)
        end
        constraints[param] = (prior, param_unit)
    end
    return constraints
end

function setup_model(model_dict::Dict{String,Any})
    model = get_model(model_dict["NAME"])
    constraints = get_constraints(model_dict["CONSTRAINTS"])
    return model(constraints)
end

end
