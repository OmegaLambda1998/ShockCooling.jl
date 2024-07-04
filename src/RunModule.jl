module RunModule

# External Packages
using Supernovae
using Unitful
using UnitfulAstro
const UNITS = [Unitful, UnitfulAstro]

# Internal Packages
include(joinpath(@__DIR__, "ModelModule.jl"))
using .ModelModule
include(joinpath(@__DIR__, "FittingModule.jl"))
using .FittingModule

# Exports
export run_ShockCooling

function modify_supernova!(supernova::Supernova, modifications::Dict{String,Any})
    @debug "Before modification there are $(length(supernova.lightcurve.observations)) observations"

    # Unpack modifications
    min_time = get(modifications, "MIN_TIME", -Inf)
    min_time_unit = uparse(get(modifications, "MIN_TIME_UNIT", "d"), unit_context=UNITS)
    min_time_range = get(modifications, "MIN_TIME_RANGE", [-Inf, Inf]) .* min_time_unit

    max_time = get(modifications, "MAX_TIME", Inf)
    max_time_unit = uparse(get(modifications, "MAX_TIME_UNIT", "d"), unit_context=UNITS)
    max_time_range = get(modifications, "MAX_TIME_RANGE", [-Inf, Inf]) .* max_time_unit

    # Get Min Time
    @debug "Min time option = $min_time"
    if min_time in ["min", "max"]
        observations = [obs for obs in supernova.lightcurve.observations if ((obs.time > min_time_range[1]) & (obs.time < min_time_range[2]))]
        best_obs = observations[1]
        for obs in observations
            # If min_time is min, get the observations with the minimum flux within the min_time_range
            # If min_time is max, get the observations with the maximum flux within the min_time_range
            if ((min_time == "min") & (obs.flux < best_obs.flux)) | ((min_time == "max") & (obs.flux > best_obs.flux))
                best_obs = obs
            end
        end
        min_time = best_obs.time
    else
        min_time = min_time * min_time_unit
    end
    @debug "Min time set to $min_time"

    # Get Max Time
    @debug "Max time option = $max_time"
    if max_time in ["min", "max"]
        observations = [obs for obs in supernova.lightcurve.observations if ((obs.time > max_time_range[1]) & (obs.time < max_time_range[2]))]
        best_obs = observations[1]
        for obs in observations
            # If max_time is min, get the observations with the minimum flux within the max_time_range
            # If max_time is max, get the observations with the maximum flux within the max_time_range
            if ((max_time == "min") & (obs.flux < best_obs.flux)) | ((max_time == "max") & (obs.flux > best_obs.flux))
                best_obs = obs
            end
        end
        max_time = best_obs.time
    else
        max_time = max_time * max_time_unit
    end
    @debug "Max time set to $max_time"

    # Apply modifications
    observations = [obs for obs in supernova.lightcurve.observations if ((obs.time > min_time) & (obs.time < max_time))]
    supernova.lightcurve.observations = observations
    @debug "After modification there are $(length(supernova.lightcurve.observations)) observations"
end

"""
    run_ShockCooling(toml::Dict{String, Any})

Main entrance function for the package

# Arguments
- `toml::Dict{String, Any}`: Input toml file containing all options for the package
"""
function run_ShockCooling(toml::Dict{String,Any})
    config = toml["GLOBAL"]

    # Load supernova data
    @info "Loading in supernova data"
    supernova = run_Supernovae(toml)
    @info "Finished loading $(supernova.name)"

    if "MODIFICATIONS" in keys(toml["DATA"])
        @info "Modifying SN lightcurve"
        modify_supernova!(supernova, toml["DATA"]["MODIFICATIONS"])
    end

    models = toml["MODEL"]
    @info "Loading in $(length(models)) models"
    models = setup_model.(models)
    models = Dict(model.class => model for model in models)
    @info "Finished loading in $(length(models)) models"

    @info "Preparing to fit"
    fitting_opts = toml["FITTING"]

    # Store results
    priors = Dict()
    chains = Dict()
    accept_ratios = Dict()
    logdensities = Dict()
    blobs = Dict()

    # Load chain
    if "CHAIN" in keys(fitting_opts)
        for model_name in keys(fitting_opts["CHAIN"])
            chain_path = abspath(joinpath(config["BASE_PATH"], fitting_opts["CHAIN"][model_name]))
            @info "Loading chain from $chain_path"
            chain, logdensity = load_chain(chain_path)
            chains[model_name] = chain
            logdensities[model_name] = logdensity
        end
    else
        for (class, model) in models
            @info "Fitting $(model.name)"
            prior, chain, accept_ratio, logdensity, blob = run_mcmc(fitting_opts, model, supernova)
            save_chain(joinpath(config["OUTPUT_PATH"], "chain_$(model.name).jld2"), chain, logdensity)
            priors[class] = prior
            chains[class] = chain
            accept_ratios[class] = accept_ratio
            logdensities[class] = logdensity
            blobs[class] = blob
        end
    end

    bestfits = Dict()
    params = Dict()
    for (model_name, chain) in chains
        @info "Analysing $model_name"
        model = models[model_name]
        llhood = likelihood_function(model, supernova)
        bestfit = get_bestfit(chain, logdensities[model_name])
        bestfits[model_name] = bestfit
        param = Dict()
        bestfit_params = []
        for (i, k) in enumerate(sort!(collect(keys(model.parameter_names))))
            min_like = isnan(bestfit[i][3]) ? NaN : round(bestfit[i][3], digits=3) 
            max_like = isnan(bestfit[i][1]) ? NaN : round(bestfit[i][1], digits=3) 
            @info "Marginalised maximum likelihood $(model.parameter_names[k]) = $(round(bestfit[i][2], digits=3))+$max_like/-$min_like $(model.constraints[uppercase(k)][2])"
            @info "Maximum likelihood $(model.parameter_names[k]) = $(round(bestfit[i][4], digits=3)) $(model.constraints[uppercase(k)][2])"
            param[k] = bestfit[i][4] * model.constraints[uppercase(k)][2]
            push!(bestfit_params, bestfit[i][4])
        end
        likelihood = llhood(bestfit_params)
        @info "Bestfit model has likelihood $likelihood"
        params[model_name] = param
    end
end

end
