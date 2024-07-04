# Fitting Module
module FittingModule

# Internal Packages 
using ..ModelModule

# External Packages 
using Unitful, UnitfulAstro
const UNITS = [Unitful, UnitfulAstro]
using KissMCMC
using Distributions
using Supernovae
using JLD2
using StatsBase
using LinearAlgebra
using Interpolations

# Exports
export run_mcmc
export save_chain, load_chain
export likelihood_function
export get_bestfit

"""
    get_prior(model::Model, numwalkers::Int64)

Generate a set of initial parameters for the MCMC walkers. The parameters are drawn from the prior distributions specified in the model.

# Arguments
- `model::Model`: The [`Model`](model) to generate priors for
- `numwalkers::Int64`: The number of walkers to generate
"""
function get_prior(model::Model, numwalkers::Int64)
    ks = sort!(collect(keys(model.constraints)))
    x0 = [[rand(model.constraints[k][1]) for k in ks] for _ in 1:numwalkers]
    return x0
end

"""
    prior_constraint(model::Model, param::Dict)

Calculate the prior probability of a set of parameters. This is the sum of the log probabilities of each parameter drawn from the prior distributions.

# Arguments
- `model::Model`: The [`Model`](model) to generate the prior probability for
- `param::Dict`: The parameters to generate the prior probability
"""
function prior_constraint(model::Model, param::Dict)
    return sum([logpdf(model.constraints[k][1], ustrip(param[k])) for k in keys(model.constraints)])
end

"""
    likelihood_function(model::Model, supernova::Supernova)

Generate a function that calculates the likelihood of a set of parameters given a supernova. This function is used by the MCMC sampler.

# Arguments
- `model::Model`: The [`Model`](model) to generate the likelihood function for
- `supernova::Supernova`: The [`Supernova`](supernova) to generate the likelihood function for
"""
function likelihood_function(model::Model, supernova::Supernova)
    function llhood(param)
        ks = sort!(collect(keys(model.constraints)))
        param = Dict(k => param[i] * model.constraints[k][2] for (i, k) in enumerate(ks))
        prior = prior_constraint(model, param)
        # If outside physical bounds, don't bother evaluating the probability
        if !isfinite(prior)
            return prior
        end
        m_absmag = run_model(model, param, supernova)
        m_mag = absmag_to_mag.(m_absmag, supernova.redshift)
        m_flux = mag_to_flux.(m_mag, supernova.zeropoint)
        chi = sum(((m_flux .- get(supernova, "flux")) ./ get(supernova, "flux_err")) .^ 2)
        r = length(m_flux) - length(ks)
        r_chi = -0.5 * chi / r
        ll = prior + r_chi
        return ll
    end
    return llhood
end

function save_chain(path::AbstractString, chain, logdensity)
    save(path, Dict("CHAIN" => chain, "LOGDENSITY" => logdensity))
end

function load_chain(path::AbstractString)
    d = load(path)
    return d["CHAIN"], d["LOGDENSITY"]
end

function run_mcmc(config::Dict{String,Any}, model::Model, supernova::Supernova)
    @debug "Running with $(Threads.nthreads()) threads"
    # Important config
    numwalkers = Int64(config["NUMWALKERS"])
    @debug "Numwalkers: $numwalkers"
    thinning = Int64(get(config, "THINNING", 1))
    @debug "Thinning: $thinning"
    burnin = Int64(config["BURNIN"])
    @debug "Burnin: $burnin"
    iterations = Int64(config["ITERATIONS"])
    @debug "Iterations: $iterations"
    llhood = likelihood_function(model, supernova)

    # Unimportant config
    use_progress_meter = get(config, "PROGRESS_METER", true)

    @info "Generating priors"
    x0 = get_prior(model, numwalkers)

    @info "Running MCMC"
    chain, accept_ratio, logdensities, blob = emcee(llhood, x0; niter=iterations, nburnin=burnin, nthin=thinning, use_progress_meter=use_progress_meter)
    @info "MCMC Finished"

    chain, accept_ratio, logdensities, blob = squash_walkers(chain, accept_ratio, logdensities, blob)
    @info "MCMC has accept ratio of $accept_ratio"

    return x0, chain, accept_ratio, logdensities, blob
end


function get_smoothed_histogram(chain::Vector{Float64}, logdens::Vector{Float64}, target::Float64=0.6827)
    h = fit(StatsBase.Histogram, chain; nbins=floor(Int64, sqrt(length(chain)) / 10))
    h = normalize(h, mode=:pdf)
    edges = collect(h.edges[1])
    hist = h.weights
    edge_centers = @. 0.5 * (edges[2:end] + edges[1:end-1])
    xs = LinRange(edge_centers[1], edge_centers[end], 10000)
    ys = linear_interpolation(edge_centers, hist)(xs)
    cs = cumsum(ys)
    cs = cs ./ maximum(cs)
    return xs, ys, cs
end

function get_bestfit(chain::Vector{Float64}, logdens::Vector{Float64}, i::Int64, target::Float64=0.6827)
    xs, ys, cs = get_smoothed_histogram(chain, logdens)
    n_pad = 1000
    x_start = xs[1] * ones(n_pad, 1)
    x_end = xs[end] * ones(n_pad, 1)
    y_start = LinRange(0, ys[1], n_pad)
    y_end = LinRange(ys[end], 0, n_pad)
    xs = vcat(x_start, xs, x_end)
    ys = vcat(y_start, ys, y_end)
    cs = cumsum(ys)
    cs = cs ./ maximum(cs)
    startIndex = argmax(ys)
    maxVal = ys[startIndex]
    minVal = 0
    mid = 0.5 * (maxVal + minVal)
    l1 = reverse(ys[1:startIndex])
    l2 = ys[startIndex:end]
    threshold = 0.01
    x1 = nothing
    i1 = nothing
    x2 = nothing
    i2 = nothing
    area = nothing
    deviation = nothing
    count = 0
    while isnothing(x1)
        mid = 0.5 * (maxVal + minVal)
        count += 1
        try
            if count > 50
                throw(ErrorException("Failed to converge"))
            end
            li1 = [i for i in 1:length(l1) if l1[i] < mid]
            i1 = startIndex - li1[1]
            li2 = [i for i in 1:length(l2) if l2[i] < mid]
            i2 = startIndex + li2[1]
            area = cs[i2] - cs[i1]
            deviation = abs(area - target)
            if deviation < threshold
                x1 = xs[i1]
                x2 = xs[i2]
            elseif area < target
                maxVal = mid
            elseif area > target
                minVal = mid
            end
        catch e
            @warn "Chain not constrained"
            @debug "Error:\n$e"
            x1 = NaN
            x2 = NaN
        end
    end
    bestfit = [x1, xs[startIndex], x2]
    return bestfit
end

function get_bestfit(chain::Vector{Vector{Float64}}, logdensity::Vector{Float64})
    maxind = argmax(logdensity)
    if length(chain) > length(chain[1])
        bestfit = [get_bestfit([c[i] for c in chain], logdensity, i) for i in 1:length(chain[1])]
    else
        bestfit = [get_bestfit(chain[i], logdensity, i) for i in 1:length(chain)]
    end
    for i in 1:length(bestfit)
        if length(chain) > length(chain[1])
            push!(bestfit[i], chain[maxind][i])
        else
            push!(bestfit[i], chain[i][maxind])
        end
        bestfit[i][1] = bestfit[i][2] - bestfit[i][1]
        bestfit[i][3] = bestfit[i][3] - bestfit[i][2]
    end
    return bestfit 
end


end
