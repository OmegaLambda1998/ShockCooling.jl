# Based of Sapir & Waxman 2017
# DOI: 10.3847/1538-4357/aa64df
export SW17_BSG
mutable struct SW17_BSG <: Model
    name::String
    class::String
    parameter_names::Dict{String,LaTeXString}
    constraints::Dict{String,Tuple{Distribution,Unitful.FreeUnits}}
end

# Default parameter names
function SW17_BSG(name::String, class::String, constraints::Dict{String,Tuple{Distribution,Unitful.FreeUnits}})
    parameter_names = Dict{String,LaTeXString}("R" => L"R_{e}~[R_{\odot}]", "M" => L"M_{e}~[M_{\odot}]", "V" => L"v_{s}~[\frac{km}{s}]", "T" => L"t_{off}~[Days]")
    return SW17_BSG(name, class, parameter_names, constraints)
end

# Default name and parameter names
function SW17_BSG(constraints::Dict{String,Tuple{Distribution,Unitful.FreeUnits}})
    name = "Sapir & Waxman (2017) BSG"
    class = "SW17_BSG"
    return SW17_BSG(name, class, constraints)
end

function bolometric_luminosity(::SW17_BSG, param::Dict{String,Unitful.Quantity{Float64}}, observation::Observation)
    # Load parameters
    r = param["R"] # Envelope radius
    m = param["M"] # Envelope mass
    v = param["V"] # Envelope velocity
    t = (observation.time - param["T"]) # Time since explosion
    t = max(0 * u"d", t)
    if t == 0 * u"d"
        return 0 * u"erg / s"
    end

    # Physical values
    k = 1 # 0.34 cm^2 / g
    mc = 1u"Msun"

    # Rescaled input parameters
    menv = m + mc
    vs = v / (10^8.5 * u"cm / s") |> NoUnits
    M0 = m / 1u"Msun" |> NoUnits
    td = t / 1u"d" |> NoUnits
    R13 = r / (10^13 * u"cm") |> NoUnits

    # n dependent values
    α = 0.73
    a = 4.57
    ϵ2 = -0.175
    fp = 0.08 * m / mc
    Lu = (1.66e42) * u"erg / s"

    # Other values
    ttr = 19.5u"d" * (k * M0 / vs)^0.5

    # Ensure within physical limits
    t_min = 0.2 * (R13 / vs) * max(0.5, (R13^0.4) / (((fp * k * M0)^0.2) * (vs^0.7))) * u"d"
    t_max = 7.4 * ((R13 / k)^0.55) * u"d"
    if (t < t_min) | (t > t_max)
        return 0 * u"erg / s"
    end

    LRW = Lu * vs * vs * R13 * ((vs * td * td / (fp * (menv / 1u"Msun") * k))^ϵ2) / k
    e = -(a * t / ttr)^α

    L = LRW * exp(e)

    return L
end

function radius(model::SW17_BSG, param::Dict{String,Unitful.Quantity{Float64}}, observation::Observation)
    L = bolometric_luminosity(model, param, observation)
    T = temperature(model, param, observation)
    σ = 5.6704e-8 * u"W / m^2 / K^4" # Stefan Boltzmann constant

    R = (L / (4π * σ * T^4))^0.5

    return R
end

function temperature(::SW17_BSG, param::Dict{String,Unitful.Quantity{Float64}}, observation::Observation)
    r = param["R"] # Envelope radius
    m = param["M"] # Envelope mass
    v = param["V"] # Envelope velocity
    t = observation.time - param["T"] # Time since explosion
    t = max(0 * u"d", t)

    # Physical values
    k = 1
    mc = 1u"Msun"
    menv = m + mc

    # Rescaled input parameters
    vs = v / (10^8.5 * u"cm / s") |> NoUnits
    td = t / 1u"d" |> NoUnits

    M0 = m / 1u"Msun" |> NoUnits

    R13 = r / (10^13 * u"cm") |> NoUnits


    # n dependent values
    Tu = 1.96e4 * u"K"
    fp = 0.08 * m / mc
    ϵ1 = 0.016

    T = Tu * ((vs * vs * td * td / (fp * (menv / 1u"Msun") * k))^ϵ1) * ((R13 / (td * td * k))^0.25)
    return T
end

function model_flux(model::SW17_BSG, param::Dict{String,Unitful.Quantity{Float64}}, observation::Observation)
    T = temperature(model, param, observation)
    return synthetic_flux(observation.filter, T)
end

function run_model(model::SW17_BSG, param::Dict{String,Unitful.Quantity{Float64}}, supernova::Supernova)
    m_flux = [model_flux(model, param, obs) for obs in supernova.lightcurve.observations] .|> u"erg / s / cm^2 / Hz"
    dist = 10u"pc"
    R = [radius(model, param, obs) for obs in supernova.lightcurve.observations]
    abs_mag = @. -48.6 - 2.5 * (log10(ustrip(m_flux)) + log10(uconvert(NoUnits, R / dist)^2))
    replace!(abs_mag, NaN => -10)
    return abs_mag * u"AB_mag"
end
