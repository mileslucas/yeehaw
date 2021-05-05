# using ArviZ
using CSV
using DataFrames
using Distributions
using NPZ
using StatsPlots
using Random
using Transits
using Turing
using NestedSamplers


rootdir(args...) = joinpath(@__DIR__, args...)
datadir(args...) = rootdir("data", args...)

data = npzread(datadir("generated_data_100ppm.npz"))
t = data["t"]
flux = data["flux"]
data = npzread(datadir("generated_data.npz"))
true_mod = data["flux"]

ground_truth = Dict(
    "period" => 3.5,
    "t0" => 1.3,
    "us" => [0.5, 0.2],
    "r" => 0.03,
    "a" => 10.0,
    "yerr" => 1e-4
)

phases = @. (t - ground_truth["t0"] + ground_truth["period"] / 2) % ground_truth["period"] - ground_truth["period"] / 2
# phases, flux_phase, mod_phase = Transits.phaseup(phases, flux, true_mod; period=ground_truth["period"], t0=ground_truth["t0"])
maxphase = 0.1
inds = @. abs(phases) < maxphase


# show data
scatter(phases[inds], flux[inds], yerr=ground_truth["yerr"], c=:black, label="data")
plot!(phases[inds], true_mod[inds], label="truth")

orbit = KeplerianOrbit(ground_truth["a"], ground_truth["period"], 0, ground_truth["t0"], 0)
law = PolynomialLimbDark(ground_truth["us"])
f = @. law(orbit, t, ground_truth["r"])
plot(phases[inds], f[inds], label="transits")
##
# using Zygote
# Turing.setadbackend(:zygote)
rng = Random.seed!(8462852)

@model function transit_model(t, flux)
    logP ~ Normal(log(ground_truth["period"]), 0.1)
    period = exp(logP)

    t0 ~ Normal(ground_truth["t0"], 10)

    us ~ Kipping13()
    r ~ Uniform(0.01, 0.1)
    loga ~ Normal(log(ground_truth["a"]), 0.1)
    a = exp(loga)

    # define orbit
    orbit = KeplerianOrbit(a, period, 0, t0, 0)
    law = PolynomialLimbDark(us)

    yerr ~ Truncated(Cauchy(0, 0.1), 0, Inf)
    for i in eachindex(flux)
        flux[i] ~ Normal(law(orbit, t[i], r), yerr)
    end
end

model = transit_model(t, flux)

chain = sample(rng, model, NUTS(5000, 0.9), MCMCThreads(), 5000, 4)