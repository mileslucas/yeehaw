# using ArviZ
using CSV
using DataFrames
using Distributions
using NPZ
using StatsPlots
using Random
using Transits
using NestedSamplers
using StatsBase


rootdir(args...) = joinpath(@__DIR__, args...)
datadir(args...) = rootdir("data", args...)

data = npzread(datadir("generated_data_100ppm.npz"))
const t = data["t"]
const flux = data["flux"]
data = npzread(datadir("generated_data.npz"))
const true_mod = data["flux"]

ground_truth = Dict(
    "period" => 3.5,
    "t0" => 1.3,
    "us" => [0.5, 0.2],
    "r" => 0.03,
    "a" => 10.0,
    "yerr" => 1e-4
)

phases = @. (t - ground_truth["t0"] + ground_truth["period"] / 2) % ground_truth["period"] - ground_truth["period"] / 2
maxphase = 0.1
inds = @. abs(phases) < maxphase


# show data
scatter(phases[inds], flux[inds], yerr=ground_truth["yerr"], c=:black, label="data")
plot!(phases[inds], true_mod[inds], label="truth")

orbit = KeplerianOrbit(ground_truth["a"], ground_truth["period"], 0, ground_truth["t0"], 0)
law = PolynomialLimbDark(ground_truth["us"])
f = @. law(orbit, t, ground_truth["r"])
plot(phases[inds], f[inds], label="transits")

# Define a function to optimize.
function loglike(X)
    logP, t0, r, loga, yerr = @view X[1:5]
    us = @view X[end-1:end]
    
    orbit = KeplerianOrbit(exp(loga), exp(logP), 0, t0, 0)
    ld = QuadLimbDark(us)
    ll = zero(eltype(X))
    for idx in eachindex(flux)
        dist = Normal(ld(orbit, t[idx], r), yerr)
        ll += logpdf(dist, flux[idx])
    end
    return ll
end

function lightcurve(X)
    logP, t0, r, loga, yerr = @view X[1:5]
    us = @view X[end-1:end]
    
    orbit = KeplerianOrbit(exp(loga), exp(logP), 0, t0, 0)
    ld = QuadLimbDark(us)
    return @. ld(orbit, t, r)
end

X0 = [log(3.5), 1.3, 0.03, log(10), 1e-4, 0.5, 0.2]

## Distributions
priors = Dict(
    :logP => Normal(log(ground_truth["period"]), 0.1),
    :t0 => Normal(ground_truth["t0"], 10),
    :r => Uniform(0.01, 0.1),
    :loga => Normal(log(ground_truth["a"]), 0.1),
    :yerr => Truncated(Cauchy(0, 0.1), 0, Inf),
)

function kiptform!(X)
    sqrtq1 = sqrt(first(X))
    twoq2 = 2 * last(X)
    X[begin] = sqrtq1 * twoq2
    X[end] = sqrtq1 * (1 - twoq2)
    return X
end

kiptform(X) = kiptform!(copy(X))

function prior_transform(X)
    out = copy(X)
    for (i, k) in enumerate((:logP, :t0, :r, :loga, :yerr))
        out[i] = quantile(priors[k], X[i])
    end
    us = @view out[end-1:end]
    kiptform!(us)
    return out
end

## sampling
model = NestedModel(loglike, prior_transform)
splr = Nested(7, 1000; proposal=Proposals.RWalk(ratio=0.9))

chain, state = sample(model, splr, param_names=["logP", "t0", "r", "loga", "yerr", "us[1]", "us[2]"])
