using CSV
using DataFrames
using Distributions
using StatsPlots
using PyCall
using Random
using Transits
using Turing

exo = pyimport("exoplanet")
fits = pyimport("astropy.io.fits")

rootdir(args...) = joinpath(@__DIR__, args...)
datadir(args...) = rootdir("data", args...)

rng = Random.seed!(8462852)

##

t = 0:0.02:80

function exo_to_arr(input)
    X = input[1].eval()
    Y = input[2].eval()
    return @. sqrt(X^2 + Y^2)
end

@model function transit_model(flux, fluxerr)
    logP ~ Normal(log(13), 0.1)
    period = exp(logP)

    t0 ~ Normal()

    us ~ Kipping13()
    r ~ Uniform(0.01, 0.1)
    b ~ Uniform(0, r)

    # define orbit
    orbit = exo.orbits.KeplerianOrbit(;period, t0, b)
    bs = exo_to_arr(orbit.get_relative_position(t))

    law = QuadLimbDark(us)
    y = @. law(bs, r)
    flux ~ MvNormal(y, fluxerr)
    return (;period, t0, us, r, b, flux, fluxerr)
end


function generate_data(yerr)
    err = yerr .* randn(rng, length(t))
    prior = transit_model(missing, err)
    return prior()
end

data = generate_data(5e-4)

phases = @. (t - data.t0 + data.period / 2) % data.period - data.period / 2
scatter(phases, data.flux, yerr=data.fluxerr, xlims=(-0.3, 0.3), c=:black, ms=3)

fits.writeto(datadir("simulated_data_0.fits"), [t data.flux data.fluxerr])
DataFrame(
    snr=5e-4,
    period=data.period,
    t0=data.t0,
    u1=data.us[1],
    u2=data.us[2],
    r=data.r,
    b=data.b,
) |> CSV.write(datadir("params_0.csv"))
