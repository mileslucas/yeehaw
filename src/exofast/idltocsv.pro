pro idltocsv, filename=filename

;; load in the MCMC output structure, called `mcmcss`
restore, filename

;; get chi2 and find good indices
chi2 = *(mcmcss.chi2)
minchi2 = min(chi2, ndx)
; reform to (nsteps, nchains) array
chi2 = reform(chi2, mcmcss.nsteps/mcmcss.nchains, mcmcss.nchains)
burnndx = getburnndx(chi2, goodchains=goodchains)

;; pick output parameters
period = mcmcss.planet[0].period
t0 = mcmcss.planet[0].tc
r = mcmcss.planet[0].p
ar = mcmcss.planet[0].ar
u1 = mcmcss.band[0].u1
u2 = mcmcss.band[0].u2
var = mcmcss.transit[0].variance

;; reform to (nsteps, nchains) arrays
period_2d = reform(period.value, mcmcss.nsteps/mcmcss.nchains, mcmcss.nchains)
t0_2d = reform(t0.value, mcmcss.nsteps/mcmcss.nchains, mcmcss.nchains)
r_2d = reform(r.value, mcmcss.nsteps/mcmcss.nchains, mcmcss.nchains)
ar_2d = reform(ar.value, mcmcss.nsteps/mcmcss.nchains, mcmcss.nchains)
u1_2d = reform(u1.value, mcmcss.nsteps/mcmcss.nchains, mcmcss.nchains)
u2_2d = reform(u2.value, mcmcss.nsteps/mcmcss.nchains, mcmcss.nchains)
sigma_2d = reform(sqrt(var.value), mcmcss.nsteps/mcmcss.nchains, mcmcss.nchains)

;; remove burn-in and any bad chains
period_trim = period_2d[burnndx:-1, goodchains]
t0_trim = t0_2d[burnndx:-1, goodchains]
r_trim = r_2d[burnndx:-1, goodchains]
ar_trim = ar_2d[burnndx:-1, goodchains]
u1_trim = u1_2d[burnndx:-1, goodchains]
u2_trim = u2_2d[burnndx:-1, goodchains]
sigma_trim = sigma_2d[burnndx:-1, goodchains]

;; reform back to 1d arrays
N = size(period_trim, /N_ELEMENTS)
period_1d = reform(period_trim, N)
t0_1d = reform(t0_trim, N)
r_1d = reform(r_trim, N)
ar_1d = reform(ar_trim, N)
u1_1d = reform(u1_trim, N)
u2_1d = reform(u2_trim, N)
sigma_1d = reform(sigma_trim, N)

;; save to CSV file
outpath = repstr(filename, ".idl", ".csv")
hdr = ["period", "t0", "r", "a", "u1", "u2", "yerr"]
write_csv, outpath, period_1d, t0_1d, r_1d, ar_1d, u1_1d, u2_1d, sigma_1d, header=hdr

end
