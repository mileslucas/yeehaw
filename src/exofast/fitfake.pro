pro fitfake, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, $
             nthread=nthread, outpath=outpath

path = filepath('',root_dir='/Users/miles/dev/yeehaw',subdir=['src','exofast'])

if n_elements(outpath) eq 0 then $
   outpath = '/Users/miles/dev/yeehaw/data/exofast' + path_sep() + 'generated_100ppm' + path_sep()

exofastv2, nplanets=1, tranpath=path+'n20210504.TESS.fake.nonoise.dat', $
           priorfile=path+'fake.priors', $
           prefix=outpath +'generated_100ppm.', maxsteps=maxsteps, $
           nthin=nthin, circular=[1], fitrv=[0],fittran=[1], $
           debug=debug, verbose=verbose, nthread=nthread, /NOMIST

end
