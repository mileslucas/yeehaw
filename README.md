# yeehaw

A round-up of probabilistic transit-fitting software.

## Software tested
* [exoplanet](https://github.com/exoplanet-dev/exoplanet)
* [BATMAN](https://github.com/lkreidberg/batman)
* [PyTransit](https://github.com/hpparvi/PyTransit)
* [Juliet](https://github.com/nespinoza/juliet)
* [EXOFASTv2](https://github.com/jdeast/EXOFASTv2)


## Setup

1. Install Python requirements (preferrably inside virtual environment)

Dependencies are resolved using [poetry](https://github.com/python-poetry/poetry)
```
$ poetry install
```
2. (Optional) install Julia dependencies
```
$ julia --project=@. -e 'using Pkg; Pkg.instantiate()'
```

Most of the code is contained in jupyter notebooks, so you will have to set up an interactive environment to interface with them.

## Getting Started

Data was simulated using exoplanet (ALFM 20) light curves and subsequently fit by various packages. The data can be found in [data](data/), which is produced by the [exoplanet notebook](https://github.com/mileslucas/yeehaw/blob/master/src/exoplanet.ipynb). For packages which do not provide statistical modeling, like BATMAN, and PyTransit, a model is built using PyMC3. This requires specifying theano ops for the light curves, which are in [src/ops.py](src/ops.py). exoplanet, Juliet, and EXOFASTv2 all have statistical frameworks built into or on top of their respective packages. The [comparison notebook](https://github.com/mileslucas/yeehaw/blob/master/src/comparison.ipynb) compares the posterior samples for all of the tested packages and models.

## License

The code produced here is MIT licensed, any package code is subject to its respective license.