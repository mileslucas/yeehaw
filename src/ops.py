import theano.tensor as tt
from theano.graph.op import Op
from theano.compile.ops import as_op

import batman
import pytransit
import numpy as np

with np.load("../data/generated_data.npz") as data:
    t = data["t"]

tm = pytransit.QuadraticModel()
tm.set_data(t)

def PyTransitlightcurve_np(t0, period, a, r, us):
       return tm.evaluate(r, us, t0, period, a, np.pi/2)

@as_op(itypes = [tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar, tt.dvector],
       otypes = [tt.dvector])
def PyTransitlightcurve(*inputs):
       np_ins = [np.atleast_1d(inp) for inp in inputs]
       return PyTransitlightcurve_np(*np_ins)


pars = batman.TransitParams()

@as_op(itypes = [tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar, tt.dvector],
       otypes = [tt.dvector])
def BATMANlightcurve(t0, period, a, r, us):
    # set up batman parameter structure
    pars.t0 = t0
    pars.per = period
    pars.a = a
    pars.rp = r
    pars.u = us
    # default orientation with circular orbit
    pars.inc = 90
    pars.ecc = 0
    pars.w = 90
    pars.limb_dark = "quadratic"
    m = batman.TransitModel(pars, t)
    return m.light_curve(pars)