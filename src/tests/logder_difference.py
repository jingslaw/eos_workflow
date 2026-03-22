from abipy.ppcodes.oncv_parser import OncvParser
import numpy as np

pspout = 'Eu-4f-o1.out'

p = OncvParser(pspout)
p.scan(verbose=1)
logders = p.atan_logders
lmax = p.lmax

for l in logders.ae:
    # diff with counting the weight on fermi dirac distribution
    f1, f2 = logders.ae[l], logders.ps[l]

    abs_diff = np.abs(f1.values - f2.values)
