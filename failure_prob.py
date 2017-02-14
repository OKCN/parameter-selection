"""Computing failure probabilities (heuristic and exact) of the key agreement protocol.

Based on the paper:
    Joppe Bos, Craig Costello, Leo Ducas, Ilya Mironov, Michael Naehrig, Valeria
    Nikolaenko, Ananth Raghunathan, Douglas Stebila.  Frodo: Take off the ring!
    Practical, quantum-secure key exchange from LWE.  In ACM Conference on Computer
    and Communications Security (CCS) 2016, ACM, October, 2016.
    DOI: http://dx.doi.org/10.1145/2976749.2978425
    Eprint http://eprint.iacr.org/2016/659

Copyright (c) 2016 Joppe Bos, Leo Ducas, Ilya Mironov, Valeria Nikolaenko,
                   Ananth Raghunathan, Douglas Stebila

Released under the MIT License; see LICENSE.txt for details.
"""

from math import sqrt, exp, log, ceil
from discrete_distr import pdf_product, std_modulo, nfoldconvolution, convolution, sym_binomial, sym_uniform, dgauss

def noise_failure_prob(noise, q, n, w, reclen, g, cut):
    m = 2**w
    d = max(filter(lambda x: m*x*2 < q*(1. - 1./g), range(q)))
    print "d = ", d

    eps = sym_uniform(2**cut)
    noise_sqr = pdf_product(noise, noise, q)                       # noise_sqr = noise * noise
    noise_eps = pdf_product(noise, convolution(noise, eps, q), q)  # noise_eps = noise * (noise + eps)
    v1 = nfoldconvolution(n, noise_sqr, q) # v1 = n * (noise * noise)
    v2 = nfoldconvolution(n, noise_eps, q) # v2 = n * (noise * (noise + eps))
    vv = convolution(convolution(v1, v2, q), noise, q)  # vv = v1 + v2 + noise

    failure_pr = reclen * sum(map(lambda x: vv[x], filter(lambda x: min(x, q - x) > d, vv.keys())))

    return failure_pr
