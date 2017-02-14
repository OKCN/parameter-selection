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

def okcn_get_d(q, n, m, g):
    d = max(filter(lambda x: m*(x*2 + 1) < q*(1. - 1./g), range(q)))
    return d

def akcn_get_d(q, n, m, g):
    d = max(filter(lambda x: m*(x*2 + 1) < q*(1. - m*1./g), range(q)))
    return d

def main():
    print "RLWE with OKCN, AKCN"

    bandwith = lambda q, n, g : (ceil(log(q, 2)) * n * 2 + ceil(log(g, 2)) * n + 256)*1. / 8

    q = 12289; n = 1024; m = 2; 
    reclen = n; 
    noise = sym_binomial(32)
    noise_sqr = pdf_product(noise, noise, q)  
    v = nfoldconvolution(2*n, noise_sqr, q)
    v = convolution(v, noise, q)  # v = 2*n * (noise * noise) + noise

    g = 2**4
    d = okcn_get_d(q, n, m, g)
    failure_pr = reclen * sum(map(lambda x: v[x], filter(lambda x: min(x, q - x) > d, v.keys())))
    print "OKCN q = {}, n = {}, m = {}, g = {}, d = {}, bandwith = {}".format(q, n, m, g, d, bandwith(q, n, g))
    print "exact pr of failure = 2^{:.2f}".format(log(failure_pr, 2))
    print 

    g = 2**6
    d = okcn_get_d(q, n, m, g)
    failure_pr = reclen * sum(map(lambda x: v[x], filter(lambda x: min(x, q - x) > d, v.keys())))
    print "OKCN q = {}, n = {}, m = {}, g = {}, d = {}, bandwith = {}".format(q, n, m, g, d, bandwith(q, n, g))
    print "exact pr of failure = 2^{:.2f}".format(log(failure_pr, 2))
    print 

    g = 2**4
    d = akcn_get_d(q, n, m, g)
    failure_pr = reclen * sum(map(lambda x: v[x], filter(lambda x: min(x, q - x) > d, v.keys())))
    print "AKCN q = {}, n = {}, m = {}, g = {}, d = {}, bandwith = {}".format(q, n, m, g, d, bandwith(q, n, g))
    print "exact pr of failure = 2^{:.2f}".format(log(failure_pr, 2))
    print 

    g = 2**6
    d = akcn_get_d(q, n, m, g)
    failure_pr = reclen * sum(map(lambda x: v[x], filter(lambda x: min(x, q - x) > d, v.keys())))
    print "AKCN q = {}, n = {}, m = {}, g = {}, d = {}, bandwith = {}".format(q, n, m, g, d, bandwith(q, n, g))
    print "exact pr of failure = 2^{:.2f}".format(log(failure_pr, 2))
    print 

if __name__ == "__main__":
    main()
