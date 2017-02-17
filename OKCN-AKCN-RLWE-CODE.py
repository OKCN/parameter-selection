"""

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

def getVariance(v, q):
    return sum(v[i] * min(i, q - i)**2 for i in v.keys())

def main():
    print "RLWE with OKCN, AKCN"

    bandwith = lambda q, n, g : (ceil(log(q, 2)) * n * 2 + ceil(log(g, 2)) * n + 256)*1. / 8

    q, n, m = 12289, 1024, 2
    reclen = n; 
    noise = sym_binomial(32)
    noise_sqr = pdf_product(noise, noise, q)  
    v = nfoldconvolution(2*n, noise_sqr, q)
    v = convolution(v, noise, q)  # v = 2*n * (noise * noise) + noise

    print getVariance(v, q), 2.*n*getVariance(noise, q)**2 + getVariance(noise, q)

    g = 2**2
    d = okcn_get_d(q, n, m, g)
    failure_pr = sum(map(lambda x: v[x], filter(lambda x: min(x, q - x) > d, v.keys())))
    print "OKCN q = {}, n = {}, m = {}, g = {}, d = {}, bandwith = {}".format(q, n, m, g, d, bandwith(q, n, g))
    print "per bit pr of failure = 2^{:.2f}".format(log(failure_pr, 2))
    print 

    g = 2**3
    d = okcn_get_d(q, n, m, g)
    failure_pr = sum(map(lambda x: v[x], filter(lambda x: min(x, q - x) > d, v.keys())))
    print "OKCN q = {}, n = {}, m = {}, g = {}, d = {}, bandwith = {}".format(q, n, m, g, d, bandwith(q, n, g))
    print "per bit pr of failure = 2^{:.2f}".format(log(failure_pr, 2))
    print 

    g = 2**4
    d = okcn_get_d(q, n, m, g)
    failure_pr = sum(map(lambda x: v[x], filter(lambda x: min(x, q - x) > d, v.keys())))
    print "OKCN q = {}, n = {}, m = {}, g = {}, d = {}, bandwith = {}".format(q, n, m, g, d, bandwith(q, n, g))
    print "per bit pr of failure = 2^{:.2f}".format(log(failure_pr, 2))
    print 

    g = 2**4
    d = akcn_get_d(q, n, m, g)
    failure_pr = sum(map(lambda x: v[x], filter(lambda x: min(x, q - x) > d, v.keys())))
    print "AKCN q = {}, n = {}, m = {}, g = {}, d = {}, bandwith = {}".format(q, n, m, g, d, bandwith(q, n, g))
    print "per bit pr of failure = 2^{:.2f}".format(log(failure_pr, 2))
    print 

    g = 2**6
    d = akcn_get_d(q, n, m, g)
    failure_pr = sum(map(lambda x: v[x], filter(lambda x: min(x, q - x) > d, v.keys())))
    print "AKCN q = {}, n = {}, m = {}, g = {}, d = {}, bandwith = {}".format(q, n, m, g, d, bandwith(q, n, g))
    print "per bit pr of failure = 2^{:.2f}".format(log(failure_pr, 2))
    print 

if __name__ == "__main__":
    main()
