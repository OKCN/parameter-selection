""" Classical, quantum, and plausible (conservative) quantum cost estimates.

    Modified from the scripts of Frodo

    Released under the MIT License; see LICENSE.txt for details.
"""

from math import pi, exp, log, sqrt, ceil, floor, erfc

log_infinity = 9999


def delta_BKZ(b):
    """Computes the root hermite factor delta of BKZ-b.

    Args:
      b (integer): blocksize

    Returns:
      delta (real)
    """
    return ((pi * b)**(1. / b) * b / (2 * pi * exp(1)))**(1. / (2. * b - 2.))


def svp_plausible(b):
    """Computes the best plausible quantum cost of SVP in dimension b.

    Args:
      b (integer): blocksize

    Returns:
      log_2 of the cost.
    """
    return log(sqrt(4. / 3), 2) * b + log(b, 2)


def svp_quantum(b):
    """Computes the best known quantum cost of SVP in dimension b.

    Args:
      b (integer): blocksize

    Returns:
      log_2 of the cost.
    """
    return log(sqrt(13. / 9), 2) * b + log(b, 2)


def svp_classical(b):
    """Computes the best known classical cost of SVP in dimension b.

    Args:
      b (integer): blocksize

    Returns:
      log_2 of the cost.
    """
    return log(sqrt(3. / 2), 2) * b + log(b, 2)


def nvec_sieve(b):
    """Computes the number of vectors output by the sieving step for SVP in dimension b.

    Args:
      b (integer): blocksize

    Returns:
      log_2 of the number of vectors (proxy for the cost).
    """
    return log(sqrt(4. / 3), 2) * b


def primal_cost(q, n, m, s, b, cost_svp=svp_classical, verbose=False):
    """Computes the cost of a primal attack using m samples and blocksize b (infinity if the attack fails).

    The attack assumes 'small secret' (i.e., the distribution of the secret is the same as errors).
    Args:
      q: LWE modulus
      n: LWE dimension
      m: number of samples
      s: standard deviation of the error distribution
      b: blocksize
      cost_svp: cost function for the SVP oracle
      verbose: spill details on stdout

    Returns:
      log_2 of the number of vectors
    """
    d = n + m
    delta = delta_BKZ(b)
    if verbose:
        print "Primal attacks uses block-size", b, "and ", m, "samples"

    if s * sqrt(b) < delta**(2. * b - d - 1) * q**(1. * m / d):
        return cost_svp(b)  # + log(n-b)/log(2)

    else:
        return log_infinity

def dual_cost(q, n, m, s, b, cost_svp=svp_classical, verbose=False):
    """Computes the cost a dual attack using m samples and blocksize b (infinity of the attacks fails).

    The attack assume 'small secret' (i.e., the distribution of the secret is the same as errors).
    Args:
      q: LWE modulus
      n: LWE dimension
      m: number of samples
      s: standard deviation of the error distribution
      b: blocksize
      cost_svp: cost function for the SVP oracle
      verbose: spill details on stdout

    Returns:
      log_2 of the cost.
    """

    d = n + m
    delta = delta_BKZ(b)
    l = delta**d * q**(1. * n / d)

    tau = l * s / q

    log2_eps = 2 + (- 2 * pi * pi * tau * tau) / log(2)

    cost = max(0, - 2 * log2_eps - nvec_sieve(b)) + cost_svp(b)
    if verbose:
        print "Dual attacks uses block-size", b, "and ", m, "samples"
        print "log2(epsilon) = ", log2_eps, "log2 nvector per run", nvec_sieve(b)
    return cost

def optimize_attack(
        q,
        n,
        max_m,
        s,
        cost_attack=primal_cost,
        cost_svp=svp_classical,
        verbose=False):
    """ Finds optimal attack parameters against an LWE instance.

      q: LWE modulus
      n: LWE dimension
      max_m: maximum number of samples
      s: standard deviation of the error distribution
      cost_attack: primal vs dual attack
      cost_svp: cost function for the SVP oracle
      verbose: spill details on stdout

    Returns:
      Best parameters (m and b) and log_2 of their cost.
    """
    best_cost = log_infinity
    best_b = 0
    best_m = 0
    for b in range(50, n + max_m + 1):
        # Try all block-size, until just one call to SVP cost more than the
        # best attack found so far
        if cost_svp(b) > best_cost:
            break
        # Try all possible number of number of sample
        for m in range(max(1, b - n), max_m):
            cost = cost_attack(q, n, m, s, b, cost_svp)
            if cost < best_cost:
                best_cost = cost
                best_m = m
                best_b = b
    if verbose:
        # Re-call the cost_attack on the best params spilling extra info
        cost_attack(q, n, best_m, s, best_b, cost_svp=svp_classical, verbose=True)
    return (best_m, best_b, best_cost)

def weighted_primal_cost(q, n, m, s_e, s_s, b, cost_svp = svp_classical):
    w = 1. * s_e / s_s
    s = s_s
    q_hat = int(round(1.*q / w))
    #print "weighted_primal q, s = ", q_hat, s
    return primal_cost(q_hat, n, m, s, b, cost_svp)

def weighted_dual_cost(q, n, m, s_e, s_s, b, cost_svp = svp_classical):
    c = 1. * s_e / s_s
    s = s_e

    d = n + m
    delta = delta_BKZ(b)
    l = delta**d * (1.*q/c)**(1. * n / d)

    tau = l * s_e / q

    log2_eps = 2 + (- 2 * pi * pi * tau * tau) / log(2)

    cost = max(0, - 2 * log2_eps - nvec_sieve(b)) + cost_svp(b)
    return cost

def weighted_optimize_attack(
        q,
        n,
        max_m,
        s_e,
        s_s,
        cost_attack,
        cost_svp):
    best_cost = log_infinity
    best_b = 0
    best_m = 0
    for b in range(50, n + max_m + 1):
        # Try all block-size, until just one call to SVP cost more than the
        # best attack found so far
        if cost_svp(b) > best_cost:
            break
        # Try all possible number of number of sample
        for m in range(max(1, b - n), max_m):
            cost = cost_attack(q, n, m, s_e, s_s, b, cost_svp)
            if cost < best_cost:
                best_cost = cost
                best_m = m
                best_b = b
    return (best_m, best_b, best_cost)

bandwidth = lambda p, q, n, l, g: (ceil(log(p, 2))*n*l*2 + ceil(log(g, 2)*l*l)) / 8000.


def main():
    #parameters = [ {'q': 2**15, 'p': 2**12, 'n': 672, 'samples': 672*8*2+64, 's_s': sqrt(2.0)},
    parameters = [ {'q': 2**15, 'p': 2**12, 'n': 680, 'samples': 680*8*2+64, 's_s': sqrt(1.7)},
                   {'q': 2**15, 'p': 2**12, 'n': 832, 'samples': 832*8*2+64, 's_s': sqrt(1.4)} ]

    for param in parameters:
        print "LWR (Rounded Gaussian)"
        q, p, n, samples, s_s = param['q'], param['p'], param['n'], param['samples'], param['s_s']

        print 'q = {}, p = {}, n = {}, samples = {}, s_s = {}'.format(q, p, n, samples, s_s)
        s_e = 1.*q / (sqrt(12.) * p)

        print "bw. = ", bandwidth(p, q, n, 8, 2**8)
        print "q_hat, s_e, s_s = ", 1.*q/(s_e / s_s), s_e, s_s

        m, b, c = weighted_optimize_attack(q, n, samples, s_e, s_s, cost_attack=weighted_primal_cost, cost_svp=svp_classical)
        print "Primal m = {}, b = {}, C = {:.2f}".format(m, b, c)
        m, b, c = weighted_optimize_attack(q, n, samples, s_e, s_s, cost_attack=weighted_dual_cost, cost_svp=svp_classical)
        print "Dual m = {}, b = {}, C = {:.2f}".format(m, b, c)

        m, b, c = weighted_optimize_attack(q, n, samples, s_e, s_s, cost_attack=weighted_primal_cost, cost_svp=svp_quantum)
        print "Primal m = {}, b = {}, Q = {:.2f}".format(m, b, c)
        m, b, c = weighted_optimize_attack(q, n, samples, s_e, s_s, cost_attack=weighted_dual_cost, cost_svp=svp_quantum)
        print "Dual m = {}, b = {}, Q = {:.2f}".format(m, b, c)

        m, b, c = weighted_optimize_attack(q, n, samples, s_e, s_s, cost_attack=weighted_primal_cost, cost_svp=svp_plausible)
        print "Primal m = {}, b = {}, P = {:.2f}".format(m, b, c)
        m, b, c = weighted_optimize_attack(q, n, samples, s_e, s_s, cost_attack=weighted_dual_cost, cost_svp=svp_plausible)
        print "Dual m = {}, b = {}, P = {:.2f}".format(m, b, c)

        print

if __name__ == "__main__":
    main()
