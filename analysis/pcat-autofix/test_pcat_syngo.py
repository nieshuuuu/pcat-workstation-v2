"""Brute-force check that the min-cut really returns the global optimum.

The Ishikawa construction has two failure modes that are completely silent:

  - scipy's maximum_flow accepts an int64 matrix and truncates it to int32 internally,
    so a sentinel like 1 << 40 becomes a capacity of ZERO and every constraint in the
    graph quietly disappears.  The solver returns a flow of 0 and no error.
  - the pairwise arcs do not realise g(a - b).  They realise g(a - b) + f(a) + f(b) for
    a separable f that has to be subtracted back out of the unary costs.  Forget it and
    the energy is wrong by a term that looks plausible.

Both were present in the first version of pcat_syngo and neither showed up as anything
other than slightly-odd contours.  Hence this file.

    python test_pcat_syngo.py
"""
import sys, itertools
import numpy as np
sys.path.insert(0, __import__("os").path.dirname(__file__))
import numpy as np
sys.path.insert(0, r"Z:\Shu Nie\shu_pcat\tools")
import pcat_syngo as sg

def energy(lab, cost, delta, gt, gz):
    Z, A, R = cost.shape
    def hub(d, gam):
        d = abs(float(d))
        return gam * (d * d / (2 * delta) if d <= delta else d - delta / 2)
    e = sum(cost[z, a, lab[z, a]] for z in range(Z) for a in range(A))
    if A >= 3:
        for z in range(Z):
            for a in range(A):
                e += hub(lab[z, a] - lab[z, (a + 1) % A], gt)
    for z in range(Z - 1):
        for a in range(A):
            e += hub(lab[z, a] - lab[z + 1, a], gz)
    return e

rng = np.random.default_rng(11)
bad = worst = 0
n = 0
for Z, A, R in [(1, 3, 6), (1, 4, 5), (1, 5, 4), (2, 3, 5), (2, 4, 4), (1, 4, 6)]:
    for delta in (1.0, 2.0, 3.0):
        for gt, gz in ((0.5, 0.1), (2.0, 1.0), (0.1, 0.9)):
            cost = rng.normal(size=(Z, A, R)) * rng.uniform(0.3, 3.0)
            lab = sg.optimal_surface(cost, gt, gz, delta)
            e_cut = energy(lab, cost, delta, gt, gz)
            eb = min(energy(np.array(c).reshape(Z, A), cost, delta, gt, gz)
                     for c in itertools.product(range(R), repeat=Z * A))
            n += 1
            d = e_cut - eb
            worst = max(worst, d)
            bad += d > 1e-6
print("%d configurations, %d mismatched, worst excess energy %.3e" % (n, bad, worst))
