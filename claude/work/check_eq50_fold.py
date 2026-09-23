"""Does the note's range error matter when eq (50) is folded with the reactor spectrum?
Model A fold (Park eq 1, sigma_tot = sigma_KN) with (i) the exact dsigma/dx on [x_-, x_+] and (ii) the note's eq (50) on
the note's range [ (m w + M^2/2)/(2w + m), w ].  1 GW, eps = 1.
Author: Eliana L. De Paoli.
Part of "Sensitivity to dark photons in a reactor antineutrino experiment with Skipper CCDs".
Written with agent assistance; provenance/ records what was derived, recalled and verified.
"""
import numpy as np, json
from scipy.integrate import quad
import check_cross_section as cc
m, alpha = cc.m, cc.alpha
S = cc.S_park

def dN_note(x, M, Emax=40.0):
    # note's range: x <= w  and  x >= (m w + M^2/2)/(2 w + m)  <=>  w >= (x m - M^2/2)/(2x - m)  for x > m/2
    wlo = max(x, (x*m - M*M/2)/(2*x - m)) if x > m/2 else x
    if wlo >= Emax: return 0.0
    f = lambda w: S(w)/cc.sigma_KN(w, m, alpha)*cc.dsig_dx_apuntes(w, x, m, M, alpha)
    return quad(f, wlo, Emax, limit=200)[0]

out = {}
print("threshold implied by the note's range (range non-empty): w >= M/2; with w' >= M imposed: w >= M. Exact: M + M^2/2m")
for M in (0.1, 0.5, 1.0):
    Eth = M + M*M/(2*m)
    # photons the note lets produce A' but that are below the true threshold
    frac_bad = quad(S, M, Eth)[0]/quad(S, M, 40)[0]
    print(f"M={M}: photons with M < w < E_th = {frac_bad*100:.1f}% of those above M")
    rows = []
    for x in (1.0, 2.0, 3.0, 5.0):
        a = cc.dN_dx_compton(x, M); n = dN_note(x, M)
        rows.append((x, a/1e21, n/1e21, n/a if a > 0 else float('nan')))
        print(f"   x={x} MeV: exact fold {a/1e21:.4f}, note-(50) fold {n/1e21:.4f}, ratio {rows[-1][3]:.3f}  [1e21/MeV/s]")
    NA = sum(quad(lambda x: cc.dN_dx_compton(x, M), lo, hi, limit=100)[0] for lo, hi in ((3, 5), (5, 8)))
    NN = sum(quad(lambda x: dN_note(x, M), lo, hi, limit=100)[0] for lo, hi in ((3, 5), (5, 8)))
    print(f"   A' per s in [3,8] MeV: exact {NA:.3e}, note {NN:.3e}, ratio {NN/NA:.3f}")
    out[f"M{M}"] = dict(frac_photons_below_true_threshold=frac_bad, rows=rows, N_window_exact=NA, N_window_note=NN, ratio=NN/NA)
json.dump(out, open("check_eq50_fold_out.json", "w"), indent=1)
