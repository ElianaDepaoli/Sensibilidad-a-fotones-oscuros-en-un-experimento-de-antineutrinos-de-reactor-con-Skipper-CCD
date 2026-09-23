#!/usr/bin/env python3
"""Where x_- and x_+ come from, and that the two routes to them agree.

Route A (the one used in report_cross_section eq (5)): two-body momentum in the CM,
boosted to the lab. Route B (the one the Spanish note takes, but done with the exact
kinematics): set cos(theta) = +-1 in the exact angle-energy relation and solve for x.

Everything is asserted, symbolically with sympy where it is an identity and numerically
where it is an inequality. Run:  work/.venv/bin/python work/check_xpm.py [out.json]

Author: Eliana L. De Paoli.
Part of "Sensitivity to dark photons in a reactor antineutrino experiment with Skipper CCDs".
Written with agent assistance; provenance/ records what was derived, recalled and verified.
"""
import json
import sys

import numpy as np
import sympy as sp

OUT = {}

m_, M_, E_, x_ = sp.symbols('m M E x', positive=True)
s_ = m_**2 + 2*m_*E_                      # (p+k)^2 with the electron at rest
lam_ = (s_ - (m_ + M_)**2)*(s_ - (m_ - M_)**2)   # Kallen lambda(s, m^2, M^2)

# The published answer, eq (5) of reports/report_cross_section_*.tex
xp_eq5 = ((E_ + m_)*(s_ + M_**2 - m_**2) + E_*sp.sqrt(lam_))/(2*s_)
xm_eq5 = ((E_ + m_)*(s_ + M_**2 - m_**2) - E_*sp.sqrt(lam_))/(2*s_)


def section(t):
    print(f"\n{t}\n" + "-"*len(t))


# ---------------------------------------------------------------------------
section("1. Route A: CM two-body momentum, boosted to the lab")
# In the CM the A' energy and momentum are fixed by s alone:
Eq_cm = (s_ + M_**2 - m_**2)/(2*sp.sqrt(s_))
q_cm = sp.sqrt(lam_)/(2*sp.sqrt(s_))
# The CM moves in the lab with beta = E/(E+m)  (total momentum E over total energy E+m):
gamma = (E_ + m_)/sp.sqrt(s_)
betagamma = E_/sp.sqrt(s_)
# x = gamma*Eq_cm + betagamma*q_cm*cos(theta_cm), extremes at cos(theta_cm) = +-1
xp_boost = gamma*Eq_cm + betagamma*q_cm
xm_boost = gamma*Eq_cm - betagamma*q_cm
assert sp.simplify(xp_boost - xp_eq5) == 0
assert sp.simplify(xm_boost - xm_eq5) == 0
# the boost factors themselves are not assumed: check gamma^2 - (betagamma)^2 = 1
assert sp.simplify(gamma**2 - betagamma**2 - 1) == 0
print("   PASS: gamma*Eq_cm +- betagamma*q_cm == eq (5), and gamma^2-(bg)^2 = 1")

# ---------------------------------------------------------------------------
section("2. Route B: cos(theta) = +-1 in the exact lab relation")
Q_ = sp.sqrt(x_**2 - M_**2)
# (p + k - q)^2 = m^2  ->  M^2 + 2m(E-x) - 2E(x - Q cos) = 0
onshell = M_**2 + 2*m_*(E_ - x_) - 2*E_*(x_ - Q_*sp.Symbol('c'))
cos_exact = sp.solve(onshell, sp.Symbol('c'))[0]
assert sp.simplify(cos_exact - ((x_ - m_)*E_ + m_*x_ - M_**2/2)/(E_*Q_)) == 0
print("   PASS: solving the on-shell condition for cos gives eq (4) of the report")

# cos = +-1 : M^2 + 2mE - 2(m+E)x = -+ 2E sqrt(x^2-M^2). Square it.
A = M_**2 + 2*m_*E_ - 2*(m_ + E_)*x_
squared = sp.expand(A**2 - 4*E_**2*(x_**2 - M_**2))
poly = sp.Poly(squared, x_)
c2, c1, c0 = poly.all_coeffs()
assert sp.simplify(c2 - 4*s_) == 0                       # coefficient of x^2 is 4s
assert sp.simplify(c1 + 4*(m_ + E_)*(s_ + M_**2 - m_**2)) == 0
# discriminant is E^2 * lambda, i.e. the square root of eq (5) is not an accident
disc = sp.simplify(c1**2 - 4*c2*c0)
assert sp.simplify(disc - 16*E_**2*lam_) == 0
roots = sp.solve(squared, x_)
assert any(sp.simplify(r - xp_eq5) == 0 for r in roots)
assert any(sp.simplify(r - xm_eq5) == 0 for r in roots)
print("   PASS: the quadratic has coefficients (4s, -4(E+m)(s+M^2-m^2), .), discriminant 16E^2 lambda,")
print("         and its two roots are exactly x_+ and x_- of eq (5)")

# ---------------------------------------------------------------------------
section("3. Which root is which: forward/backward, and the forward cone")
m = 0.51099895
def xpm(E, M):
    s = m*m + 2*m*E
    lam = (s - (m + M)**2)*(s - (m - M)**2)
    lam = max(lam, 0.0)
    base = (E + m)*(s + M*M - m*m)
    r = E*np.sqrt(lam)
    return (base - r)/(2*s), (base + r)/(2*s)

def cos_theta(E, x, M):
    return ((x - m)*E + m*x - M*M/2)/(E*np.sqrt(x*x - M*M))

# Whether x_- is backward emission depends on whether the A' outruns the CM:
#   beta*_A' = sqrt(lambda)/(s+M^2-m^2)   vs   beta_cm = E/(E+m).
# Sympy reduces beta*_A' > beta_cm to a linear condition in E:
beta_star = sp.sqrt(lam_)/(s_ + M_**2 - m_**2)
beta_cm = E_/(E_ + m_)
cond = sp.simplify(sp.expand(((beta_star**2 - beta_cm**2)*(E_ + m_)**2*(s_ + M_**2 - m_**2)**2)))
# cond > 0  <=>  (2mE + M^2)^2 > 4 M^2 (E+m)^2  <=>  2E(m-M) > M(2m-M)
assert sp.simplify(cond - s_*((2*m_*E_ + M_**2)**2 - 4*M_**2*(E_ + m_)**2)) == 0
assert sp.factor(sp.simplify((2*m_*E_ + M_**2)**2 - 4*M_**2*(E_+m_)**2)
                 - (2*m_*E_ + M_**2 - 2*M_*(E_+m_))*(2*m_*E_ + M_**2 + 2*M_*(E_+m_))) == 0
print("   PASS (sympy): beta*_A' > beta_cm  <=>  2E(m-M) > M(2m-M),")
print("         i.e. only for M < m and E > E_b = M(2m-M)/(2(m-M)); for M >= m, never.")

def E_back(M):
    return M*(2*m - M)/(2*(m - M)) if M < m else np.inf

grid = [(E, M) for M in (0.01, 0.1, 0.5, 1.0, 1.4) for E in np.linspace(1.001, 60.0, 60)
        if E > M + M*M/(2*m)]
worst = 0.0
n_back = n_cone = 0
for E, M in grid:
    lo, hi = xpm(E, M)
    worst = max(worst, abs(cos_theta(E, hi, M) - 1))            # x_+ is always theta = 0
    Ap = M*M + 2*m*E - 2*(m + E)*hi
    Am = M*M + 2*m*E - 2*(m + E)*lo
    assert Ap < 0, (E, M, Ap)                                    # x_+ always on the cos = +1 branch
    if E > E_back(M):
        assert Am > 0, (E, M, Am)
        worst = max(worst, abs(cos_theta(E, lo, M) + 1)); n_back += 1
    else:
        assert Am < 0, (E, M, Am)
        worst = max(worst, abs(cos_theta(E, lo, M) - 1)); n_cone += 1
assert worst < 1e-7   # cancellation in |k| near E_b, where x_- turns from forward to backward
OUT["xpm_costheta_residual"] = worst
OUT["E_backward_M0.5"] = E_back(0.5)
print(f"   PASS: |cos(theta)| = 1 at both endpoints to {worst:.1e} over {len(grid)} points")
print(f"         ({n_back} points with x_- backward, {n_cone} with the emission confined to a forward cone)")
print(f"         E_b(M=0.1) = {E_back(0.1):.4f} MeV (vs E_th = {0.1+0.01/(2*m):.4f});  E_b(M=0.5) = {E_back(0.5):.3f} MeV")

# In the cone regime cos(theta) has an interior minimum: that is the cone half-angle.
for (E, M) in ((3.0, 0.5), (3.0, 1.0)):
    lo, hi = xpm(E, M)
    xs = np.linspace(lo, hi, 200001)[1:-1]
    cmin = cos_theta(E, xs, M).min()
    OUT[f"theta_max_deg_E{E}_M{M}"] = float(np.degrees(np.arccos(cmin)))
    print(f"   E={E} MeV, M={M} MeV: emission confined to theta <= {np.degrees(np.arccos(cmin)):.2f} deg")

# outside the interval the exact relation refuses: |cos| > 1
bad = []
for E, M in grid:
    lo, hi = xpm(E, M)
    for x in (M + 0.3*(lo - M), lo - 1e-3*(hi - lo), hi + 1e-3*(hi - lo), hi*1.2):
        if x <= M:
            continue
        if abs(cos_theta(E, x, M)) <= 1:
            bad.append((E, M, x, cos_theta(E, x, M)))
assert not bad, bad[:5]
print(f"   PASS: |cos(theta)| > 1 at every sampled x outside [x_-, x_+] ({len(grid)*4} probes)")

# ---------------------------------------------------------------------------
section("4. The M -> 0 limit is the note's own range (and Compton)")
lim_p = sp.limit(xp_eq5, M_, 0)
lim_m = sp.limit(xm_eq5, M_, 0)
assert sp.simplify(lim_p - E_) == 0
assert sp.simplify(lim_m - E_/(1 + 2*E_/m_)) == 0
print("   PASS: x_+ -> E  and  x_- -> E/(1+2E/m), the Compton backscatter energy = note eq (25)")
print("         so the note's range is the M=0 limit of eq (5), not an approximation of it")

# ---------------------------------------------------------------------------
section("5. Threshold: lambda >= 0")
thr = sp.solve(sp.Eq(lam_, 0), E_)
Eth = [t for t in thr if sp.simplify(t - (M_ + M_**2/(2*m_))) == 0]
assert Eth, thr
print("   PASS: lambda(s,m^2,M^2) = 0  <=>  E = M + M^2/(2m)  (the other root is E = -M + M^2/2m < 0)")
# at threshold the two endpoints merge at the value the CM at rest gives
x_at_thr = sp.simplify(xp_eq5.subs(E_, M_ + M_**2/(2*m_)))
assert sp.simplify(x_at_thr - M_*(M_ + M_**2/(2*m_) + m_)/(m_ + M_)) == 0
print("   PASS: at threshold x_+ = x_- = M(E_th+m)/(m+M)  (A' at rest in the CM)")
for M in (0.1, 0.5, 1.0):
    Eth_n = M + M*M/(2*m)
    lo, hi = xpm(Eth_n, M)
    OUT[f"x_at_threshold_M{M}"] = hi
    print(f"     M={M}: E_th={Eth_n:.4f} MeV, x_-=x_+={hi:.4f} MeV")

# ---------------------------------------------------------------------------
section("6. The u-channel pole x = M^2/2m never lies inside [x_-, x_+]")
worst = 1e9
for M in (0.05, 0.1, 0.3, 0.5, 0.8, 1.0, 1.4):
    pole = M*M/(2*m)
    Es = np.linspace(M + M*M/(2*m), 60.0, 4000)[1:]
    ratios = np.array([xpm(E, M)[0] for E in Es])/pole
    worst = min(worst, ratios.min())
    OUT[f"min_xminus_over_pole_M{M}"] = float(ratios.min())
assert worst > 1
OUT["min_xminus_over_pole_all"] = float(worst)
print(f"   PASS: min over M,E of x_-/(M^2/2m) = {worst:.3f} > 1  (E up to 60 MeV, M up to 1.4 MeV)")

# ---------------------------------------------------------------------------
section("7. Numbers quoted in the write-up")
for M in (0.1, 0.5, 1.0):
    lo, hi = xpm(3.0, M)
    OUT[f"xrange_E3_M{M}"] = [lo, hi]
    print(f"   E=3 MeV, M={M}: [x_-, x_+] = [{lo:.4f}, {hi:.4f}] MeV")
lo1, hi1 = xpm(3.0, 1.0)
OUT["xminus_E3_M1"] = lo1
# the note's own range at the same point, for the contrast
note_lo = 3.0/(1 + 2*3.0/m)
OUT["note_xminus_E3"] = note_lo
print(f"   note's range at E=3 MeV (its eq (24)-(25), M-independent): [{note_lo:.4f}, 3.0000] MeV")

print("\nALL ASSERTIONS PASSED")
if len(sys.argv) > 1:
    with open(sys.argv[1], "w") as f:
        json.dump(OUT, f, indent=1)
    print(f"wrote {sys.argv[1]}")
