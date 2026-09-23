#!/usr/bin/env python3
"""Independent check of the dark-Compton cross section  gamma(k) + e(p) -> A'(q) + e(p').

Every block ends in an assert. What is derived here from scratch:
  * the squared amplitude, from explicit 4x4 Dirac matrices and explicit
    polarization vectors (two for the photon, three for the massive A'),
    evaluated at random on-shell lab-frame points;
  * the exact lab kinematics (threshold, x_min, x_max, cos(theta)(x));
  * the transverse / longitudinal content of the produced A';
  * the ratio sigma(A'e->gamma e)/sigma_KN for a transverse vs an unpolarized A'.
What is recalled (textbook, not derived here): dsigma/dt = |M|^2 / (16 pi lambda),
Tr identities are NOT used (the matrices are multiplied numerically), Klein-Nishina
total, and the spin-completeness sum u ubar = pslash + m.

Formulas under test are transcribed from the documents in apuntes/ and papers_reactor_dark_photon/
and labelled by source and equation number.

Author: Eliana L. De Paoli.
Part of "Sensitivity to dark photons in a reactor antineutrino experiment with Skipper CCDs".
Written with agent assistance; provenance/ records what was derived, recalled and verified.
"""
import json
import sys
import numpy as np
import sympy as sp
from scipy.integrate import quad

rng = np.random.default_rng(20260921)
OUT = {}

# ----------------------------------------------------------------------------
# 0. Dirac matrices (Dirac representation), metric (+,-,-,-)
# ----------------------------------------------------------------------------
I2 = np.eye(2); Z2 = np.zeros((2, 2))
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], complex)
g0 = np.block([[I2, Z2], [Z2, -I2]]).astype(complex)
gs = [np.block([[Z2, s], [-s, Z2]]) for s in (sx, sy, sz)]
GAM = [g0] + gs                      # gamma^mu, upper index
METRIC = np.diag([1., -1., -1., -1.])
GAM_LOW = [METRIC[m, m] * GAM[m] for m in range(4)]   # gamma_mu
I4 = np.eye(4, dtype=complex)


def slash(v):
    return sum(GAM_LOW[m] * v[m] for m in range(4))


def mdot(a, b):
    return a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3]


def bar(M):
    return g0 @ M.conj().T @ g0


# ----------------------------------------------------------------------------
# 1. Lab kinematics, exact.  p=(m,0), k=(E,0,0,E), q=(x, Q sin, 0, Q cos), Q=sqrt(x^2-M^2)
# ----------------------------------------------------------------------------
def kallen(s, a, b):
    return (s - (np.sqrt(a) + np.sqrt(b))**2) * (s - (np.sqrt(a) - np.sqrt(b))**2)


def x_pm(E, m, M):
    """Exact endpoints of the A' lab energy (derived: boost of the CM two-body momentum)."""
    s = m*m + 2*m*E
    lam = kallen(s, m*m, M*M)
    root = E*np.sqrt(lam)
    base = (E + m)*(s + M*M - m*m)
    return (base - root)/(2*s), (base + root)/(2*s)


def cos_theta(E, x, m, M):
    """Derived from (p+k-q)^2 = m^2:  M^2 + 2m(E-x) - 2E(x - Q cos) = 0."""
    Q = np.sqrt(x*x - M*M)
    return ((x - m)*E + m*x - M*M/2)/(E*Q)


def momenta(E, x, m, M):
    Q = np.sqrt(x*x - M*M)
    c = cos_theta(E, x, m, M)
    s_ = np.sqrt(1 - c*c)
    p = np.array([m, 0, 0, 0.])
    k = np.array([E, 0, 0, E])
    q = np.array([x, Q*s_, 0, Q*c])
    pp = p + k - q
    return p, k, q, pp


def polarizations_A(q, M):
    """Three orthonormal polarization vectors of a massive vector with momentum q in the x-z plane."""
    x, Q = q[0], np.hypot(q[1], q[3])
    s_, c = q[1]/Q, q[3]/Q
    e1 = np.array([0, c, 0, -s_])
    e2 = np.array([0, 0, 1, 0.])
    eL = np.array([Q, x*s_, 0, x*c])/M
    return [e1, e2, eL]


# ----------------------------------------------------------------------------
# 2. Squared amplitude from explicit matrices.  M = eps e^2 ubar(p') Gamma^{nu mu} u(p) eps_mu(k) xi*_nu(q)
#    Gamma^{nu mu} = gamma^nu (rs+m) gamma^mu / a + gamma^mu (ru+m) gamma^nu / b,  rs=p+k, ru=p-q
#    The factor eps^2 e^4 is stripped: we compute F_num = (1/4) sum_pol |ubar Gamma u|^2
# ----------------------------------------------------------------------------
def vertex(p, k, q, m, eps_vec, xi_vec):
    rs, ru = p + k, p - q
    a = mdot(rs, rs) - m*m
    b = mdot(ru, ru) - m*m
    xs, es = slash(xi_vec), slash(eps_vec)
    return (xs @ (slash(rs) + m*I4) @ es)/a + (es @ (slash(ru) + m*I4) @ xs)/b


def F_numeric(E, x, m, M, which="all"):
    """(1/4) sum over spins and over the listed A' polarizations of |M|^2 / (eps e^2)^2."""
    p, k, q, pp = momenta(E, x, m, M)
    eps_ph = [np.array([0, 1, 0, 0.]), np.array([0, 0, 1, 0.])]
    xis = polarizations_A(q, M)
    sel = {"all": xis, "T": xis[:2], "L": xis[2:]}[which]
    tot = 0.
    for ep in eps_ph:
        for xi in sel:
            G = vertex(p, k, q, m, ep, np.conj(xi))
            # sum over spins:  Tr[(p'+m) G (p+m) Gbar]
            tot += np.trace((slash(pp) + m*I4) @ G @ (slash(p) + m*I4) @ bar(G)).real
    return tot/4


def F_covariant(E, x, m, M):
    """Same trace, but with the polarization sums replaced by -g_{mu mu'} -g_{nu nu'} (Ward identity assumed)."""
    p, k, q, pp = momenta(E, x, m, M)
    rs, ru = p + k, p - q
    a = mdot(rs, rs) - m*m
    b = mdot(ru, ru) - m*m
    tot = 0.
    for mu in range(4):
        for nu in range(4):
            G = (GAM[nu] @ (slash(rs) + m*I4) @ GAM[mu])/a + (GAM[mu] @ (slash(ru) + m*I4) @ GAM[nu])/b
            Gl = (GAM_LOW[nu] @ (slash(rs) + m*I4) @ GAM_LOW[mu])/a + (GAM_LOW[mu] @ (slash(ru) + m*I4) @ GAM_LOW[nu])/b
            tot += np.trace((slash(pp) + m*I4) @ G @ (slash(p) + m*I4) @ bar(Gl)).real
    return tot/4


# ----------------------------------------------------------------------------
# 3. Published / note formulas, transcribed.  a = s-m^2, b = u-m^2, c = m^2, d = M^2
# ----------------------------------------------------------------------------
def F_v2(a, b, c, d):
    """apuntes/report_park_production_v2_en.pdf eq (36)."""
    S = 1/a + 1/b
    return -2*(a/b + b/a) + 4*(2*c + d)*(S + c*S**2 - d/(a*b))


def F_apuntes(a, b, c, d):
    """apuntes/Dark_Photon_cross_section-1.pdf eq (40) = Dark_Photon_cross_Section.pdf eq (41),
    |M|^2 = 2 eps^2 e^4 [ -u/s - s/u + 2M^2(1/s+1/u) - M^4 (1/s+1/u)^2 su/4 ] with s,u the propagator denominators.
    Returned in the same normalisation as F (i.e. divided by eps^2 e^4)."""
    S = 1/a + 1/b
    return 2*(-b/a - a/b + 2*d*S - d*d*S**2*a*b/4)


def M2_deNiverville(s, u, m, M):
    """papers 03 deNiverville-Lee-Lee eq (18)-(20): |M|^2 = 32 pi^2 alpha^2 eps^2 (A+B)/((m^2-s)^2 (m^2-u)^2).
    Returned divided by (eps e^2)^2 = 16 pi^2 alpha^2 eps^2, i.e. 2(A+B)/(...)."""
    A = 6*m**8 - 2*M**4*(m*m - s)*(m*m - u) - s*u*(s*s + u*u) + m*m*(s + u)*(s*s + 6*s*u + u*u)
    B = -m**4*(3*s*s + 14*s*u + 3*u*u) + 2*M*M*(-4*m*m*s*u + m**4*(s + u) + s*u*(s + u))
    return 2*(A + B)/((m*m - s)**2*(m*m - u)**2)


def dsig_dEX_smirnov(E, x, m, M, alpha):
    """papers 04 Smirnov et al eq (3) with Z=1, g_X = eps e (eps stripped), then dsigma/dE_X = (1/E) dsigma/dx_S.
    x_S = 1 - E_X/E + M^2/(2 E m)."""
    s = m*m + 2*E*m
    xS = 1 - x/E + M*M/(2*E*m)
    gX2 = 4*np.pi*alpha
    pre = alpha*gX2/(2*(s - m*m)**3*(1 - xS))
    br = (-2*m*m*s*(xS*xS + 2) + 2*M**4 - 2*M*M*s*xS + s*s*(xS*xS - 2*xS + 2)
          - (m**4*(xS**3 - 3*xS*xS - 2) + 2*m*m*M*M*(xS - 2))/(1 - xS))
    return pre*br/E


def dsig_dx_exact(E, x, m, M, alpha):
    """dsigma/dx = eps^2 e^4 F /(32 pi m E^2)  (eps stripped).  Recalled: dsigma/dt=|M|^2/(16 pi (s-m^2)^2); derived: dt/dx = 2m."""
    a = 2*m*E
    b = M*M - 2*m*x
    return (4*np.pi*alpha)**2*F_v2(a, b, m*m, M*M)/(32*np.pi*m*E*E)


def sigma_KN(E, m, alpha):
    y = E/m
    re = alpha/m
    return 2*np.pi*re*re*((1 + y)/y**2*(2*(1 + y)/(1 + 2*y) - np.log(1 + 2*y)/y) + np.log(1 + 2*y)/(2*y) - (1 + 3*y)/(1 + 2*y)**2)


# ----------------------------------------------------------------------------
# 4. Run the checks
# ----------------------------------------------------------------------------
m = 0.51099895   # MeV
alpha = 1/137.035999


def dsig_dx_apuntes(E, x, m, M, alpha):
    """doc1 eq (50): pi eps^2 alpha^2/(m w^2) (1 + M^2/(2 m w)) [w'/w + w/w' - sin^2 + M^2 f], f = 2/(s+u)(1 - m^2 M^2/(s u)),
    s,u the propagator denominators, cos(theta) from eq (20)."""
    c = 1 - m*(1/x - 1/E) - M*M/(2*E*x)
    a, b = 2*m*E, M*M - 2*m*x
    f = 2/(a + b)*(1 - m*m*M*M/(a*b))
    return np.pi*alpha**2/(m*E*E)*(1 + M*M/(2*m*E))*(x/E + E/x - (1 - c*c) + M*M*f)

def dsig_dx_doc2(E, x, m, M, alpha):
    """doc2 eq (54): same bracket, prefactor pi eps^2 alpha^2/(m w^2) without the (1+M^2/2mw) factor."""
    c = 1 - m*(1/x - 1/E) - M*M/(2*E*x)
    a, b = 2*m*E, M*M - 2*m*x
    f = 2/(a + b)*(1 - m*m*M*M/(a*b))
    return np.pi*alpha**2/(m*E*E)*(x/E + E/x - (1 - c*c) + M*M*f)



# ----------------------------------------------------------------------------
# Reactor spectra (Park eq 3, P = 1 GW, eps = 1)
# ----------------------------------------------------------------------------
def S_park(E, P_MW=1000.0):
    """Park eq (3): photons per s per MeV."""
    return 0.58e18*P_MW*np.exp(-E/0.91)


def dN_dx_compton(x, M, Emax=40.0):
    """Model A: Park eq (1) with sigma_tot = sigma_KN and the exact dsigma/dx, per s per MeV at 1 GW, eps=1.
    Only photons with x_-(E) <= x <= x_+(E) contribute: E runs from E_L (x_+(E_L) = x) to E_U (x_-(E_U) = x, or Emax)."""
    from scipy.optimize import brentq
    Eth = M + M*M/(2*m)
    E0 = Eth*(1 + 1e-9)
    if x_pm(Emax, m, M)[1] <= x:
        return 0.0
    EL = brentq(lambda E: x_pm(E, m, M)[1] - x, E0, Emax) if x_pm(E0, m, M)[1] < x else E0
    g = lambda E: x_pm(E, m, M)[0] - x
    EU = brentq(g, EL, Emax) if g(EL) < 0 < g(Emax) else (Emax if g(Emax) < 0 else EL)
    if EU <= EL:
        return 0.0
    f = lambda E: S_park(E)/sigma_KN(E, m, alpha)*dsig_dx_exact(E, x, m, M, alpha)
    return quad(f, EL, EU, limit=200)[0]


def dN_dx_oscillation(x, M):
    """Model B: Danilov eq (6) in the limit m_X >> m_gamma, m_X^2 >> E Gamma: P = eps^2, E_X = E_gamma, hard cutoff x > M."""
    return S_park(x) if x > M else 0.0


def main():

    print("== 1. explicit-matrix |M|^2 vs published formulas at random on-shell points ==")
    worst = {"ward": 0, "v2": 0, "deNiv": 0, "smirnov": 0}
    long_frac = []
    for _ in range(60):
        M = rng.uniform(0.05, 1.5)
        Eth = M + M*M/(2*m)
        E = rng.uniform(Eth*1.02, 10)
        xlo, xhi = x_pm(E, m, M)
        x = rng.uniform(xlo + 1e-6*(xhi - xlo), xhi - 1e-6*(xhi - xlo))
        p, k, q, pp = momenta(E, x, m, M)
        assert abs(mdot(pp, pp) - m*m) < 1e-9 and abs(mdot(q, q) - M*M) < 1e-9
        Fn = F_numeric(E, x, m, M)
        Fc = F_covariant(E, x, m, M)
        a, b = 2*m*E, M*M - 2*m*x
        s, u = m*m + a, m*m + b
        worst["ward"] = max(worst["ward"], abs(Fn/Fc - 1))
        worst["v2"] = max(worst["v2"], abs(Fn/F_v2(a, b, m*m, M*M) - 1))
        worst["deNiv"] = max(worst["deNiv"], abs(Fn/M2_deNiverville(s, u, m, M) - 1))
        worst["smirnov"] = max(worst["smirnov"], abs(dsig_dx_exact(E, x, m, M, alpha)/dsig_dEX_smirnov(E, x, m, M, alpha) - 1))
        long_frac.append((E, M, F_numeric(E, x, m, M, "L")/Fn))
    for kk, v in worst.items():
        print(f"   max |ratio-1| vs {kk:8s}: {v:.2e}")
        assert v < 1e-9, kk
    OUT["max_rel_dev_matrix_vs_v2_F"] = worst["v2"]
    OUT["max_rel_dev_matrix_vs_deNiverville"] = worst["deNiv"]
    OUT["max_rel_dev_exact_vs_smirnov"] = worst["smirnov"]
    OUT["max_rel_dev_explicit_pol_vs_covariant"] = worst["ward"]
    print("   PASS: explicit-polarization sum == covariant (-g) sum: Ward identity holds numerically")
    print("   PASS: v2 eq(36), deNiverville eq(18-20) and Smirnov eq(3) are the same cross section")

    print("== 2. the Spanish notes' |M|^2 (eq 40/41) against the same points ==")
    devs = []
    for _ in range(60):
        M = rng.uniform(0.05, 1.5); Eth = M + M*M/(2*m); E = rng.uniform(Eth*1.02, 10)
        xlo, xhi = x_pm(E, m, M); x = rng.uniform(xlo*1.001, xhi*0.999)
        a, b = 2*m*E, M*M - 2*m*x
        devs.append(F_apuntes(a, b, m*m, M*M)/F_v2(a, b, m*m, M*M) - 1)
    devs = np.array(devs)
    print(f"   ratio apuntes/exact - 1 : min {devs.min():+.3f}, max {devs.max():+.3f}")
    assert np.abs(devs).max() > 0.05
    OUT["apuntes_M2_rel_dev_min"] = float(devs.min()); OUT["apuntes_M2_rel_dev_max"] = float(devs.max())
    # symbolic: difference at M=0 is exactly the electron-mass terms of Klein-Nishina
    A, B, C, D = sp.symbols('a b c d')
    Sb = 1/A + 1/B
    Fv2 = -2*(A/B + B/A) + 4*(2*C + D)*(Sb + C*Sb**2 - D/(A*B))
    Fap = 2*(-B/A - A/B + 2*D*Sb - D**2*Sb**2*A*B/4)
    diff0 = sp.simplify((Fv2 - Fap).subs(D, 0) - (8*C*Sb + 8*C**2*Sb**2))
    assert diff0 == 0
    print("   PASS (sympy): at M=0, exact - apuntes = 8 m^2 (1/a+1/b) + 8 m^4 (1/a+1/b)^2  -> the notes drop the m_e terms of Klein-Nishina")
    # KN limit of the exact formula in lab variables
    E_, x_, th = sp.symbols('E x theta', positive=True)
    FKN_lab = 2*(E_/x_ + x_/E_ - sp.sin(th)**2)
    xKN = E_/(1 + E_/sp.sqrt(C)*(1 - sp.cos(th)))
    FKN_inv = Fv2.subs(D, 0).subs({A: 2*sp.sqrt(C)*E_, B: -2*sp.sqrt(C)*x_})
    assert sp.simplify((FKN_inv - FKN_lab).subs(x_, xKN)) == 0
    print("   PASS (sympy): exact F at M=0 == 2[w'/w + w/w' - sin^2] (Klein-Nishina) with the Compton relation")

    print("== 3. kinematics claimed in the notes ==")
    Ms = [0.1, 0.5, 1.0]
    thr = {M: M + M*M/(2*m) for M in Ms}
    print("   threshold E_th = M + M^2/2m :", {M: round(v, 6) for M, v in thr.items()},
          "   (notes say 'omega -> m_A'')")
    OUT["Eth_M1MeV"] = thr[1.0]
    E = 3.0; M = 1.0
    xlo, xhi = x_pm(E, m, M)
    # doc1 eq (24)-(25) and doc2 eq (58)-(59)
    d1max, d1min = E, (m*E + M*M/2)/(2*E + m)
    d2max, d2min = E/(1 + M*M/(2*m*E)), E/(1 + 2*E/m + M*M/(2*m*E))
    print(f"   E=3 MeV, M=1 MeV: exact x in [{xlo:.4f},{xhi:.4f}] ; doc1 [{d1min:.4f},{d1max:.4f}] ; doc2 [{d2min:.4f},{d2max:.4f}]")
    OUT["xrange_exact_E3_M1"] = [xlo, xhi]; OUT["xrange_doc1_E3_M1"] = [d1min, d1max]; OUT["xrange_doc2_E3_M1"] = [d2min, d2max]
    assert abs(xhi - E) > 1e-3 and abs(d2max - xhi) > 1e-3
    # doc1 eq (19) vs doc2 eq (20): both claim to solve the same approximate equation (18); they differ
    w, wp, M_, c_ = sp.symbols('omega omegap M costh', positive=True)
    eq18 = M_**2 + 2*sp.sqrt(C)*(w - wp) - 2*w*wp*(1 - c_)
    sol = sp.solve(eq18, wp)[0]
    doc1 = (w + M_**2/(2*sp.sqrt(C)))/(1 + w/sp.sqrt(C)*(1 - c_))
    doc2 = w/(1 + w/sp.sqrt(C)*(1 - c_) + M_**2/(2*sp.sqrt(C)*w))
    assert sp.simplify(sol - doc1) == 0 and sp.simplify(sol - doc2) != 0
    print("   PASS (sympy): doc1 eq(19) solves its eq(18); doc2 eq(20) does NOT solve the same equation (and doc2 eq(19) has the sign of M^2 flipped)")
    # the exact relation has |q| = sqrt(x^2-M^2), not x.  Size of the |q|~x approximation:
    c_exact = cos_theta(E, 2.0, m, M)
    c_doc1 = 1 - m*(1/2.0 - 1/E) - M*M/(2*E*2.0)
    print(f"   E=3, x=2, M=1: cos(theta) exact {c_exact:.4f}, doc1 eq(20) {c_doc1:.4f}")
    OUT["costheta_exact_E3_x2_M1"] = c_exact; OUT["costheta_doc1_E3_x2_M1"] = c_doc1

    print("== 4. the notes' final dsigma/domega' (doc1 eq 50 / doc2 eq 54) against the exact one ==")
    tab = []
    for M in Ms:
        E = 3.0
        xlo, xhi = x_pm(E, m, M)
        xs = np.linspace(xlo + 1e-4, xhi - 1e-4, 7)
        r1 = dsig_dx_apuntes(E, xs, m, M, alpha)/dsig_dx_exact(E, xs, m, M, alpha)
        r2 = dsig_dx_doc2(E, xs, m, M, alpha)/dsig_dx_exact(E, xs, m, M, alpha)
        tab.append((M, r1.min(), r1.max(), r2.min(), r2.max()))
        print(f"   E=3 MeV M={M}: doc1/exact in [{r1.min():.3f},{r1.max():.3f}] ; doc2/exact in [{r2.min():.3f},{r2.max():.3f}]")
    OUT["doc_over_exact_E3"] = tab
    assert abs(tab[0][1] - 1) < 0.05 and abs(tab[2][2] - 1) > 0.05  # fine at 0.1 MeV, off at 1 MeV
    # also the total integrated over the (wrong) limits
    tot = []
    for M in Ms:
        E = 3.0
        xlo, xhi = x_pm(E, m, M)
        ex = quad(lambda x: dsig_dx_exact(E, x, m, M, alpha), xlo, xhi)[0]
        d1 = quad(lambda x: dsig_dx_apuntes(E, x, m, M, alpha), xlo, xhi)[0]
        pole = M*M/(2*m)
        inside = (m*E + M*M/2)/(2*E + m) < pole < E
        tot.append((M, d1/ex, bool(inside)))
        print(f"   E=3 MeV M={M}: sigma_doc1(over exact range)/sigma_exact = {d1/ex:.3f}; pole u-m^2=0 at E_A'={pole:.3f} inside doc1's own range: {inside}")
    OUT["sigma_doc1_over_exact_E3"] = tot

    print("== 5. totals ==")
    # exact total at M=0 equals Klein-Nishina
    for E in (1., 3., 8.):
        xlo, xhi = x_pm(E, m, 0.0)
        num = quad(lambda x: dsig_dx_exact(E, x, m, 1e-9, alpha), xlo, xhi)[0]
        assert abs(num/sigma_KN(E, m, alpha) - 1) < 1e-6
    print("   PASS: quadrature of exact dsigma/dx at M->0 == Klein-Nishina total (E = 1, 3, 8 MeV)")
    # v2 eq (41)-(42): dG/db = F symbolically
    H = 2*C + D
    G = (-B**2/A + 4*H*(1/A + C/A**2)*B + (-2*A + 4*H*(1 + (2*C - D)/A))*sp.log(-B/A) - 4*H*C/B)
    assert sp.simplify(sp.diff(G, B) - Fv2) == 0
    print("   PASS (sympy): v2 eq(41) primitive satisfies dG/db = F")
    # doc1 eq (52) vs Klein-Nishina eq (53)
    kap = sp.symbols('kappa', positive=True)
    xx = sp.symbols('xx', positive=True)
    integrand52 = xx + (1 - 2*kap - 2*kap**2)/xx + kap**2/xx**2 + kap*(2 + kap)
    I52 = sp.integrate(integrand52, (xx, kap/(2 + kap), 1))
    claimed52 = sp.Rational(1, 2) + 4*kap - kap**2/(2*(2 + kap)**2) + (1 - 2*kap - 2*kap**2)*sp.log(1 + 2/kap)
    # doc1 says sigma = pi eps^2 alpha^2/(m w) * [claimed52]; KN with x=1/kappa, r_e=alpha/m:
    y = 1/kap
    KN = 2*sp.pi*(1/sp.sqrt(C))**2*((1 + y)/y**2*(2*(1 + y)/(1 + 2*y) - sp.log(1 + 2*y)/y) + sp.log(1 + 2*y)/(2*y) - (1 + 3*y)/(1 + 2*y)**2)  # /alpha^2
    doc1_tot = sp.pi/(sp.sqrt(C)*(sp.sqrt(C)/kap))*claimed52   # /alpha^2, omega = m/kappa
    d_int = sp.simplify(I52 - claimed52)
    d_kn = [float((doc1_tot - KN).subs({C: 1, kap: kv})/KN.subs({C: 1, kap: kv})) for kv in (0.1, 0.5, 2.0)]
    print(f"   doc1 eq(52): integral of its integrand - its printed result = {d_int}")
    print(f"   doc1 eq(52) vs KN eq(53): relative difference at kappa=0.1,0.5,2 = {[round(v,4) for v in d_kn]}")
    OUT["doc1_eq52_vs_KN_reldiff"] = d_kn
    # doc2 eq (64) asymptotic vs KN asymptotic  sigma_KN -> pi r_e^2 / y (ln 2y + 1/2)
    yv = 200.
    kn_as = np.pi*(alpha/m)**2/yv*(np.log(2*yv) + 0.5)
    doc2_as = alpha**2/(m*m*yv)*(2*np.log(2*yv) + 0.5)
    print(f"   y=E/m=200: KN exact/asymptotic = {sigma_KN(yv*m, m, alpha)/kn_as:.4f}; doc2 eq(64)/KN exact = {doc2_as/sigma_KN(yv*m, m, alpha):.4f}")
    OUT["doc2_eq64_over_KN_y200"] = doc2_as/sigma_KN(yv*m, m, alpha)

    print("== 6. polarization of the produced A' and the 2/3 detection factor ==")
    lf = np.array([[e, M, f] for e, M, f in long_frac])
    print(f"   longitudinal fraction of produced A' over the random sample: max {lf[:,2].max():.3e}")
    for E, M in ((3., 0.1), (3., 0.5), (3., 1.0), (8., 1.0)):
        xlo, xhi = x_pm(E, m, M)
        fT = quad(lambda x: F_numeric(E, x, m, M, "T"), xlo, xhi)[0]
        fL = quad(lambda x: F_numeric(E, x, m, M, "L"), xlo, xhi)[0]
        print(f"   E={E} M={M}: sigma_L/sigma_tot = {fL/(fT+fL):.4f}")
        OUT[f"long_frac_E{E}_M{M}"] = fL/(fT + fL)

    # inverse process A'(q)+e(p) -> gamma(k)+e(p'):  same vertex with the roles crossed.
    def inverse_sigma(EA, m, M, alpha, which):
        """sigma(A' e -> gamma e) for a transverse ('T'), longitudinal ('L') or unpolarized ('all', average over 3) A'.
        Lab: p=(m,0), q=(EA,0,0,Q).  Integrate over the CM angle with dsigma/dOmega* = |M|^2 /(64 pi^2 s) * (k*/q*)."""
        s = m*m + M*M + 2*m*EA
        Q = np.sqrt(EA*EA - M*M)
        qstar = np.sqrt(kallen(s, m*m, M*M))/(2*np.sqrt(s))
        kstar = (s - m*m)/(2*np.sqrt(s))
        p = np.array([m, 0, 0, 0.]); q = np.array([EA, 0, 0, Q])
        xis = [np.array([0, 1, 0, 0.]), np.array([0, 0, 1, 0.]), np.array([Q, 0, 0, EA])/M]
        sel = {"T": xis[:2], "L": xis[2:], "all": xis}[which]
        navg = {"T": 2, "L": 1, "all": 3}[which]

        def integrand(cth):
            sth = np.sqrt(1 - cth*cth)
            # photon in CM, boost to lab along z
            beta = Q/(EA + m); gam = (EA + m)/np.sqrt(s)
            kst = np.array([kstar, kstar*sth, 0, kstar*cth])
            k = np.array([gam*(kst[0] + beta*kst[3]), kst[1], 0, gam*(kst[3] + beta*kst[0])])
            pp = p + q - k
            assert abs(mdot(pp, pp) - m*m) < 1e-8
            # photon polarizations orthogonal to k
            kh = k[1:]/np.linalg.norm(k[1:])
            e1 = np.cross(kh, [0, 1, 0]); e1 /= np.linalg.norm(e1)
            e2 = np.cross(kh, e1)
            tot = 0.
            for ep in (np.r_[0, e1], np.r_[0, e2]):
                for xi in sel:
                    # incoming A' polarization xi (not conjugated), outgoing photon eps* (real here)
                    rs, ru = p + q, p - k
                    a = mdot(rs, rs) - m*m; b = mdot(ru, ru) - m*m
                    G = (slash(ep) @ (slash(rs) + m*I4) @ slash(xi))/a + (slash(xi) @ (slash(ru) + m*I4) @ slash(ep))/b
                    tot += np.trace((slash(pp) + m*I4) @ G @ (slash(p) + m*I4) @ bar(G)).real
            M2 = tot/(2*navg)   # average over electron spin (2) and A' polarizations
            return M2*(4*np.pi*alpha)**2/(64*np.pi**2*s)*(kstar/qstar)*2*np.pi
        return quad(integrand, -1, 1)[0]

    for EA, M in ((3., 0.01), (3., 0.1), (5., 0.1)):
        sT = inverse_sigma(EA, m, M, alpha, "T")
        sU = inverse_sigma(EA, m, M, alpha, "all")
        sL = inverse_sigma(EA, m, M, alpha, "L")
        kn = sigma_KN(EA, m, alpha)
        print(f"   E_A'={EA} M={M}: sigma_T/sigma_KN = {sT/kn:.5f}, sigma_unpol/sigma_KN = {sU/kn:.5f}, sigma_L/sigma_KN = {sL/kn:.2e}")
        OUT[f"inverse_T_over_KN_E{EA}_M{M}"] = sT/kn
        OUT[f"inverse_unpol_over_KN_E{EA}_M{M}"] = sU/kn
        if M == 0.01:
            assert abs(sT/kn - 1) < 1e-3 and abs(sU/kn - 2/3) < 1e-3
    print("   PASS: transverse A' -> sigma_KN (factor 1); unpolarized A' -> (2/3) sigma_KN. The 2/3 is the unpolarized average.")
    print("   polarization-resolved detection factors at the masses where the produced A' is not purely transverse:")
    for EA, M in ((3., 0.5), (5., 0.5), (8., 0.5), (3., 1.0), (5., 1.0), (8., 1.0)):
        sT = inverse_sigma(EA, m, M, alpha, "T"); sL = inverse_sigma(EA, m, M, alpha, "L"); kn = sigma_KN(EA, m, alpha)
        print(f"   E_A'={EA} M={M}: sigma_T/sigma_KN = {sT/kn:.4f}, sigma_L/sigma_KN = {sL/kn:.4f}, unpolarized/sigma_KN = {(2*sT+sL)/(3*kn):.4f}")
        OUT[f"inverse_T_over_KN_E{EA}_M{M}"] = sT/kn; OUT[f"inverse_L_over_KN_E{EA}_M{M}"] = sL/kn

    print("== 7. Park's approximation dsigma/dE_A' ~ eps^2 dsigma_C/dE_r (eq 2) ==")
    def dsigC_dEr(E, x, m, alpha):
        return dsig_dx_exact(E, x, m, 1e-9, alpha)
    for M in Ms:
        E = 3.0
        xlo, xhi = x_pm(E, m, M)
        ex = quad(lambda x: dsig_dx_exact(E, x, m, M, alpha), xlo, xhi)[0]
        kn = sigma_KN(E, m, alpha)
        print(f"   E=3 MeV M={M}: sigma_exact/sigma_KN = {ex/kn:.4f}")
        OUT[f"sigma_exact_over_KN_E3_M{M}"] = ex/kn

    print("== 8. approach A (Compton-like, Park eq 1 with sigma_tot = sigma_KN) vs approach B (oscillation, Danilov eq 6 with E_X = E_gamma) ==")
    # Park eq (3): dN/dE = 0.58e18 (P/MW) exp(-E/0.91) per s per MeV; P = 1 GW = 1000 MW
    S = lambda E: 0.58e18*1000*np.exp(-E/0.91)
    NB = quad(S, 3, 8)[0]                       # approach B, eps=1, m_X >> m_gamma, no absorption term: P = eps^2, E_X = E_gamma
    for M in (0.1, 0.5, 1.0):
        def inner(E):
            xlo, xhi = x_pm(E, m, M)
            lo, hi = max(3.0, xlo), min(8.0, xhi)
            if hi <= lo:
                return 0.0
            return S(E)/sigma_KN(E, m, alpha)*quad(lambda x: dsig_dx_exact(E, x, m, M, alpha), lo, hi)[0]
        NA = quad(inner, 3.0, 40.0, limit=200)[0]
        print(f"   M={M} MeV, 1 GW, eps=1, A' energies in [3,8] MeV: approach A = {NA:.3e} /s ; approach B = {NB:.3e} /s ; ratio A/B = {NA/NB:.3f}")
        OUT[f"N_window_A_M{M}"] = NA; OUT[f"N_window_B"] = NB; OUT[f"ratio_A_over_B_M{M}"] = NA/NB
    Ntot = quad(S, 1, np.inf)[0]
    print(f"   check: integral of Park's spectrum above 1 MeV at 1 GW = {Ntot:.4e} /s (Park prints 1.76e20)")
    assert abs(Ntot/1.76e20 - 1) < 0.01
    OUT["Ngamma_above_1MeV_1GW"] = Ntot

    print("== 9. dN/dE_A' at the reactor (1 GW, eps=1): model A vs model B vs Park's Figure 1 read by eye ==")
    park_fig1_by_eye = {1.0: 0.40, 2.0: 0.065, 3.0: 0.022, 4.0: 0.007}   # 1e21/MeV/s, M=0.1 MeV curve, read from the raster (recalled, +-30%)
    for x in (1.0, 2.0, 3.0, 4.0):
        A = dN_dx_compton(x, 0.1)/1e21; B = dN_dx_oscillation(x, 0.1)/1e21
        print(f"   x={x} MeV: A (eq 1, exact sigma) = {A:.4f} ; B (= S_gamma) = {B:.4f} ; Park Fig 1 ~ {park_fig1_by_eye[x]:.3f} ; Park/A = {park_fig1_by_eye[x]/A:.1f}, Park/B = {park_fig1_by_eye[x]/B:.2f}")
        OUT[f"dNdx_A_M0.1_x{x}"] = A; OUT[f"dNdx_B_x{x}"] = B
    # flat-KN sanity estimate of A at x=2: dN/dx ~ <r> * int_x^inf S dE with r = (1/sigma_KN) dsigma/dx ~ 0.15-0.42 /MeV
    est_lo, est_hi = 0.15*quad(S_park, 2, 40)[0]/1e21, 0.42*quad(S_park, 2, 40)[0]/1e21
    A2 = OUT["dNdx_A_M0.1_x2.0"]
    print(f"   sanity: flat-KN bounds at x=2: [{est_lo:.4f}, {est_hi:.4f}] contain A = {A2:.4f}: {est_lo < A2 < est_hi}")
    assert est_lo < A2 < est_hi
    assert park_fig1_by_eye[2.0]/A2 > 2.5 and abs(park_fig1_by_eye[2.0]/OUT["dNdx_B_x2.0"] - 1) < 0.3
    print("   Park's Fig 1 tail is the photon spectrum S_gamma(E), not the Compton-like convolution of his eq (1)")
    # normalisation check: integral of model A over x equals the number of photons (M -> 0)
    # photons below 0.02 MeV are not in the spectrum's stated domain; A' from them fall below x = 0.02 and are left out on both sides
    NA_tot = sum(quad(lambda x: dN_dx_compton(x, 1e-3), lo, hi, limit=200)[0] for lo, hi in ((0.02, 0.3), (0.3, 1.0), (1.0, 12.0)))
    Ng = quad(S_park, 0.02, 40)[0]
    print(f"   integral of A over x>0.02 MeV (M=1 keV) = {NA_tot:.4e} /s ; photons above 0.02 MeV = {Ng:.4e} /s ; ratio {NA_tot/Ng:.4f}")
    assert abs(NA_tot/Ng - 1) < 0.01
    OUT["NA_tot_over_Ngamma"] = NA_tot/Ng

    json.dump(OUT, open(sys.argv[1] if len(sys.argv) > 1 else "check_cross_section_out.json", "w"), indent=1, default=float)
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
