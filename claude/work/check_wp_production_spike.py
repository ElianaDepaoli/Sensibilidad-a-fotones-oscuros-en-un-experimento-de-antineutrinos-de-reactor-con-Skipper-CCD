#!/usr/bin/env python3
"""Why the Fig. 1 curves spike at w' -> mA' once wp_production_new_f is used with the OLD Rw_f.

Reproduces, in python, the three pieces of mycodes/DP3_Danilov.C as they stand after the user
took only the *_new_f kinematics: the new production grid (starts at w' = mA'), the old
integration domain Rw_f (built from invwpmin_f, the inverse of the APPROXIMATE w'min), and the
old integrand ds_dwp_f (eq 50), which contains cos^2(theta) explicitly and is now fed the exact
cos(theta).

Run: work/.venv/bin/python work/check_wp_production_spike.py [out.json]

Author: Eliana L. De Paoli.
Part of "Sensitivity to dark photons in a reactor antineutrino experiment with Skipper CCDs".
Written with agent assistance; provenance/ records what was derived, recalled and verified.
"""
import json
import sys

import numpy as np

OUT = {}
me = 0.51099895          # MeV
alfa = 1/137.035999084
hbarc = 1.973269804e-11  # cm MeV
WMAX = 10.0              # MeV, as in the macro


# --- the macro's functions -------------------------------------------------
def cos_exact(M, w, wp):                       # cos_theta_f, already corrected in DP3_Danilov.C
    return ((wp - me)*w + me*wp - M*M/2)/(w*np.sqrt(wp*wp - M*M))


def fsu(M, w, wp):                             # eq (46)
    s, u = 2*me*w, M*M - 2*me*wp
    return (2/(s + u))*(1 - (me*M)**2/(s*u))


def ds_dwp_50(M, w, wp, eps=1.0):              # eq (50) = ds_dwp_f
    ket = wp/w + w/wp + M*M*fsu(M, w, wp) + cos_exact(M, w, wp)**2 - 1
    return np.pi*alfa**2/me*(1 + M*M/(2*me*w))*(eps/w)**2*ket*hbarc**2


def F_new(M, w, wp):
    a, b, c, d = 2*me*w, M*M - 2*me*wp, me*me, M*M
    S = 1/a + 1/b
    return -2*(a/b + b/a) + 4*(2*c + d)*(S + c*S*S - d/(a*b))


def ds_dwp_new(M, w, wp, eps=1.0):             # ds_dwp_new_f: no angle anywhere
    return np.pi*alfa**2*eps**2*F_new(M, w, wp)/(2*me*w*w)*hbarc**2


def xpm(M, w):
    s = me*me + 2*me*w
    lam = max((s - (me + M)**2)*(s - (me - M)**2), 0.0)
    base, r = (w + me)*(s + M*M - me*me), w*np.sqrt(lam)
    return (base - r)/(2*s), (base + r)/(2*s)


def invwpmin(M, wp):                           # invwpmin_f: inverse of the APPROXIMATE w'min
    return np.inf if wp >= me/2 else (me*wp - M*M/2)/(me - 2*wp)


def Rw_old(M, wp, wmax=WMAX):                  # Rw_f: [w', min(invwpmin, wmax)], empty if that is <= w'
    hi = min(invwpmin(M, wp), wmax)
    return None if (wp <= M/2 or hi <= wp) else (wp, hi)


def Rw_new(M, wp, wmax=WMAX, n=200000):        # Rw_new_f: {w : x_-(w) <= w' <= x_+(w)}
    Eth = M + M*M/(2*me)
    ws = np.linspace(Eth*(1 + 1e-12), wmax, n)
    lo = np.array([xpm(M, w)[0] for w in ws])
    hi = np.array([xpm(M, w)[1] for w in ws])
    sel = (lo <= wp) & (wp <= hi)
    return None if not sel.any() else (ws[sel].min(), ws[sel].max())


def section(t):
    print(f"\n{t}\n" + "-"*len(t))


# ---------------------------------------------------------------------------
section("1. Eq (50) with the exact cos(theta) diverges as w' -> mA'")
# cos(theta) ~ 1/|k'| = 1/sqrt(w'^2-M^2), so cos^2 ~ 1/(w'-M): eq (50) has a 1/(w'-M) pole at
# the very point the new grid now starts from. The exact dsigma/dw' has no angle in it at all.
for M in (0.1, 0.5, 1.0):
    w = 3.0
    eps_rel = np.array([1e-2, 1e-3, 1e-4, 1e-5, 1e-8, 1e-9])
    v50 = np.array([ds_dwp_50(M, w, M*(1 + e)) for e in eps_rel])
    vnew = np.array([ds_dwp_new(M, w, M*(1 + e)) for e in eps_rel])
    # slope from the last two points: the pole only dominates the bracket once cos^2 >> w'/w + w/w'
    slope = np.log(v50[-1]/v50[-2])/np.log(eps_rel[-1]/eps_rel[-2])
    slope_new = np.log(vnew[-1]/vnew[-2])/np.log(eps_rel[-1]/eps_rel[-2])
    flat = vnew.max()/vnew.min()
    print(f"   mA'={M}: eq(50) grows x{v50[3]/v50[0]:.0f} between w'=1.01mA' and 1.0001mA', and as"
          f" (w'-mA')^{slope:+.3f} asymptotically; la exacta tiende a una constante (exponente {slope_new:+.1e})")
    assert slope < -0.9, (M, slope)          # eq (50) has a pole at w' = mA'
    assert abs(slope_new) < 1e-3, (M, slope_new)   # the exact one has none: it tends to a constant
    # (for mA'=1 MeV the exact value still varies by ~2 over that span, but from the u-channel pole
    #  at w' = mA'^2/2me = 0.978 MeV nearby, not from w' -> mA'; that region is outside [x_-,x_+])
    OUT[f"ds50_slope_wp_to_mAp_M{M}"] = float(slope)
    OUT[f"ds50_growth_M{M}"] = float(v50[3]/v50[0])
print("   PASS: eq (50) ~ (w'-mA')^-1 at fixed w=3 MeV; the exact dsigma/dw' is flat there")

# ---------------------------------------------------------------------------
section("2. The old domain Rw_f does not protect that region; the exact one does")
for M in (0.1, 0.5, 1.0):
    wp = M*(1 + 1e-4)          # ~ the first point of wp_production_new_f (delta_frac = 1e-4)
    old, new = Rw_old(M, wp), Rw_new(M, wp)
    print(f"   mA'={M}, w'={wp:.6f}:")
    print(f"      Rw_f (viejo):  {'vacio' if old is None else f'w in [{old[0]:.3f}, {old[1]:.3f}] MeV'}")
    print(f"      Rw_new_f:      {'VACIO -> cero' if new is None else f'w in [{new[0]:.3f}, {new[1]:.3f}] MeV'}")
    if new is not None:
        ws = np.linspace(*new, 2001)
        cmax = max(abs(cos_exact(M, w, wp)) for w in ws)
        assert cmax <= 1 + 1e-9, (M, cmax)
        OUT[f"max_abs_cos_in_exact_domain_M{M}"] = float(cmax)
        frac = (new[1] - new[0])/(old[1] - old[0])
        OUT[f"domain_width_ratio_M{M}"] = float(frac)
        print(f"      |cos| <= {cmax:.4f} en el dominio exacto; su ancho es {frac:.3f} del viejo")
    else:
        OUT[f"domain_width_ratio_M{M}"] = 0.0
assert Rw_new(1.0, 1.0001) is None and Rw_old(1.0, 1.0001) is not None
print("   PASS: |cos(theta)| <= 1 inside the exact domain (no pole can be reached from there);")
print("         for mA'=1 MeV the exact domain is empty up to min_w x_-(w) = 1.2547 MeV, the old one is not")

# ---------------------------------------------------------------------------
section("3. Where the old grid used to start, and why nothing showed before")
for M in (0.1, 0.5, 1.0):
    old_start = M/2 + 1e-4*(WMAX - M/2)        # wp_production_OLD_f
    new_start = M + 1e-4*(xpm(M, WMAX)[1] - M)  # wp_production_new_f
    guarded = Rw_old(M, old_start) is None
    nan = old_start < M
    OUT[f"old_grid_start_M{M}"] = float(old_start)
    OUT[f"new_grid_start_M{M}"] = float(new_start)
    print(f"   mA'={M}: grilla vieja arrancaba en w'={old_start:.4f} MeV (< mA': {nan}) -> "
          f"{'Rw_f devolvia ceros' if guarded else 'cos_theta_f = NaN (raiz de negativo)'}")
    print(f"            grilla nueva arranca en w'={new_start:.6f} MeV, justo sobre el polo de eq (50)")
    assert nan                                  # every old start is below the mass
print("   PASS: every point of the old grid below mA' was either NaN or killed by Rw_f's guard,")
print("         so the 1/(w'-mA') pole of eq (50) was never sampled until the grid moved up to mA'")

# ---------------------------------------------------------------------------
section("4. What each combination gives at the first grid point (w = 3 MeV, eps = 1)")
for M in (0.1, 0.5, 1.0):
    wp = M*(1 + 1e-4)
    far = ds_dwp_new(M, 3.0, min(2.0, 0.9*xpm(M, 3.0)[1]))
    print(f"   mA'={M}: eq(50) = {ds_dwp_50(M,3.0,wp):.3e} cm2/MeV   exacta = {ds_dwp_new(M,3.0,wp):.3e}"
          f"   (exacta lejos del borde, w'~2 MeV: {far:.3e})")

section("5. The height of the spike is set by delta_frac, which is the signature of a pole")
# A pole means the first grid point's value is not physics but a choice: it scales like 1/delta.
for M in (0.1, 0.5, 1.0):
    hs = []
    for dfrac in (1e-3, 1e-4, 1e-5, 1e-7, 1e-8):
        wp = M + dfrac*(xpm(M, WMAX)[1] - M)
        hs.append(ds_dwp_50(M, 3.0, wp))
    r1, r2, rasym = hs[1]/hs[0], hs[2]/hs[1], hs[4]/hs[3]
    OUT[f"spike_ratio_per_decade_M{M}"] = float(r2)
    OUT[f"spike_ratio_per_decade_asymptotic_M{M}"] = float(rasym)
    print(f"   mA'={M}: eq(50) en el primer punto = {hs[0]:.3e}, {hs[1]:.3e}, {hs[2]:.3e} cm2/MeV"
          f" para delta_frac = 1e-3, 1e-4, 1e-5  (x{r1:.1f}, x{r2:.1f} por decada; x{rasym:.1f} ya en el regimen del polo, 1e-7 -> 1e-8)")
    assert rasym > 9, (M, rasym)
    assert hs[2] > hs[1] > hs[0], (M, hs)
print("   PASS: dividing delta_frac by 10 multiplies the first point by ~10. Physics does not do that;")
print("         it is the 1/(w'-mA') pole of eq (50) being sampled closer in.")

print("\nALL ASSERTIONS PASSED")
if len(sys.argv) > 1:
    json.dump(OUT, open(sys.argv[1], "w"), indent=1)
    print(f"wrote {sys.argv[1]}")
