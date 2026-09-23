# Author: Eliana L. De Paoli.
# Part of "Sensitivity to dark photons in a reactor antineutrino experiment with Skipper CCDs".
# Written with agent assistance; provenance/ records what was derived, recalled and verified.
# marimo notebook: figures of the exact production cross section.

import marimo

__generated_with = "0.19.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.integrate import quad
    import check_cross_section as cc

    mo.md("# Dark-Compton cross section: exact result vs the formulas in `apuntes/`")
    return cc, mo, np, plt, quad


@app.cell
def _(cc, np, plt):
    # Figure 1: differential cross section at E = 3 MeV, exact vs note 1
    def fig_dsigma_dx(E=3.0, masses=(0.1, 0.5, 1.0), out="fig_dsigma_dx_E3.pdf"):
        m, alpha = cc.m, cc.alpha
        HC2 = (1.973269804e-11)**2 * 1e27   # MeV^-2 -> 1e-27 cm^2
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharey=False)
        for ax, M in zip(axes, masses):
            xlo, xhi = cc.x_pm(E, m, M)
            xs = np.linspace(xlo + 1e-5, xhi - 1e-5, 400)
            ex = cc.dsig_dx_exact(E, xs, m, M, alpha)*HC2
            ax.plot(xs, ex, "k-", lw=2, label="exact tree level")
            d1lo, d1hi = (m*E + M*M/2)/(2*E + m), E
            xs1 = np.linspace(d1lo + 1e-5, d1hi - 1e-5, 400)
            ax.plot(xs1, cc.dsig_dx_apuntes(E, xs1, m, M, alpha)*HC2, "r--", lw=1.5, label="note 1, eq. (50)")
            ax.set_title(f"$E_\\gamma$ = {E} MeV, $m_{{A'}}$ = {M} MeV")
            ax.set_xlabel("$E_{A'}$ [MeV]")
            ax.set_xlim(0, E*1.02)
            ax.set_ylim(0, 1.6*ex.max())   # the notes' formula has a pole at E_A' = M^2/2m (u - m^2 = 0) inside their own range
            ax.axvspan(0, xlo, color="0.9"); ax.axvspan(xhi, E*1.02, color="0.9")
        axes[0].set_ylabel(r"$\epsilon^{-2}\,d\sigma/dE_{A'}$ [$10^{-27}$ cm$^2$/MeV]")
        axes[0].legend(fontsize=8, loc="upper left")
        fig.tight_layout()
        fig.savefig(out)
        return fig

    _f1 = fig_dsigma_dx()
    _f1
    return (fig_dsigma_dx,)


@app.cell
def _(cc, np, plt, quad):
    # Figure 2: total production cross section vs photon energy, exact vs Klein-Nishina and vs the notes
    def fig_sigma_total(masses=(0.1, 0.5, 1.0), out="fig_sigma_total.pdf"):
        m, alpha = cc.m, cc.alpha
        HC2 = (1.973269804e-11)**2 * 1e24   # MeV^-2 -> barn
        Es = np.linspace(0.3, 10, 120)
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(Es, [cc.sigma_KN(E, m, alpha)*HC2 for E in Es], "k-", lw=2, label="$m_{A'}=0$ (Klein-Nishina)")
        cols = ("C0", "C1", "C3")
        for M, c in zip(masses, cols):
            Eth = M + M*M/(2*m)
            Ev = Es[Es > Eth*1.001]
            ex = [quad(lambda x: cc.dsig_dx_exact(E, x, m, M, alpha), *cc.x_pm(E, m, M))[0]*HC2 for E in Ev]
            ax.plot(Ev, ex, "-", color=c, lw=1.8, label=f"exact, $m_{{A'}}$ = {M} MeV")
            # note 1 integrated over the EXACT range: over the note's own range the integrand has a pole at E_A' = M^2/2m
            d1 = [quad(lambda x: cc.dsig_dx_apuntes(E, x, m, M, alpha), *cc.x_pm(E, m, M))[0]*HC2 for E in Ev]
            ax.plot(Ev, d1, "--", color=c, lw=1.2, label=f"note 1 (exact range), $m_{{A'}}$ = {M} MeV")
        ax.set_yscale("log")
        ax.set_xlabel("$E_\\gamma$ [MeV]")
        ax.set_ylabel(r"$\epsilon^{-2}\,\sigma(\gamma e\to A'e)$ [barn / electron]")
        ax.set_ylim(1e-3, 1)
        ax.legend(fontsize=7.5, ncol=2)
        fig.tight_layout()
        fig.savefig(out)
        return fig

    _f2 = fig_sigma_total()
    _f2
    return (fig_sigma_total,)


@app.cell
def _(cc, np, plt, quad):
    # Figure 3: transverse / longitudinal content of the produced A' (Compton production)
    def fig_long_fraction(masses=(0.1, 0.3, 0.5, 1.0), out="fig_longitudinal_fraction.pdf"):
        m = cc.m
        fig, ax = plt.subplots(figsize=(6, 4))
        for M in masses:
            Eth = M + M*M/(2*m)
            Es = np.linspace(Eth*1.02, 10, 40)
            fr = []
            for E in Es:
                xlo, xhi = cc.x_pm(E, m, M)
                fT = quad(lambda x: cc.F_numeric(E, x, m, M, "T"), xlo, xhi)[0]
                fL = quad(lambda x: cc.F_numeric(E, x, m, M, "L"), xlo, xhi)[0]
                fr.append(fL/(fT + fL))
            ax.plot(Es, fr, lw=1.8, label=f"$m_{{A'}}$ = {M} MeV")
        ax.set_xlabel("$E_\\gamma$ [MeV]")
        ax.set_ylabel(r"$\sigma_L/\sigma_{\rm tot}$ (longitudinal fraction of produced $A'$)")
        ax.set_ylim(0, 0.35)
        ax.legend()
        fig.tight_layout()
        fig.savefig(out)
        return fig

    _f3 = fig_long_fraction()
    _f3
    return (fig_long_fraction,)


@app.cell
def _(cc, np, plt, quad):
    # Figure 4: dN_A'/dE_A' at the reactor for the two production models, Park spectrum, 1 GW, eps = 1
    def fig_spectra_two_models(masses=(0.1, 0.5, 1.0), out="fig_dN_dE_two_models.pdf"):
        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
        ax = axes[0]
        xs = np.concatenate([np.linspace(0.05, 0.6, 60), np.linspace(0.6, 8, 120)[1:]])
        cols = ("C0", "C1", "C3")
        for M, c in zip(masses, cols):
            yA = np.array([cc.dN_dx_compton(x, M) for x in xs])/1e21
            yB = np.array([cc.dN_dx_oscillation(x, M) for x in xs])/1e21
            ax.plot(xs, np.where(yA > 0, yA, np.nan), "-", color=c, lw=1.8, label=f"A: Compton-like, $m_{{A'}}$ = {M} MeV")
            ax.plot(xs, np.where(yB > 0, yB, np.nan), "--", color=c, lw=1.2, label=f"B: oscillation, $m_{{A'}}$ = {M} MeV")
        ax.set_yscale("log"); ax.set_ylim(1e-4, 3); ax.set_xlim(0, 8)
        ax.set_xlabel("$E_{A'}$ [MeV]"); ax.set_ylabel(r"$\epsilon^{-2}\,dN_{A'}/dE_{A'}$ [$10^{21}$ MeV$^{-1}$ s$^{-1}$]")
        ax.axvspan(3, 8, color="0.92", zorder=0); ax.text(5.5, 1.2, "TEXONO window\n3-8 MeV", fontsize=8, ha="center")
        ax.legend(fontsize=7, loc="lower left"); ax.set_title("A' spectrum at the reactor: 1 GW, Park spectrum, $\\epsilon=1$", fontsize=10)

        ax = axes[1]
        Ms = np.geomspace(0.05, 2.5, 30)
        NA = [sum(quad(lambda x: cc.dN_dx_compton(x, M), lo, hi, limit=100)[0] for lo, hi in ((3, 5), (5, 8))) for M in Ms]
        NB = [quad(cc.S_park, max(3, M), 8)[0] for M in Ms]
        ax.plot(Ms, np.array(NA)/1e18, "k-", lw=1.8, label="A: Compton-like (Park eq. 1, exact $\\sigma$)")
        ax.plot(Ms, np.array(NB)/1e18, "k--", lw=1.4, label="B: oscillation (Danilov, $P=\\epsilon^2$)")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(0.1, 50)
        ax.set_xlabel("$m_{A'}$ [MeV]"); ax.set_ylabel(r"$\epsilon^{-2}\,N_{A'}$ with $3\leq E_{A'}\leq 8$ MeV [$10^{18}$ s$^{-1}$]")
        ax.legend(fontsize=8, loc="lower left"); ax.set_title("A' per second in the TEXONO window", fontsize=10)
        fig.tight_layout(); fig.savefig(out)
        return fig

    _f4 = fig_spectra_two_models()
    _f4
    return (fig_spectra_two_models,)

if __name__ == "__main__":
    app.run()
