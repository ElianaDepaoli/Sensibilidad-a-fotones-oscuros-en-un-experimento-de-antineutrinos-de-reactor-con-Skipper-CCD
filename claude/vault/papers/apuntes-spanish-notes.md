# The Spanish note in apuntes/ (author's name omitted)
`apuntes/Dark_Photon_cross_section-1.pdf`. (A second, anonymous file `Dark_Photon_cross_Section.pdf` with the same skeleton exists in the folder; per the user's instruction of 2026-09-21 it is not reviewed.)

**Claims.** Tree-level γe→A'e: |M|² eq 40, dσ/dΩ eq 45, dσ/dω' eq 50, range eqs 24–25, M=0 total eq 52.
**Verdict** (`work/check_cross_section.py`): |M|² wrong (−39%…+33%; drops KN's m_e terms); |q|≈ω' approximation; threshold wrong (E=M); range wrong; dσ/dω' off up to 2× and has a pole at ω'=M²/2m inside its own range at M=1 MeV; eq 52 correct (= KN).
**Fixes** (Section 6 of `reports/report_cross_section_en.pdf`): cosθ from the exact relation with |q|=√(ω'²−M²); |M|² → F of [[dark-compton-cross-section]]; dσ/dω' from dσ/dt with dt=2m dω' (no angular Jacobian); range → x±; threshold → M+M²/2m. See [[cross-section-audit-2026-09-21]].

**Which error matters for N_A'** (`work/check_eq50_fold.py`): folded with Park's spectrum, the note's eq 50 overshoots the 3–8 MeV rate by 1.05/1.46/1.99 (M=0.1/0.5/1 MeV), all from the amplitude+Jacobian; the threshold/range error contributes nothing there (Table 1 is not used in eq 50; the range eqs 24–25 carry the same error) and matters only below E_A' ≈ 1.5 MeV, where the note also integrates across its pole.
