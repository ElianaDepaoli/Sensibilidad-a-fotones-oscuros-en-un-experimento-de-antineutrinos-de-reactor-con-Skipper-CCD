# Danilov, Demidov and Gorbunov (2019)

Correction from the v3 redo: Eq. (5) contains coefficient 2, not 4. Eq. (9)
neglects detector absorption away from resonance and must not be promoted to an
exact resonance formula. The paper adopts an anti-Compton efficiency estimate
equivalent to one sixth for TEXONO. Its plotted line is now extracted from the
local vector PDF by `code/darkphoton_v3/redo.cpp::DanilovPoints`; independent
predictions retain a normalization discrepancy. [[../exercises/darkphoton-v3]]

Local source: `papers/darkphoton_v3/danilov2019.pdf`. The v3 calculation uses the paper's plasma mass, in-medium mass splitting, absorption-suppressed production probability, and detector reconversion probability. The implementation keeps the epsilon dependence in both denominators and uses a configurable constant absorption length; it does not reproduce the paper's material profile or detector efficiencies.
