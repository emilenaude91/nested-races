# Replication code for "Nested Races" (v2, September 2026)

- `fast.py` — model primitives, exact laboratory best responses, iterated-best-response equilibrium solver for all four regulatory profiles, interstate continuation values, taxonomy, and thresholds.
- `run1.py` — baseline thresholds and the 400-draw random search (no Stag Hunt under gamma = 1).
- `run2.py` — taxonomy sweeps for gamma in {0.5, 1, 1.5, 2, 3}.
- `run3.py` — numerical verification of the cross-partial and symmetric FOC formulas, n-laboratory comparative statics, implementation thresholds, risk-dominance band, and the numbers in Tables 1–2.
- `run4.py` — thresholds as functions of monitoring precision (Table 3) and both figures.
- `run5.py` — feasible implementation example (rho = 16) and the precipice threshold.

Requires numpy, scipy, matplotlib. Each script runs in under a minute.
- `robust.py`, `rob1.py`–`rob3.py` — robustness: Tullock exponent r, prize-survives variant, unimodality checks for gamma > 1, falling-behind behavioural extension (eta), and sensitivity of thresholds to epsilon, d, rho.
- `wmap.py`, `llcheck.py` — strategic-form map over (W/ℓ, q0) for average, weakest-link, and leader-link risk (Figure 3), and the decisive-contest check for leader-link.
- `lead.py`–`lead5.py`, `leadfig.py` — head-start model: difference-form contest with lead δ and noise ς, leader-link risk; static asymmetric game (`lead.py`), period-2 value functions on a lead grid (`lead2.py`, `lead4.py`), two-period lead-as-stock game with regime commitment (`lead3.py`), robustness and who-paces map (`lead5.py`), Figure 4 (`leadfig.py`).
- `rfront.py`, `rfront2.py` — continuous regulatory intensity: sustainable symmetric partial-restraint levels (the "ladder of agreements") under average and weakest-link risk.

## v4 (13 Sep 2026, after consolidated review)
- `core4.py` — unified solver: average risk, genuine soft-min weakest-link risk (tau = 0.05), increasing-returns variant, and leader-link risk with winner-causes-catastrophe accounting (prize π_j(1−q_j)). All four regulatory profiles solved jointly; single-peakedness recorded; identity check at every point.
- `v4a.py` (taxonomy sweeps and thresholds), `v4b.py` (risk-dominance thresholds, decisiveness, belief map → Fig. 3), `v4c.py`/`v4d.py` (divisible restraint with all deviations and gains, Appendix D), `lead4c.py` (two-period lead model with corrected accounting and inertia selection, Appendix C), `figs4.py` (Figs. 2–4).
- Earlier scripts are retained for the record; v4 results supersede them where they differ (γ-based "weakest-link" results are now the increasing-returns remark).
- `lead5c.py` — period-2 regime re-choice (states re-select the regime given the realised lead; Pareto-best pure NE, else mutual acceleration); tests time consistency of mutual regulation.
- `advocacy.py` — laboratories' payoffs under each regulatory regime, by lead and by internalisation (Section 8).
- `asym.py` — Appendix E: asymmetric enforcement (m_U = 1, m_C varied) under all three aggregation rules; classifies each state's disposition and lists pure equilibria.
- `edge.py` — asymmetric capability (U's laboratory has a persistent speed multiplier); computed cases referenced in Section 7.5 but not reported in the paper.
