# Research gap: how should a quantum error correction (QEC) error budget combine noise sources?

## State of the art (from the supplied literature)

Papers are cited by first author, year and digital object identifier (DOI).

Recent experiments report logical performance of small codes on superconducting hardware and then explain it by *attributing* the logical error to individual noise sources. Harper et al. (2026, doi:10.1038/s41467-026-71773-6) vary parameters in circuit simulations of a heavy-hex code and conclude that measurement noise dominates logical-gate performance, and that idling during readout is highly detrimental to a memory. Caune et al. (2026, doi:10.1038/s41467-026-73331-6) run an 8-qubit stability experiment with a field-programmable gate array (FPGA) decoder, observing logical error suppression as the number of decoding rounds grows. Landmark repetition and surface-code memories (Kelly et al., 2015, doi:10.1038/nature14270; Krinner et al., 2022, doi:10.1038/s41586-022-04566-8) use matching decoders and compare to numerical models of the device. Rispler et al. (2026, doi:10.1038/s41534-026-01271-6) show that the *threshold* itself depends on how a decoder accounts for correlated error types under circuit-level noise. Kuo and Ouyang (2026, doi:10.1038/s41534-026-01212-3) study decoders for erasure, a different noise source with its own capacity limit.

## Related work

| Paper | Year | Role | What it does | What it does NOT cover |
|---|---|---|---|---|
| Harper et al., doi:10.1038/s41467-026-71773-6 | 2026 | recent | Identifies dominant noise sources of heavy-hex logic by varying simulation parameters | No rule for combining per-source contributions into a total, nor when such a rule fails |
| Caune et al., doi:10.1038/s41467-026-73331-6 | 2026 | recent | Low-latency real-time decoding, logical error versus rounds | Does not attribute the logical error to sources |
| Vezvaee et al., doi:10.1038/s41467-026-70011-3 | 2026 | recent | Dynamical decoupling with a [[4,2,2]] code on transmons | Error budget of the logical qubit not analysed |
| Rispler et al., doi:10.1038/s41534-026-01271-6 | 2026 | recent | Thresholds with Y-error-aware decoding under circuit noise | Sub-threshold budgeting not addressed |
| Kuo and Ouyang, doi:10.1038/s41534-026-01212-3 | 2026 | recent | Belief-propagation decoders for erasure | Bit-flip circuit noise with several sources not covered |
| Meyer et al., doi:10.1088/2058-9565/ae98b1 | 2026 | recent | Learns codes tailored to structured noise | Budget of a fixed code unaddressed |
| Shao et al., doi:10.1038/s42005-026-02681-w | 2026 | recent | Variational codes adapted to device noise models | Same |
| Kelly et al., doi:10.1038/nature14270 | 2015 | landmark | Repetition-code memory with repeated error detection; error suppressed with size | Per-source decomposition not tested |
| Krinner et al., doi:10.1038/s41586-022-04566-8 | 2022 | landmark | Distance-3 surface code with matching decoder | Budget compared with numerical model only qualitatively |
| Cong et al., doi:10.1038/s41567-019-0648-8 | 2019 | landmark | Quantum convolutional neural network (QCNN) that learns a code for an unknown error model | Assumes one noise channel, no source budget |
| Terhal et al., doi:10.1088/2058-9565/ab98a5 | 2020 | review | Scalable bosonic QEC and decoding with finite squeezing | Different code family |
| Harrington et al., doi:10.1038/s42254-022-00494-8 | 2022 | review | Dissipation engineering in QEC and measurement | No budgeting method |
| Verma and Rahwani, doi:10.1007/978-3-032-30963-1_32 | 2027 | adjacent | Classical simulation cost of quantum circuits across platforms | Does not target QEC budgets |

## The gap

The standard practice behind these attributions is the *one-at-a-time* budget: switch every source off except one, measure its logical error rate, and sum the rates (or read the largest as "dominant"). For a code of distance *d* a logical error needs at least *t* = (*d* + 1)/2 faults, so the combined logical error scales as a polynomial of degree *t* in the source rates, and its cross terms are missing from a sum of single-source rates. The literature supplied here reports source rankings but gives no tested rule for the *combination*, nor a map of the noise regimes where such a rule can be trusted, nor its sensitivity to decoder calibration and Monte Carlo cost. Open-source QEC toolkits expose per-source switch-off experiments, so a validated combination rule is directly usable there.

**Novelty (one sentence):** this study tests a calibration-light, parameter-free *power-mean* budget, ε ≈ (Σₛ εₛ^(1/t))^t, against the additive budget, a lumped classical majority-vote model and direct Monte Carlo, and maps where it is accurate for a circuit-level repetition-code memory with gate, readout and idle noise.

## Hypotheses and success metrics (falsifiable)

* **H1.** At the same operating point, the additive budget under-predicts the combined logical error by a factor that grows with *t*, while the power-mean budget is within ±25 % of Monte Carlo. *Falsified if* the median power-mean ratio (prediction / Monte Carlo) lies outside 0.75 to 1.25 at distance 3 and below-threshold scales, or if its median absolute log10 error is not smaller than that of the additive budget.
* **H2.** One calibration (single-source rates at one reference point) extrapolates by rescaling, ε(κ) = (Σₛ κₛ εₛ^(1/t))^t, across noise scale and readout-to-gate ratio. *Success metric:* fraction of grid cells within ±25 % of Monte Carlo, reported separately for scales well below threshold.
* **H3.** At equal Monte Carlo cost the power-mean budget is on the Pareto front of cost versus error (it beats direct simulation of the combined noise) in the rare-event regime.

## Feasibility and scale

Everything is classical simulation: bit-flip noise on a repetition code (distance 3 and 5, rounds = distance) is exactly a linear fault-propagation problem, so a detector error model is built without a stabiliser simulator and sampled in NumPy. No quantum advantage claim is made or relevant; the outcome concerns error attribution, which transfers to hardware where the same single-source switch-off measurements are available.
