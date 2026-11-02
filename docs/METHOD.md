# Method

## System under study

A bit-flip repetition-code memory of distance *d* ∈ {3, 5} data qubits with *d* − 1 check ancillas and *d* syndrome-extraction rounds. One round has two layers of controlled-NOT (CNOT) gates (data to ancilla), an idle period on the data qubits during readout, and an ancilla readout. After the last round all data qubits are read out. Initial state and ancilla reset are ideal.

Three independent noise sources, each a bit flip (X fault) with probability set per location:

| Source | Locations | Probability |
|---|---|---|
| gate | both outputs of every CNOT | *p*<sub>g</sub> |
| readout | every ancilla measurement and the final data readout | *p*<sub>m</sub> |
| idle | every data qubit once per round, during readout | *p*<sub>i</sub> |

The default mix is *p*<sub>g</sub> = *s*, *p*<sub>m</sub> = 2*s*, *p*<sub>i</sub> = *s* with noise scale *s*. All faults are single-qubit X, so fault effects are exactly linear: the fault table (`circuit.py`) propagates each fault through the CNOT schedule and records which detectors and which logical observable it flips. Merging identical effects with the exclusive-or rule gives the detector error model (DEM), which is sampled in NumPy (`dem.py`). Every mechanism flips at most two detectors (asserted in the tests).

## Decoder

Minimum-weight matching on the DEM graph (`decoder.py`): shortest paths between detectors and a boundary node via Dijkstra, then an exact subset dynamic programme (up to 14 detection events; networkx matching beyond). Edge weights are ln((1 − *p*)/*p*) of the model the decoder assumes.

## Budget predictors

Let ε<sub>s</sub> be the logical failure probability with only source *s* on, ε the value with all on, and *t* = (*d* + 1)/2 the minimum number of faults that causes a logical error.

* **Additive (baseline, standard practice):** Σ<sub>s</sub> ε<sub>s</sub>.
* **Lumped classical majority vote (classical baseline):** each data bit flips with probability *q* = rounds·(*p*<sub>i</sub> + 2*p*<sub>g</sub>) + *p*<sub>m</sub>; the block fails if at least *t* of *d* bits flip. It ignores syndrome-measurement noise.
* **Power-mean (proposed):** (Σ<sub>s</sub> ε<sub>s</sub><sup>1/*t*</sup>)<sup>*t*</sup>. Rescaling source *s* by κ<sub>s</sub> gives (Σ<sub>s</sub> κ<sub>s</sub> ε<sub>s</sub><sup>1/*t*</sup>)<sup>*t*</sup>, versus Σ<sub>s</sub> κ<sub>s</sub><sup>*t*</sup> ε<sub>s</sub> for the additive budget.
* **Reference:** direct Monte Carlo of the combined noise with the matching decoder, the standard decoder baseline in the cited memory experiments.

## Metrics

Ratio prediction / Monte Carlo, and log10 of it (per seed); fraction of cells within ±25 %; root-mean-square (RMS) log10 error for the cost study.

## Experiments

1. **Accuracy.** Distance 3 at *s* ∈ {0.2, 0.4, 0.8, 1.6} % and distance 5 at *s* ∈ {1.2, 1.6, 2.0} % (single-source runs need enough events); 5 seeds; 100 000 (d=3) or 60 000 (d=5) shots per run.
2. **Validity map.** Calibrate once at *s* = 0.4 %, readout/gate ratio 2 (distance 3), extrapolate by rescaling over *s* ∈ {0.2 … 2.4} % × ratio ∈ {0.25 … 8}; compare with Monte Carlo (3 seeds, 100 000 shots per cell).
3. **Ablation.** (a) Decoder weights: matched to each run (single-source runs use their own model), one fixed full-noise model, or unit weights. (b) Exponent: *t* = 1 (additive), *t* = (*d* + 1)/2, and a pooled best-fit *t*.
4. **Robustness sweeps.** (a) Decoder mis-calibration: one fixed decoder, built from the full noise model with readout scaled, is used for every run; it assumes the readout probability κ ∈ {0.25 … 4} times the truth. (b) Calibration-shot sweep and Pareto front of cost versus accuracy against a 10⁶-shot reference. (c) Logical error against *s* up to 4 % for both distances, to show where the polynomial regime ends.

## Seeds and compute

Seeds 0–4 (3 for the map), NumPy `SeedSequence` streams keyed by experiment and parameters. The full run takes about 16 minutes (982 s) on a laptop central processing unit; no graphics processing unit is used.
