# Three Taylor–Green steps → 81 Clay/Clasy registers

This repo's original D3Q27 finite lattice-gas simulation was evaluated against the **actual** Clay/Clasy 81-mode Rust shell library, with a prediction-first finite-grid observable conversion. The complete implementation, CI and independently checked SHA-256 leaves live in the [Clay bridge directory](https://github.com/JesseBrown1980/continuum-Clay-variation/tree/main/bridge/taylor-green).

The two [pre-registered bridge predictions](PREDICTION.hbp) and [three-leaf extension](PREDICTION-3-LEAVES.hbp) were filed **before** the executable was run. The external comparator reads the original source at Taylor–Green commit ee71e0eea7fd61bfe5a45ebf2f590311742161f7 and checks the three SHA-pinned initial files (two trit states plus 288 lookup gates).

## Exact finite-grid observable conversion

The 24³ lattice is partitioned into 27 disjoint 8³ spatial blocks. For each block, its three **signed integer momentum sums** become three Clay field elements modulo 1,000,999: 27 × 3 = 81. Per-block momentum magnitude ≤9,216, so centered residues recover those 81 signed observables without wrap. Their sum equals total Taylor–Green momentum summed over axes, modulo the field prime.

This is **not** a lossless translation of 373,248 original trit occupancies. For example, adding a particle in the stationary zero-velocity channel changes density without changing any of the 81 momentum observables.

## Observed three-step movement

| Step | TGV trits changed | 81 Clay registers changed | Registers disagreeing with Clay's own update | TGV WJ | TGV amplitude A |
|---|---:|---:|---:|---:|---:|
| 0 | — | — | — | 108,032 | 35,328 |
| 1 | 169,912 | 81 | **81** | 1,225,664 | 37,128 |
| 2 | 237,204 | 80 | **81** | 1,876,856 | 34,600 |
| 3 | 244,184 | 80 | **81** | 2,038,608 | 33,688 |

The SHA-pinned gates reproduce the existing Taylor–Green KJ, WJ, E2 and amplitude records for t=0..3. Global momentum, occupancy mass and total channel kinetic energy remain conserved; the NULL arm is stationary, with zero projected defects. The scrambled control also develops large WJ, so these results do not attribute local noise to coherent vortex growth.

The exact three-leaf SHA256 chain is:

- Step 1: ee5036f22e94990de1992e4649302d603aa5a2891222b9333a512959c95b40be
- Step 2: 67c1ceb79d45761b68835910e0a05c76564785b77b487bff3fe6d5a49901b11e
- Step 3: c638926523a7b45e44e1396320a4f3a640c17633716e65824351a8a543819c9f

**Verdict:** The signed **81-observable conversion** is exact and verified on real inputs, with a mathematically proved range bound. But its nonlinear evolution **does not commute** with the Clay/Clasy shell update: all 81 components fail the comparison on each tested TGV step. That is a registered negative result, not a failed verification gate. Other maps, physical field reconstructions and continuum correspondence remain OPEN.

The receipts and code are linked through the [Clay 3-step bridge](https://github.com/JesseBrown1980/continuum-Clay-variation/blob/main/bridge/taylor-green/README.md). The receipt verifier independently checks the 81 signed-to-field conversions, SHA256 chain, reproduction values and deliberately corrupted data controls. No native Taylor–Green simulation code was modified by this documentation.
