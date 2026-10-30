# Forced Taylor–Green ("solar wind") with pressure (2026-10-30)

Pre-registered follow-on to the noise test, on the same simulator, gates and pinned start state.
Rust 1.81, std only, integer and ternary, clippy `-D warnings` exit 0, HBI → HBP → SHA → SH → HASH,
json=0. Seal `TGW.seal` = `4085db07…`.

## Forcing

A push moves one particle from a face channel into the opposite face channel of the same speed, so
mass and energy stay exact and only momentum is injected. One cell in eight is pushed per tick.

- U: unforced.
- F1 vortex drive: pushes in the Taylor–Green pattern.
- F2 solar wind: a uniform +x stream with gusts (double pushes for 20 of every 60 ticks).
- F3 sheared wind: +x in one half of the box and −x in the other.

## Results (`RESULT.hbp`, ticks 240..480)

| arm | vortex amplitude A (start 35,328) | large-scale flow CK (noise ≈ 0.50 M) | net x-momentum at 480 | registered result |
|---|---|---|---|---|
| U | ≈ 0 | 0.50 M | 0 | fades (Q2 PASS) |
| F1 | 30,478 mean, 26,608 minimum | 3.1 M | 0 | **SUSTAINED** (Q3) |
| F2 | ≈ 0 | 165 M | +203,660, rising every tick | Q4 PASS |
| F3 | ≈ 0 | 8.0 M | 0 | Q5 PASS |

- Q1 PASS: occupancy 0..2, exact mass and exact energy at every tick of every arm.
- A drive in the vortex's own pattern holds the Taylor–Green vortex at 86% indefinitely.
- In the periodic box, momentum from a uniform wind has nowhere to go, so the whole box drifts. A
  through-flowing wind needs inflow and outflow boundaries.
- A sheared wind builds strong large-scale flow in its own pattern.
- Local enstrophy stays at the noise floor of about 2.0 million in every arm. The large-scale
  measures A and CK carry the flow signal.

## Pressure

Pressure is taken as in effective-viscosity-and-pressure-from-the-discrete-dynamics-in-trits:
P = 2E/3V, so 3P per cell = ∑ nᵢ\|cᵢ\|², an exact integer. Pressure fluctuations (PV) are about
1.15 million in U, F1 and F3 and half that under the uniform wind. The registered pressure-mode
measures, energy-weighted (Q7) and density-weighted (Q9, sealed in `ADDENDUM2.md` before its run),
read exactly 0 at every tick of every arm: the Taylor–Green symmetry kept by the lattice cancels a
weight of the shape cos 2kx + cos 2ky. Both stay FAIL as registered. The next pressure test is one
that compares vortex cores with stagnation points directly.

## Files

`PREREG.md` (`65f87a95…`), `ADDENDUM.md` (`e573d681…`, pressure) and `ADDENDUM2.md`
(`570bdd20…`, density mode) were each sealed before the run they govern. `RESULT.run1.hbp` is run 1;
`RESULT.hbp` is run 2, which reproduces every run-1 row byte for byte and adds the addendum-2 rows.
`run.sh` rebuilds `tgwind/` and checks both.
