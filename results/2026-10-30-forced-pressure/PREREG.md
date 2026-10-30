# Forced Taylor–Green lattice gas ("solar wind"), cloud seat. Pre-registration (frozen before code)

Same subject, gates, streaming and pinned start states as the unforced test (seal a83e8f18…).
Rust 1.81, std only, integer and ternary (trits 0..2), json=0. 480 ticks per arm.

## Forcing (a push)
A push in direction +d (d = x or y) at a cell moves one particle from the face channel -d into the
face channel +d, only if n[-d] >= 1 and n[+d] <= 1. Mass and kinetic-energy sums are unchanged
exactly (both channels have |c|^2 = 1); momentum changes by 2 along d. Pushes happen after
collision and before streaming. A cell is selected at tick t when (7*cell + 13*t) mod 8 = 0.

## Arms
- U: unforced (as before).
- F1 vortex drive: selected cells push x by sign sx(x)cy(y)cz(z) and push y by sign -cx(x)sy(y)cz(z)
  (the Taylor–Green pattern, the standard way a vortex is kept going in turbulence studies).
- F2 solar wind: a uniform +x stream with gusts. Selected cells push +x; during gusts
  ((t / 20) mod 3 = 0) a selected cell pushes twice.
- F3 sheared wind: as F2, but +x for y < 12 and -x for y >= 12 (a shear layer, as where the solar
  wind meets a magnetosphere).
All arms start from the Taylor–Green state.

## Integer diagnostics per tick
E2 (energy), N (mass), occupancy check, PX and PY (net momentum), A (coherent Taylor–Green
amplitude, as before), WJ (|curl j|^2 sum, as before) and CK = sum over the 216 blocks of 4x4x4 cells
of |block momentum|^2 (large-scale flow; cell noise averages down in it).

## Predictions
- Q1 every arm, every tick: occupancy in 0..2, E2 and N equal their tick-0 values.
- Q2 U: A falls to at most A(0)/10 by tick 100; PX = PY = 0 every tick.
- Q3 F1: SUSTAINED if A(t) > A(0)/10 for every tick 100..480, else FADED. Expected SUSTAINED.
- Q4 F2: PX is non-decreasing every tick (collisions and streaming conserve momentum; the wind only
  adds). A falls to at most A(0)/10 by tick 100 (a uniform wind drives a drift, not this vortex).
- Q5 F3: PX stays within one tick's pushes of 0 (the two halves push oppositely); CK averaged over
  ticks 240..480 exceeds U's CK average by more than 2x (the shear builds large-scale flow).
- Q6 measured, not predicted: steady-state A for F1, CK averages for every arm, and WJ for every
  arm against U's noise floor.
Any build failure or missing row is INCONCLUSIVE.
