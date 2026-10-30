# Noise test: Taylor–Green enstrophy against a shuffled control (2026-10-30)

Pre-registered test of the repository's simulator, run unchanged. Rust 1.81, std only, integer and
ternary (trits 0..2), clippy `-D warnings` exit 0, HBI → HBP → SHA → SH → HASH, json=0.
Seal `TG.seal` = `a83e8f18…`.

## Arms (480 ticks each, the repository's gates and streaming)

- TGV: the repository's Taylor–Green start state.
- SCR: the same cell states shuffled across the grid (same density, energy and zero net momentum,
  no vortex).
- NULL: the uniform state.

## Results (`RESULT.hbp`)

| check | result |
|---|---|
| R0: the repository's binary reproduces its committed CSVs | byte-identical |
| R1: Rust gates equal the repository's 288 gates; Rust output equals the first C++ harness | byte-identical |
| P1: occupancy 0..2 and exact energy, every tick, every arm | PASS |
| P2: NULL stays at zero | PASS |
| P3: vortex present at the start, absent in SCR | PASS (A = 35,328 against 404) |
| P4: TGV enstrophy against SCR, ticks 5..480 | **NOISE**: within 10% on 476 of 476 ticks |

| tick | TGV enstrophy (∑\|curl j\|²) | SCR enstrophy | TGV vortex amplitude A |
|---|---|---|---|
| 0 | 108,032 | 424,624 | 35,328 |
| 1 | 1,225,664 | 1,470,888 | 37,128 |
| 15 | 1,988,256 | 2,021,212 | 17,408 |
| 48 | 1,974,704 | 1,990,136 | 3,360 |
| 480 | 2,056,952 | 1,998,084 | ≈ 0 |

- The occupancy bound holds at every tick: velocity and vorticity stay bounded, as the repository's
  boundedness theorem states.
- The enstrophy rise happens in one tick and reaches the same level, about 2.0 million, in the
  shuffled arm, which has no vortex. It is the lattice thermalising to its noise floor.
- The coherent vortex decays smoothly: A halves by tick 15 and reaches about 0 by tick 100.

## Files

`PREREG.md` (`ccc4db5d…`) was frozen before any code. `ADDENDUM.md` (`90322ba5…`) moved the
registered implementation to Rust 1.81 before the Rust run; the first C++ run is kept as
`RESULT.cpp.hbp`. The repository's float initialiser ran once in `dump_init.cpp` to export the
start states (`INIT_*.trits`) and gate tables (`GATES.lut`), pinned in `INIT.sha256`. `run.sh`
rebuilds `tgrust/` and checks byte identity with `RESULT.hbp`.
