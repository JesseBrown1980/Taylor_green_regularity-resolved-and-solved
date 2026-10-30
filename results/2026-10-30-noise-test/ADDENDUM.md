# Addendum 1 (sealed before the Rust run; the C++ harness had already run once)

Operator rule: the test runs in Rust 1.81 with clippy -D warnings, integer and ternary only,
HBI -> HBP -> SHA -> SH -> HASH, json=0, no JS, no Node, no V8.

- The first run used a C++ harness (`harness.cpp`, output kept as `RESULT.cpp.hbp`). It stays on
  record, and it is not the registered implementation.
- The registered implementation is now `tgrust/` (Rust 1.81, no dependencies, integer only). It
  rebuilds the repo's 288 gates itself, steps the dynamics, and computes every diagnostic.
- The repo's initialiser computes the Taylor–Green start state with floating-point sin/cos. That
  step belongs to the subject, not to the test. `dump_init.cpp` runs it once and exports the start
  states as trit files (`INIT_TGV.trits`, `INIT_NULL.trits`, one byte per channel, values 0..2) and
  the repo's gate tables (`GATES.lut`). These files are pinned by SHA-256 before the Rust run.
- Added prediction R1: the Rust gate tables equal `GATES.lut` byte for byte, and the Rust output equals
  `RESULT.cpp.hbp` byte for byte.
- P1 to P5 and the decision rule are unchanged, and the Rust output decides them.
