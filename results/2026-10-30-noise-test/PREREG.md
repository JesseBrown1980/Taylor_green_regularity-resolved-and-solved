# Taylor–Green lattice-gas test, cloud seat. Pre-registration (frozen before any code)

Subject: JesseBrown1980/Taylor_green_regularity-resolved-and-solved at e7c5df1, simulator
`asolaria_taylor_green_regularity.cpp` (24^3 D3Q27, occupancy 0..2, 288 involutive gates, 32 per cell
per tick, gate index (t*32 + cell*17 + j*37) mod 288, then streaming). Its dynamics are reused
unchanged. Every diagnostic here is an exact integer; json=0.

## Arms (same simulator, same gates)
- TGV: the repo's initialiser, U0 = 0.20.
- SCR: the TGV initial state with whole cell states permuted between cells by a fixed integer
  Fisher–Yates shuffle (LCG x' = (1103515245 x + 12345) mod 2^31, seed 20261030). Same multiset of
  cell states, so same density (27 everywhere), same energy, zero net momentum, no spatial coherence.
- NULL: U0 = 0 (uniform state).
Run 480 ticks each.

## Integer diagnostics per tick
- rho, j = sum n_i c_i per cell (integers). E2 = sum n_i |c_i|^2 over the grid.
- KJ = sum over cells |j|^2.
- WJ = sum over cells |curl j|^2 with undivided central differences (2x the derivative), an integer
  proxy for enstrophy (rho = 27 at tick 0 in all arms).
- A = coherent Taylor–Green amplitude: sum of jx*sx(x)cy(y)cz(z) - jy*cx(x)sy(y)cz(z), where the
  sign patterns are the signs of sin/cos(2 pi (i + 1/2)/24): sx = +1 for i < 12 else -1;
  cy = +1 for i < 6 or i >= 18 else -1.

## Predictions
- R0 reproduction: the repo's own binary, compiled as documented, writes CSVs byte-identical to the
  committed taylor_green_regularity_results.csv and null_control_results.csv.
- P1 occupancy stays in 0..2 and E2 is constant, every tick, every arm.
- P2 NULL: KJ = WJ = A = 0 every tick.
- P3 A_TGV(0) > 0 and |A_SCR(0)| <= A_TGV(0)/10.
- P4 decision rule on WJ over ticks 5..480 (476 ticks):
    NOISE   if |WJ_TGV - WJ_SCR| * 10 <= WJ_SCR on at least 90% of those ticks (within 10%),
    VORTEX  if WJ_TGV * 10 > WJ_SCR * 11 on at least 90% of those ticks,
    MIXED   otherwise.
  Expected: NOISE. Reported whatever it is.
- P5 measured, not predicted: the first tick at which A_TGV <= A_TGV(0)/2, and the tick-0 to tick-1
  ratios of WJ in TGV and SCR.
Any build failure or missing row is INCONCLUSIVE.
