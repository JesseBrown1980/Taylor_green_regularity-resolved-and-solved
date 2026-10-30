# Addendum 2 (sealed before the run it governs; the first run is kept as RESULT.hbp)

Run 1 measured B = 0 exactly at every tick of every arm while PV was large, so B as defined cannot
detect a pressure mode on this lattice. Q7 stays FAIL as registered (a measurement-design failure).

Run 2 adds, without changing anything in run 1:
- Bx = sum P3 * s2(x) and By = sum P3 * s2(y) separately, to locate the degeneracy.
- D = sum over cells rho * (s2(x) + s2(y)), the Taylor–Green pressure mode carried by density, using
  the isothermal relation P = c_s^2 rho from the same repository (c_s^2 = 2/3 there).
- Q9: for F1, mean D over ticks 240..480 is > 0 and more than 3x the absolute mean D of U.
Run 2 must reproduce every run-1 row and summary byte for byte; its new values are on extra rows.
