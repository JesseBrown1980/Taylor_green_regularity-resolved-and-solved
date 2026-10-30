# Addendum 1: pressure (sealed before any code ran)

Pressure as defined in JesseBrown1980/effective-viscosity-and-pressure-from-the-discrete-dynamics-in-trits:
P = 2E/(3V) = Tr(Pi)/3. Per cell, as an exact integer: P3(x) = 3P(x) = sum_i n_i |c_i|^2.
Its grid mean is E2/V = 54 exactly, since E2 is conserved.

Added diagnostics per tick:
- PV = sum over cells (P3 - 54)^2, the pressure-fluctuation power.
- B = sum over cells P3 * (s2(x) + s2(y)), the Taylor–Green pressure mode. Navier–Stokes gives
  p proportional to +(cos 2kx + cos 2ky)(cos 2kz + 2) for this velocity field. Here
  s2(i) = +1 if (i mod 12) is in {0,1,2,9,10,11}, else -1 (the sign of cos 2k(i + 1/2)).

Added predictions:
- Q7 F1 (sustained vortex): mean B over ticks 240..480 is > 0 (the sign Navier–Stokes predicts) and
  more than 3x the absolute mean B of U over the same ticks.
- Q8 measured, not predicted: PV per arm against U, and B per arm.
