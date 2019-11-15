---
# 2d Cell with Out-of-Plane Wavevector
---

A 2d cell with a `k_point` that has a *non-zero* $z$ component $\beta$ (e.g., planewaves incident on a 2d surface from an out-of-plane direction, propagating modes of fiber waveguides with 2d cladding cross section, etc.) would normally result in a 3d simulation with complex fields. However, Meep can model the $e^{i \beta z}$ dependence using a modified 2d cell with real fields via `kz_2d="real/imag"`  which improves performance with practically no loss in accuracy.

Mathematically, the $e^{i \beta z}$ term of the fields adds an $i\beta\hat{z}$ $\times$ cross-product to the curls of Maxwell's equations, which couples the in- (TE) and out-of-plane (TM) polarizations. To avoid complex fields (which doubles the floating-point storage requirements), in the case of real fields Meep implicitly stores $i$\*(TM fields) rather than (TM fields), in which case the $i$'s cancel in the update equations for time-stepping. This simple reformulation is equivalent to looking at the superposition of the fields at $\beta$ and the time-reversed fields at $-\beta$. Since most calculations of the fields such as flux, energy, etc are insensitive to this implicit $i$ factor, the fact that the TE fields are real but the TM fields are "imaginary" will not be noticeable to users but nevertheless important to keep in mind.

By default, however, `kz_2d="complex"` so that all fields are complex in which case the $i\beta$ term is implemented directly.
