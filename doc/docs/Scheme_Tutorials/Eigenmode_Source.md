---
# Eigenmode Source
---

This example demonstrates using the [`eigenmode-source`](../Scheme_User_Interface.md#eigenmode-source) to couple exclusively into a single waveguide mode. The structure, shown in the schematic below, is a 2d dielectric waveguide with ε=12 and out-of-plane electric field E<sub>z</sub>. The dispersion relation ω(k) for modes with *even* mirror symmetry in the y-direction is computed using [MPB](https://mpb.readthedocs.io) and also shown in the schematic as blue lines. Using this waveguide configuration, we will investigate two cases: (1) frequency of 0.15 (normalized) which is single mode, and (2) frequency of 0.35 which is multi mode. We will use the eigenmode source feature to excite just a single mode with unidirectional propagation in each case (labeled **A** and **B** in the dispersion relation) and compare the results to using a constant-amplitude source. Results will also be presented for a rotated waveguide axis.

<center>
![](../images/eigenmode_source.png)
</center>

The simulation script is in [examples/oblique-source.py](https://github.com/NanoComp/meep/blob/master/scheme/examples/oblique-source.ctl).

For the single-mode case of `fsrc=0.15`, a constant-amplitude current source excites both the waveguide mode and radiating fields (which lie within the light cone) in both directions. This is shown in the main inset of the figure above. The `eigenmode-source`, which excites the labeled **A** by setting `kx=0.4` and `bnum=1`, excites only the right-going waveguide mode. There are four key parameters to `eigenmode-source`: `direction`, `eig-kpoint`, `eig-band`, and `eig-parity`. Note that `eigenmode-source` is a line segment centered at the origin extending the length of the entire cell. The parameter `rot-angle` specifies the rotation angle of the waveguide axis and is initially 0. This enables `eig-parity` to include `EVEN-Y` and the cell to include a mirror symmetry plane in the y direction.

For the multi-mode case of `fsrc=0.35`, a constant-amplitude current source excites a superposition of both waveguide modes in addition to the radiating field. The `eigenmode-source` can excites only a given right-going mode: **A** (`kx=0.4`, `bnum=2`) or **B** (`kx=1.2`, `bnum=1`).

```scm
(set-param! resolution 50) ; pixels/μm

(set! geometry-lattice (make lattice (size 14 14 no-size)))

(set! pml-layers (list (make pml (thickness 2))))

(define-param rot-angle 0) ; rotation angle (in degrees) of waveguide, CCW around z-axis
(set! rot-angle (deg->rad rot-angle))

(set! geometry (list (make block
                       (center 0 0 0)
                       (size infinity 1 infinity)
                       (e1 (rotate-vector3 (vector3 0 0 1) rot-angle (vector3 1 0 0)))
                       (e2 (rotate-vector3 (vector3 0 0 1) rot-angle (vector3 0 1 0)))
                       (material (make medium (epsilon 12))))))

(define-param fsrc 0.15) ; frequency of eigenmode or continuous-wave (CW) source
(define-param kx 0.4)    ; initial guess for wavevector in x-direction of eigenmode 
(define-param bnum 1)    ; band index of eigenmode

(define-param eig-src? true)

(set! sources (list
               (if eig-src?
                   (make eigenmode-source
                     (src (make continuous-src (frequency fsrc)))
                     (center 0 0 0)
                     (size 0 14 0)
                     (direction (if (= rot-angle 0) AUTOMATIC NO-DIRECTION))
                     (eig-kpoint (rotate-vector3 (vector3 0 0 1) rot-angle (vector3 kx 0 0)))
                     (eig-band bnum)
                     (eig-parity (if (= rot-angle 0) (+ EVEN-Y ODD-Z) ODD-Z))
                     (eig-match-freq? true))
                   (make source
                     (src (make continuous-src (frequency fsrc)))
                     (center 0 0 0)
                     (size 0 2 0)
                     (component Ez)))))

(if (= rot-angle 0)
    (set! symmetries (list (make mirror-sym (direction Y)))))

(run-until 100 (in-volume (volume (center 0 0 0) (size 10 10 0))
                          (at-beginning output-epsilon)
                          (at-end output-efield-z)))
```

There are numerical dispersion artifacts due to the FDTD spatial and temporal discretizations which create negligible backwards (leftward-propagating) wave artifacts by eigenmode current source, carrying approximately 10<sup>-5</sup> of the power of the desired rightward-propagating mode.

We can also demonstrate the eigenmode source for a rotated waveguide axis for both the single- and multi-mode case. The results are shown in the two figures below. There is one subtlety: for mode **A** in the multi-mode case, the `bnum` parameter be set to 3 rather than 2. This is because non-zero rotation angles break the symmetry in the y direction and therefore preclude the use of `EVEN-Y` in `eig-parity`. Without any parity specified for the y direction, the second band corresponds to odd modes which is why we must select the third band which contains even modes.

<center>
![](../images/oblique_source_singlemode.png)
</center>

<center>
![](../images/oblique_source_multimode.png)
</center>
