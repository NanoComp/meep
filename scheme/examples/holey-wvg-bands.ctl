; Meep Tutorial: Hz-polarized transmission and reflection through a cavity
; formed by a periodic sequence of holes in a dielectric waveguide,
; with a defect formed by a larger spacing between one pair of holes.

; This structure is based on one analyzed in:
;    S. Fan, J. N. Winn, A. Devenyi, J. C. Chen, R. D. Meade, and
;    J. D. Joannopoulos, "Guided and defect modes in periodic dielectric
;    waveguides," J. Opt. Soc. Am. B, 12 (7), 1267-1272 (1995).

; Some parameters to describe the geometry:
(define-param eps 13) ; dielectric constant of waveguide
(define-param w 1.2) ; width of waveguide
(define-param r 0.36) ; radius of holes

; The cell dimensions
(define-param sy 12) ; size of cell in y direction (perpendicular to wvg.)
(define-param dpml 1) ; PML thickness (y direction only!)

(set! geometry-lattice (make lattice (size 1 sy no-size)))

(set! geometry
       (list (make block (center 0 0) (size infinity w infinity)
		   (material (make dielectric (epsilon eps))))
	      (make cylinder (center 0 0) (radius r) (height infinity) (material air))))

(set! pml-layers (list (make pml (direction Y) (thickness dpml))))
(set-param! resolution 20)

(define-param fcen 0.25) ; pulse center frequency
(define-param df 1.5) ; pulse freq. width: large df = short impulse

(set! sources (list
	       (make source
		 (src (make gaussian-src (frequency fcen) (fwidth df)))
		 (component Hz) (center 0.1234 0))))

(set! symmetries (list (make mirror-sym (direction Y) (phase -1))))

(define-param kx false) ; if true, do run at specified kx and get fields
(define-param k-interp 19) ; # k-points to interpolate, otherwise

(if kx
    (begin
      (set! k-point (vector3 kx))
      (run-sources+
       300 (at-beginning output-epsilon)
       (after-sources (harminv Hz (vector3 0.1234 0) fcen df)))
      (run-until (/ 1 fcen) (at-every (/ 1 fcen 20) output-hfield-z)))
    (run-k-points 300 (interpolate k-interp (list (vector3 0) (vector3 0.5)))))
