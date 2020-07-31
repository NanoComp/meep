; Calculating 2d ring-resonator modes, from the Meep tutorial.

(define-param n 3.4) ; index of waveguide
(define-param w 1) ; width of waveguide
(define-param r 1) ; inner radius of ring

(define-param pad 4) ; padding between waveguide and edge of PML
(define-param dpml 2) ; thickness of PML

(define sxy (* 2 (+ r w pad dpml))) ; cell size
(set! geometry-lattice (make lattice (size sxy sxy no-size)))

; Create a ring waveguide by two overlapping cylinders - later objects
; take precedence over earlier objects, so we put the outer cylinder first.
; and the inner (air) cylinder second.
(set! geometry (list
		(make cylinder (center 0 0) (height infinity)
		      (radius (+ r w)) (material (make dielectric (index n))))
		(make cylinder (center 0 0) (height infinity)
		      (radius r) (material air))))

(set! pml-layers (list (make pml (thickness dpml))))
(set-param! resolution 10)

; If we don't want to excite a specific mode symmetry, we can just
; put a single point source at some arbitrary place, pointing in some
; arbitrary direction.  We will only look for Ez-polarized modes.

(define-param fcen 0.15) ; pulse center frequency
(define-param df 0.1)  ; pulse width (in frequency)
(set! sources (list
               (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (component Ez) (center (+ r 0.1) 0))))

; exploit the mirror symmetry in structure+source:
(set! symmetries (list (make mirror-sym (direction Y))))

(run-sources+ 300
	      (at-beginning output-epsilon)
	      (after-sources (harminv Ez (vector3 (+ r 0.1)) fcen df)))

; Output fields for one period at the end.  (If we output
; at a single time, we might accidentally catch the Ez field when it is
; almost zero and get a distorted view.)
(run-until (/ 1 fcen) (at-every (/ 1 fcen 20) output-efield-z))
