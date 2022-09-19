; Material dispersion example, from the Meep tutorial.  Here, we simply
; simulate homogenous space filled with a dispersive material, and compute
; its modes as a function of wavevector k.  Since omega/c = k/n, we can
; extract the dielectric function epsilon(omega) = (ck/omega)^2.

(set! geometry-lattice (make lattice (size no-size no-size no-size)))
(set-param! resolution 20)

; We'll use a dispersive material with two polarization terms, just for
; illustration.  The first one is a strong resonance at omega=1.1,
; which leads to a polaritonic gap in the dispersion relation.  The second
; one is a weak resonance at omega=0.5, whose main effect is to add a
; small absorption loss around that frequency.

(set! default-material
      (make dielectric (epsilon 2.25)
	    (polarizations
	     (make polarizability
	       (omega 1.1) (gamma 1e-5) (sigma 0.5))
	     (make polarizability
	       (omega 0.5) (gamma 0.1) (sigma 2e-5))
	     )))

(define-param fcen 1.0)
(define-param df 2.0)
(set! sources (list (make source
		      (src (make gaussian-src (frequency fcen) (fwidth df)))
		      (component Ez) (center 0 0 0))))

(define-param kmin 0.3)
(define-param kmax 2.2)
(define-param k-interp 99)
(define kpts (interpolate k-interp (list (vector3 kmin) (vector3 kmax))))

(define all-freqs (run-k-points 200 kpts)) ; a list of lists of frequencies

(map (lambda (kx fs)
       (map (lambda (f)
	      (print "eps:, " (real-part f) ", " (imag-part f)
		     ", " (sqr (/ kx f)) "\n"))
	    fs))
     (map vector3-x kpts) all-freqs)
