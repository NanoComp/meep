;; From the Meep tutorial: plotting Faraday rotation of a linearly polarized plane wave

;; Parameters for a gyrotropic Lorentzian medium
(define-param epsn 1.5)    ; background permittivity
(define-param f0 1.0)      ; natural frequency
(define-param g0 1e-6)     ; damping rate
(define-param sn 0.1)      ; sigma parameter
(define-param b0 0.15)     ; magnitude of bias vector

(set! default-material
      (make dielectric
	(epsilon epsn)
	(E-susceptibilities
	 (make gyrotropic-lorentzian-susceptibility
	   (frequency f0)
	   (sigma sn)
	   (gamma g0)
	   (bias (vector3 0 0 b0))))))

;; Set up and run the Meep simulation:
(define-param tmax 100)
(define-param L 20.0)
(define-param fsrc 0.8)
(define-param src-z -8.5)
(set-param! resolution 50)

(set! geometry-lattice (make lattice (size 0 0 L)))

(set! pml-layers (list (make pml (thickness 1.0) (direction Z))))

(set! sources (list
	       (make source
		 (src (make continuous-src (frequency fsrc)))
		 (component Ex)
		 (center (vector3 0 0 src-z)))))

(run-until tmax
	   (to-appended "efields"
			(at-end output-efield-x)
			(at-end output-efield-y)))
