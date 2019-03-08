(set-param! resolution 128)  ; pixels/μm

(define Si (make dielectric (index 3.45)))

(define-param syz 10)
(set! geometry-lattice (make lattice (size no-size syz syz)))

(define-param a 1.0) ; waveguide width
(define-param s 1.0) ; waveguide separation distance

(set! geometry (list
		(make block (center 0 (* -0.5 (+ s a)) 0)
		      (size infinity a a) (material Si))
		(make block (center 0 (* 0.5 (+ s a)) 0)
		      (size infinity a a) (material Si))))

(set! k-points (list (vector3 0.5 0 0)))

(set-param! num-bands 1)
(set-param! tolerance 1e-9)

(define-param yodd? true)
(if yodd? (run-yodd-zodd) (run-yeven-zodd))

(print "data:, " s ", " (list-ref freqs 0) ", "
       (list-ref (compute-group-velocity-component (vector3 1 0 0)) 0) "\n")
