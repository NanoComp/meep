(set-param! resolution 50)   ; pixels/μm

(define-param sxy 4)
(define-param dpml 1)

(set! geometry-lattice (make lattice (size (+ sxy (* 2 dpml)) (+ sxy (* 2 dpml)) no-size)))

(set! pml-layers (list (make pml (thickness dpml))))

(define-param fcen 1.0)
(define-param df 0.4)
(define-param src-cmpt Ez)
(set! sources (list (make source
                      (src (make gaussian-src (frequency fcen) (fwidth df)))
                      (center 0)
                      (component src-cmpt))))

(if (= src-cmpt Ex)
    (set! symmetries (list (make mirror-sym (direction X) (phase -1))
                           (make mirror-sym (direction Y) (phase +1)))))
(if (= src-cmpt Ey)
    (set! symmetries (list (make mirror-sym (direction X) (phase +1))
                           (make mirror-sym (direction Y) (phase -1)))))
(if (= src-cmpt Ez)
    (set! symmetries (list (make mirror-sym (direction X) (phase +1))
                           (make mirror-sym (direction Y) (phase +1)))))

(define nearfield-box
  (add-near2far fcen 0 1
		(make near2far-region (center 0 (* +0.5 sxy)) (size sxy 0) (weight +1))
		(make near2far-region (center 0 (* -0.5 sxy)) (size sxy 0) (weight -1))
		(make near2far-region (center (* +0.5 sxy) 0) (size 0 sxy) (weight +1))
		(make near2far-region (center (* -0.5 sxy) 0) (size 0 sxy) (weight -1))))

(define flux-box
  (add-flux fcen 0 1
	    (make flux-region (center 0 (* +0.5 sxy)) (size sxy 0) (weight +1))
	    (make flux-region (center 0 (* -0.5 sxy)) (size sxy 0) (weight -1))
	    (make flux-region (center (* +0.5 sxy) 0) (size 0 sxy) (weight +1))
	    (make flux-region (center (* -0.5 sxy) 0) (size 0 sxy) (weight -1))))

(run-sources+ (stop-when-fields-decayed 50 src-cmpt (vector3 0 0) 1e-8))

(print "near-flux:, " (list-ref (get-fluxes flux-box) 0) "\n")

(define-param r (/ 1000 fcen)) ; half side length of far-field square box OR radius of far-field circle
(define-param res-ff 1)        ; resolution of far fields (points/μm)
(define far-flux (+ (list-ref (flux nearfield-box Y (volume (center 0 r 0) (size (* 2 r) 0 0)) res-ff) 0)
                    (- (list-ref (flux nearfield-box Y (volume (center 0 (- r) 0) (size (* 2 r) 0 0)) res-ff) 0))
                    (list-ref (flux nearfield-box X (volume (center r 0 0) (size 0 (* 2 r) 0)) res-ff) 0)
                    (- (list-ref (flux nearfield-box X (volume (center (- r) 0 0) (size 0 (* 2 r) 0)) res-ff) 0))))

(print "far-flux-box:, " far-flux "\n")

(define-param npts 100)  ; number of points in [0,2*pi) range of angles
(map (lambda (n)
       (let ((ff (get-farfield nearfield-box (vector3 (* r (cos (* 2 pi (/ n npts)))) (* r (sin (* 2 pi (/ n npts)))) 0))))
	 (print "farfield:, " n ", " (* 2 pi (/ n npts)))
	 (map (lambda (m)
		(print ", " (list-ref ff m)))
	      (arith-sequence 0 1 6))
	 (print "\n")))
         (arith-sequence 0 1 npts))
