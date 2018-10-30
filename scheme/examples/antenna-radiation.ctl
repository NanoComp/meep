(set-param! resolution 50)
(define-param sxy 4)
(define-param dpml 1)
(set! geometry-lattice (make lattice (size (+ sxy (* 2 dpml)) (+ sxy (* 2 dpml)) no-size)))
(set! pml-layers (list (make pml (thickness dpml))))

(define-param fcen 1.0)
(define-param df 0.4)
(define-param src-cmpt Ez)
(set! sources (list (make source (src (make gaussian-src (frequency fcen) (fwidth df))) (center 0) (component src-cmpt))))

(if (= src-cmpt Ex)
    (set! symmetries (list (make mirror-sym (direction Y)))))
(if (= src-cmpt Ey)
    (set! symmetries (list (make mirror-sym (direction X)))))
(if (= src-cmpt Ez)
    (set! symmetries (list (make mirror-sym (direction X)) (make mirror-sym (direction Y)))))

(define nearfield
  (add-near2far fcen 0 1
		(make near2far-region (center 0 (* 0.5 sxy)) (size sxy 0))
		(make near2far-region (center 0 (* -0.5 sxy)) (size sxy 0) (weight -1))
		(make near2far-region (center (* 0.5 sxy) 0) (size 0 sxy))
		(make near2far-region (center (* -0.5 sxy) 0) (size 0 sxy) (weight -1))))

(run-sources+ (stop-when-fields-decayed 50 src-cmpt (vector3 0 0) 1e-8))

(define-param r (/* 1000 fcen))    ; 1000 wavelengths out from the source
(define-param npts 100)            ; number of points in [0,2*pi) range of angles
(map (lambda (n)
       (let ((ff (get-farfield nearfield (vector3 (* r (cos (* 2 pi (/ n npts)))) (* r (sin (* 2 pi (/ n npts)))) 0))))
	 (print "farfield:, " (number->string n) ", " (number->string (* 2 pi (/ n npts))))
	 (map (lambda (m)
		(print ", " (number->string (list-ref ff m))))
	      (arith-sequence 0 1 6))
	 (print "\n")))
         (arith-sequence 0 1 npts))
