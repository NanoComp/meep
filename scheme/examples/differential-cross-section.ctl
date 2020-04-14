(define-param r 1.0) ;; radius of sphere

(define-param frq-cen 1.0)

(set-param! resolution 20) ;; pixels/um

(define dpml 0.5)
(define dair 1.5) ;; at least 0.5/frq_cen padding between source and near-field monitor

(define boundary-layers (list (make pml (thickness dpml))))
(set! pml-layers boundary-layers)

(define s (* 2 (+ dpml dair r)))
(define cell (make lattice (size s s s)))
(set! geometry-lattice cell)

;; circularly-polarized source with propagation axis along x
;; (is-integrated? true) necessary for any planewave source extending into PML
(define circ-pol-src (list
                      (make source
                       (src (make gaussian-src (frequency frq-cen) (fwidth (* 0.2 frq-cen)) (is-integrated? true)))
                       (center (+ (* -0.5 s) dpml) 0 0)
                       (size 0 s s)
                       (component Ez))
                      (make source
                       (src (make gaussian-src (frequency frq-cen) (fwidth (* 0.2 frq-cen)) (is-integrated? true)))
                       (center (+ (* -0.5 s) dpml) 0 0)
                       (size 0 s s)
                       (component Ey)
                       (amplitude 0+1i))))

(set! sources circ-pol-src)

(set! k-point (vector3 0))

(define box-flux (add-flux frq-cen 0 1
                  (make flux-region (center (- (* 2 r)) 0 0) (size 0 (* 4 r) (* 4 r)))))

(define nearfield-box (add-near2far frq-cen 0 1
                       (make near2far-region (center (- (* 2 r)) 0 0) (size 0 (* 4 r) (* 4 r)) (weight +1))
                       (make near2far-region (center (+ (* 2 r)) 0 0) (size 0 (* 4 r) (* 4 r)) (weight -1))
                       (make near2far-region (center 0 (- (* 2 r)) 0) (size (* 4 r) 0 (* 4 r)) (weight +1))
                       (make near2far-region (center 0 (+ (* 2 r)) 0) (size (* 4 r) 0 (* 4 r)) (weight -1))
                       (make near2far-region (center 0 0 (- (* 2 r))) (size (* 4 r) (* 4 r) 0) (weight +1))
                       (make near2far-region (center 0 0 (+ (* 2 r))) (size (* 4 r) (* 4 r) 0) (weight -1))))

(run-sources+ 10)

(display-fluxes box-flux)

(save-near2far "nearfield-box-n2f" nearfield-box)

(reset-meep)

(define nsphere 2.0)
(set! geometry (list
                (make sphere
                  (material (make medium (index nsphere)))
                  (radius r)
                  (center 0))))

(set! geometry-lattice cell)

(set! pml-layers boundary-layers)

(set! sources circ-pol-src)

(set! k-point (vector3 0))

(define nearfield-box (add-near2far frq-cen 0 1
                       (make near2far-region (center (- (* 2 r)) 0 0) (size 0 (* 4 r) (* 4 r)) (weight +1))
                       (make near2far-region (center (+ (* 2 r)) 0 0) (size 0 (* 4 r) (* 4 r)) (weight -1))
                       (make near2far-region (center 0 (- (* 2 r)) 0) (size (* 4 r) 0 (* 4 r)) (weight +1))
                       (make near2far-region (center 0 (+ (* 2 r)) 0) (size (* 4 r) 0 (* 4 r)) (weight -1))
                       (make near2far-region (center 0 0 (- (* 2 r))) (size (* 4 r) (* 4 r) 0) (weight +1))
                       (make near2far-region (center 0 0 (+ (* 2 r))) (size (* 4 r) (* 4 r) 0) (weight -1))))

(load-minus-near2far "nearfield-box-n2f" nearfield-box)

(run-sources+ 100)

(define-param npts 100)           ;; number of points in [0,pi) range of polar angles to sample far fields along semi-circle

(define-param ff-r (* 10000 r))

(map (lambda (n)
       (let ((ff (get-farfield nearfield-box (vector3* ff-r (vector3 (cos (* pi (/ n npts))) 0 (sin (* pi (/ n npts))))))))
        (print "farfield:, " n ", " (* pi (/ n npts)))
        (map (lambda (m)
              (print ", " (list-ref ff m)))
         (arith-sequence 0 1 6))
        (print "\n")))
 (arith-sequence 0 1 npts))
