(set-param! resolution 100) ;; pixels/um

(define-param dpml 1.0)
(set! pml-layers (list (make pml (thickness dpml))))

(define-param r 1.0)    ;; radius of cylinder
(define-param dair 2.0) ;; air padding thickness

(define s (* 2 (+ dpml dair r)))
(set! geometry-lattice (make lattice (size s s no-size)))

(define-param wvl 1.0)
(define fcen (/ wvl))

;; (is-integrated? true) necessary for any planewave source extending into PML
(set! sources (list (make source
                      (src (make gaussian-src (frequency fcen) (fwidth (* 0.1 fcen)) (is-integrated? true)))
                      (center (+ (* -0.5 s) dpml) 0)
                      (size 0 s)
                      (component Ez))))

(set! symmetries (list (make mirror-sym (direction Y))))

(set! geometry (list (make cylinder
                       (material SiO2)
                       (center 0 0)
                       (radius r)
                       (height infinity))))

(set! k-point (vector3 0 0 0))

(define dft-fields (add-dft-fields (list Dz Ez) fcen fcen 1 #:yee-grid true (volume (center 0 0 0) (size (* 2 r) (* 2 r)))))

(define flux-box (add-flux fcen 0 1
                           (make flux-region (center (- r) 0) (size 0 (* 2 r)) (weight +1))
                           (make flux-region (center (+ r) 0) (size 0 (* 2 r)) (weight -1))
                           (make flux-region (center 0 (+ r)) (size (* 2 r) 0) (weight -1))
                           (make flux-region (center 0 (- r)) (size (* 2 r) 0) (weight +1))))

(run-sources+ 100)

(output-dft dft-fields "dft-fields-cylinder")

(display-fluxes flux-box)
