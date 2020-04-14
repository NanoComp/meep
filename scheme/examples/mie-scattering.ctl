(define-param r 1.0) ;; radius of sphere

(define wvl-min (/ (* 2 pi r) 10))
(define wvl-max (/ (* 2 pi r) 2))

(define frq-min (/ wvl-max))
(define frq-max (/ wvl-min))
(define frq-cen (* 0.5 (+ frq-min frq-max)))
(define dfrq (- frq-max frq-min))
(define nfrq 100)

;; at least 8 pixels per smallest wavelength, i.e. (floor (/ 8 wvl-min))
(set-param! resolution 25)

(define dpml (* 0.5 wvl-max))
(define dair (* 0.5 wvl-max))

(define boundary-layers (list (make pml (thickness dpml))))
(set! pml-layers boundary-layers)

(define symm (list (make mirror-sym (direction Y))
                   (make mirror-sym (direction Z) (phase -1))))
(set! symmetries symm)

(define s (* 2 (+ dpml dair r)))
(define cell (make lattice (size s s s)))
(set! geometry-lattice cell)

;; (is-integrated? true) necessary for any planewave source extending into PML
(define pw-src (make source
                 (src (make gaussian-src (frequency frq-cen) (fwidth dfrq) (is-integrated? true)))
                 (center (+ (* -0.5 s) dpml) 0 0)
                 (size 0 s s)
                 (component Ez)))
(set! sources (list pw-src))

(set! k-point (vector3 0))

(define box-x1 (add-flux frq-cen dfrq nfrq
                         (make flux-region (center (- r) 0 0) (size 0 (* 2 r) (* 2 r)))))
(define box-x2 (add-flux frq-cen dfrq nfrq
                         (make flux-region (center (+ r) 0 0) (size 0 (* 2 r) (* 2 r)))))
(define box-y1 (add-flux frq-cen dfrq nfrq
                         (make flux-region (center 0 (- r) 0) (size (* 2 r) 0 (* 2 r)))))
(define box-y2 (add-flux frq-cen dfrq nfrq
                         (make flux-region (center 0 (+ r) 0) (size (* 2 r) 0 (* 2 r)))))
(define box-z1 (add-flux frq-cen dfrq nfrq
                         (make flux-region (center 0 0 (- r)) (size (* 2 r) (* 2 r) 0))))
(define box-z2 (add-flux frq-cen dfrq nfrq
                         (make flux-region (center 0 0 (+ r)) (size (* 2 r) (* 2 r) 0))))

(run-sources+ 10)

(display-fluxes box-x1)

(save-flux "box-x1-flux" box-x1)
(save-flux "box-x2-flux" box-x2)
(save-flux "box-y1-flux" box-y1)
(save-flux "box-y2-flux" box-y2)
(save-flux "box-z1-flux" box-z1)
(save-flux "box-z2-flux" box-z2)

(reset-meep)

(define nsphere 2.0)
(set! geometry (list
                (make sphere
                  (material (make medium (index nsphere)))
                  (radius r)
                  (center 0))))

(set! geometry-lattice cell)

(set! pml-layers boundary-layers)

(set! symmetries symm)

(set! sources (list pw-src))

(set! k-point (vector3 0))

(define box-x1 (add-flux frq-cen dfrq nfrq
                         (make flux-region (center (- r) 0 0) (size 0 (* 2 r) (* 2 r)))))
(define box-x2 (add-flux frq-cen dfrq nfrq
                         (make flux-region (center (+ r) 0 0) (size 0 (* 2 r) (* 2 r)))))
(define box-y1 (add-flux frq-cen dfrq nfrq
                         (make flux-region (center 0 (- r) 0) (size (* 2 r) 0 (* 2 r)))))
(define box-y2 (add-flux frq-cen dfrq nfrq
                         (make flux-region (center 0 (+ r) 0) (size (* 2 r) 0 (* 2 r)))))
(define box-z1 (add-flux frq-cen dfrq nfrq
                         (make flux-region (center 0 0 (- r)) (size (* 2 r) (* 2 r) 0))))
(define box-z2 (add-flux frq-cen dfrq nfrq
                         (make flux-region (center 0 0 (+ r)) (size (* 2 r) (* 2 r) 0))))

(load-minus-flux "box-x1-flux" box-x1)
(load-minus-flux "box-x2-flux" box-x2)
(load-minus-flux "box-y1-flux" box-y1)
(load-minus-flux "box-y2-flux" box-y2)
(load-minus-flux "box-z1-flux" box-z1)
(load-minus-flux "box-z2-flux" box-z2)

(run-sources+ 100)

(display-fluxes box-x1 box-x2 box-y1 box-y2 box-z1 box-z2)
