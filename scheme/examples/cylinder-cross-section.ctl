(define-param r 0.7) ;; radius of cylinder
(define-param h 2.3) ;; height of cylinder

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
(define dair (* 1.0 wvl-max))

(define boundary-layers (list (make pml (thickness dpml))))
(set! pml-layers boundary-layers)

(define sr (+ r dair dpml))
(define sz (+ dpml dair h dair dpml))
(define cell (make lattice (size sr 0 sz)))
(set! geometry-lattice cell)
(set! dimensions CYLINDRICAL)
(set-param! m -1)

;; (is-integrated? true) necessary for any planewave source extending into PML
(define circ-src (list (make source
                         (src (make gaussian-src (frequency frq-cen) (fwidth dfrq) (is-integrated? true)))
                         (center (* 0.5 sr) 0 (+ (* -0.5 sz) dpml))
                         (size sr 0 0)
                         (component Er))
                       (make source
                         (src (make gaussian-src (frequency frq-cen) (fwidth dfrq) (is-integrated? true)))
                         (center (* 0.5 sr) 0 (+ (* -0.5 sz) dpml))
                         (size sr 0 0)
                         (component Ep)
                         (amplitude 0-1i))))

(set! sources circ-src)

(define box-z1 (add-flux frq-cen dfrq nfrq
                         (make flux-region (center (* 0.5 r) 0 (* -0.5 h)) (size r 0 0))))
(define box-z2 (add-flux frq-cen dfrq nfrq
                         (make flux-region (center (* 0.5 r) 0 (* +0.5 h)) (size r 0 0))))
(define box-r (add-flux frq-cen dfrq nfrq
                         (make flux-region (center r 0 0) (size 0 0 h))))

(run-sources+ 10)

(display-fluxes box-z1)

(save-flux "box-z1-flux" box-z1)
(save-flux "box-z2-flux" box-z2)
(save-flux "box-r-flux" box-r)

(reset-meep)

(define ncyl 2.0)
(set! geometry (list
                (make block
                  (material (make medium (index ncyl)))
                  (size r 0 h)
                  (center (* 0.5 r) 0 0))))

(set! geometry-lattice cell)

(set! pml-layers boundary-layers)

(set! sources circ-src)

(define box-z1 (add-flux frq-cen dfrq nfrq
                         (make flux-region (center (* 0.5 r) 0 (* -0.5 h)) (size r 0 0))))
(define box-z2 (add-flux frq-cen dfrq nfrq
                         (make flux-region (center (* 0.5 r) 0 (* +0.5 h)) (size r 0 0))))
(define box-r (add-flux frq-cen dfrq nfrq
                         (make flux-region (center r 0 0) (size 0 0 h))))

(load-minus-flux "box-z1-flux" box-z1)
(load-minus-flux "box-z2-flux" box-z2)
(load-minus-flux "box-r-flux" box-r)

(run-sources+ 100)

(display-fluxes box-z1 box-z2 box-r)
