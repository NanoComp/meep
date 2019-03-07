;; compute group velocity of a waveguide mode using two different methods
;; (1) ratio of Poynting flux to energy density
;; (2) via MPB from get-eigenmode-coefficients

(set-param! resolution 20)

(set! geometry-lattice (make lattice (size 10 5 no-size)))

(set! geometry (list (make block
                       (center 0 0 0)
                       (size infinity 1 infinity)
                       (material (make medium (epsilon 12))))))

(set! pml-layers (list (make pml (thickness 1))))

(define-param fsrc 0.15)

(set! sources (list (make eigenmode-source
                      (src (make gaussian-src (frequency fsrc) (fwidth (* 0.2 fsrc))))
                      (center -3 0 0)
                      (size 0 5 0)
                      (eig-band 1)
                      (eig-parity (+ ODD-Z EVEN-Y))
                      (eig-match-freq? true))))

(set! symmetries (list (make mirror-sym (direction Y))))

(define flux (add-flux fsrc 0 1 (make flux-region (center 3 0 0) (size 0 5 0))))
(define energy (add-energy fsrc 0 1 (make energy-region (center 3 0 0) (size 0 5 0))))
(run-sources+ 100)

(define res (get-eigenmode-coefficients flux (list 1) #:eig-parity (+ ODD-Z EVEN-Y)))
(define mode-vg (array-ref (list-ref res 1) 0 0))

(define poynting-flux (list-ref (get-fluxes flux) 0))
(define e-energy (list-ref (get-electric-energy energy) 0))
(define ratio-vg (/ (* 0.5 poynting-flux) e-energy))

(print "group-velocity:, " ratio-vg ", " mode-vg "\n")
