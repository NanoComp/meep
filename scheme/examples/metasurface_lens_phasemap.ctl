;; phase map of a binary-phase grating unit cell from the Meep tutorial

(set-param! resolution 50)  ; pixels/μm

(define-param dpml 1.0)     ; PML thickness
(define-param dsub 2.0)     ; substrate thickness
(define-param dpad 2.0)     ; padding between grating and PML
(define-param gp 0.3)       ; grating period
(define-param gh 1.8)       ; grating height
(define-param gdc 0.5)      ; grating duty cycle

(define sx (+ dpml dsub gh dpad dpml))
(define sy gp)

(define cell (make lattice (size sx sy no-size)))
(set! geometry-lattice cell)

(define boundary-layers (list (make pml (thickness dpml) (direction X))))
(set! pml-layers boundary-layers)

(define-param lcen 0.5)      ; center wavelength
(define fcen (/ lcen))       ; center frequency
(define df (* 0.2 fcen))     ; frequency width

(define pulse-src (list (make source
                          (src (make gaussian-src (frequency fcen) (fwidth df)))
                          (component Ez)
                          (center (+ (* -0.5 sx) dpml (* 0.5 dsub)) 0 0)
                          (size 0 sy 0))))

(set! sources pulse-src)

(set! k-point (vector3 0 0 0))

(define glass (make medium (index 1.5)))
(set! default-material glass)

(define symm (list (make mirror-sym (direction Y))))
(set! symmetries symm)

(define mon-pt (vector3 (- (* 0.5 sx) dpml (* 0.5 dpad)) 0 0))
(define flux-obj (add-flux fcen 0 1 (make flux-region (center mon-pt) (size 0 sy 0))))

(run-sources+ 50)

(define input-flux (get-fluxes flux-obj))

(reset-meep)

(set! geometry-lattice cell)

(set! pml-layers boundary-layers)

(set! sources pulse-src)

(set! k-point (vector3 0 0 0))

(set! default-material air)

(set! geometry (list (make block
                       (material glass)
                       (size (+ dpml dsub) infinity infinity)
                       (center (+ (* -0.5 sx) (* 0.5 (+ dpml dsub))) 0 0))
                     (make block
                       (material glass)
                       (size gh (* gdc gp) infinity)
                       (center (+ (* -0.5 sx) dpml dsub (* 0.5 gh)) 0 0))))

(set! symmetries symm)

(set! flux-obj (add-flux fcen 0 1 (make flux-region (center mon-pt) (size 0 sy 0))))

(run-sources+ 200)

(define res (get-eigenmode-coefficients flux-obj (list 1) #:eig-parity (+ ODD-Z EVEN-Y)))
(define coeffs (list-ref res 0))

(define mode-tran (/ (sqr (magnitude (array-ref coeffs 0 0 0))) (list-ref input-flux 0)))
(define mode-phase (angle (array-ref coeffs 0 0 0)))
(if (> mode-phase 0) (set! mode-phase (- mode-phase (* 2 pi))))

(print "mode:, " mode-tran ", " mode-phase "\n")
