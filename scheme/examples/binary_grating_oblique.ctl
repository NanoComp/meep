;; reflectance and transmittance spectra for planewave at oblique incidence
;; of a binary phase grating from the Meep tutorial

(set-param! resolution 50)  ; pixels/μm

(define-param dpml 1.0)     ; PML thickness
(define-param dsub 3.0)     ; substrate thickness
(define-param dpad 3.0)     ; padding between grating and PML
(define-param gp 10.0)      ; grating period
(define-param gh 0.5)       ; grating height
(define-param gdc 0.5)      ; grating duty cycle

(define sx (+ dpml dsub gh dpad dpml))
(define sy gp)

(define cell (make lattice (size sx sy no-size)))
(set! geometry-lattice cell)

(define boundary-layers (list (make pml (thickness dpml) (direction X))))
(set! pml-layers boundary-layers)

(define-param wvl-cen 0.5)  ; center wavelength
(define fcen (/ wvl-cen))   ; center frequency
(define df (* 0.05 fcen))   ; frequency width

(define ng 1.5)
(define glass (make medium (index ng)))
(set! default-material glass)

(define-param use-cw-solver? false)       ; CW solver or time stepping?
(define-param cw-solver-tol 1e-6)         ; CW solver tolerance
(define-param cw-solver-max-iters 2000)   ; CW solver max iterations
(define-param cw-solver-L 10)             ; CW solver L

; rotation angle of incident planewave; counter clockwise (CCW) about Z axis, 0 degrees along +X axis
(define-param theta-in 10.7)
(set! theta-in (deg->rad theta-in))

; k (in source medium) with correct length (plane of incidence: XY)
(define k (rotate-vector3 (vector3 0 0 1) theta-in (vector3 (* fcen ng) 0 0)))

(define symm '())
(define eig-parity ODD-Z)
(if (= theta-in 0)
    (begin
      (set! k (vector3 0 0 0))
      (set! symm (list (make mirror-sym (direction Y))))
      (set! eig-parity (+ eig-parity EVEN-Y))))

(set! k-point k)

(set! symmetries symm)

(define (pw-amp k x0)
  (lambda (x)
    (exp (* 0+1i 2 pi (vector3-dot k (vector3- x x0))))))

(define src-pt (vector3 (+ (* -0.5 sx) dpml (* 0.3 dsub)) 0 0))
(define pw-src (list (make source
                       (if use-cw-solver?
                           (src (make continuous-src (frequency fcen) (fwidth df)))
                           (src (make gaussian-src (frequency fcen) (fwidth df))))
                       (component Ez)
                       (center src-pt)
                       (size 0 sy 0)
                       (amp-func (pw-amp k src-pt)))))

(set! sources pw-src)

(define refl-pt (vector3 (+ (* -0.5 sx) dpml (* 0.5 dsub)) 0 0))
(define refl-flux (add-flux fcen 0 1 (make flux-region (center refl-pt) (size 0 sy 0))))

(if use-cw-solver?
    (begin
      (init-fields)
      (meep-fields-solve-cw fields cw-solver-tol cw-solver-maxiters cw-solver-L))
    (run-sources+ 100))

(save-flux "flux" refl-flux)
(define input-flux (get-fluxes refl-flux))
(define freqs (get-flux-freqs refl-flux))

(reset-meep)

(set! geometry-lattice cell)

(set! pml-layers boundary-layers)

(set! sources pw-src)

(set! k-point k)

(set! symmetries symm)

(set! default-material air)

(set! geometry (list (make block
                       (material glass)
                       (size (+ dpml dsub) infinity infinity)
                       (center (+ (* -0.5 sx) (* 0.5 (+ dpml dsub))) 0 0))
                     (make block
                       (material glass)
                       (size gh (* gdc gp) infinity)
                       (center (+ (* -0.5 sx) dpml dsub (* 0.5 gh)) 0 0))))

(set! refl-flux (add-flux fcen 0 1 (make flux-region (center refl-pt) (size 0 sy 0))))
(load-minus-flux "flux" refl-flux)

(define tran-pt (vector3 (- (* 0.5 sx) dpml (* 0.5 dpad)) 0 0))
(define tran-flux (add-flux fcen 0 1 (make flux-region (center tran-pt) (size 0 sy 0))))

(if use-cw-solver?
    (begin
      (init-fields)
      (meep-fields-solve-cw fields cw-solver-tol cw-solver-maxiters cw-solver-L))
    (run-sources+ 200))

; number of reflected orders
(define nm-r (- (floor (* (- (* fcen ng) (vector3-y k)) gp)) (ceiling (* (- (- (* fcen ng)) (vector3-y k)) gp))))
(if (= theta-in 0) (set! nm-r (* 0.5 nm-r)))

(define res (get-eigenmode-coefficients refl-flux (arith-sequence 1 1 nm-r) #:eig-parity eig-parity))
(define r-coeffs (list-ref res 0))
(define kdom (list-ref res 3))

(define Rsum 0)
(define r-angle 0)
(map (lambda (nm)
       (let ((r-kdom (list-ref kdom nm))
             (Rmode (/ (sqr (magnitude (array-ref r-coeffs nm 0 1))) (list-ref input-flux 0))))
         (set! r-angle (* (if (positive? (vector3-y r-kdom)) +1 -1) (acos (/ (vector3-x r-kdom) (* ng fcen)))))
         (print "refl:, " nm ", " (rad->deg r-angle) ", " Rmode "\n")
         (set! Rsum (+ Rsum Rmode))))
     (arith-sequence 0 1 (- nm-r 1)))

; number of transmitted orders
(define nm-t (- (floor (* (- fcen (vector3-y k)) gp)) (ceiling (* (- (- fcen) (vector3-y k)) gp))))
(if (= theta-in 0) (set! nm-t (* 0.5 nm-t)))

(set! res (get-eigenmode-coefficients tran-flux (arith-sequence 1 1 nm-t) #:eig-parity eig-parity))
(define t-coeffs (list-ref res 0))
(set! kdom (list-ref res 3))

(define Tsum 0)
(define t-angle 0)
(map (lambda (nm)
       (let ((t-kdom (list-ref kdom nm))
             (Tmode (/ (sqr (magnitude (array-ref t-coeffs nm 0 0))) (list-ref input-flux 0))))
         (set! t-angle (* (if (positive? (vector3-y t-kdom)) +1 -1) (acos (/ (vector3-x t-kdom) fcen))))
         (print "tran:, " nm ", " (rad->deg t-angle) ", " Tmode "\n")
         (set! Tsum (+ Tsum Tmode))))
     (arith-sequence 0 1 (- nm-t 1)))

(print "mode-coeff:, " Rsum ", " Tsum ", " (+ Rsum Tsum) "\n")

(define r-flux (get-fluxes refl-flux))
(define t-flux (get-fluxes tran-flux))

(define Rflux (/ (- (list-ref r-flux 0)) (list-ref input-flux 0)))
(define Tflux (/ (list-ref t-flux 0) (list-ref input-flux 0)))

(print "poynting-flux:, "  Rflux ", " Tflux ", " (+ Rflux Tflux) "\n")
