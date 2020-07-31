;; diffraction spectra of a binary phase grating from the Meep tutorial

(set-param! resolution 60)  ; pixels/μm

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

(define-param wvl-min 0.4)          ; min wavelength
(define-param wvl-max 0.6)          ; max wavelength
(define fmin (/ wvl-max))           ; min frequency
(define fmax (/ wvl-min))           ; max frequency
(define fcen (* 0.5 (+ fmin fmax))) ; pulse frequency center
(define df (- fmax fmin))           ; pulse frequency width

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

(define-param nfreq 21)
(define mon-pt (vector3 (- (* 0.5 sx) dpml (* 0.5 dpad)) 0 0))
(define flux-mon (add-flux fcen df nfreq (make flux-region (center mon-pt) (size 0 sy 0))))

(run-sources+ (stop-when-fields-decayed 50 Ez mon-pt 1e-9))

(define input-flux (get-fluxes flux-mon))
(define freqs (get-flux-freqs flux-mon))

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

(define mode-mon (add-flux fcen df nfreq (make flux-region (center mon-pt) (size 0 sy 0))))

(run-sources+ (stop-when-fields-decayed 50 Ez mon-pt 1e-9))

(define-param nmode 10)
(define res (get-eigenmode-coefficients mode-mon (arith-sequence 1 1 nmode) #:eig-parity (+ ODD-Z EVEN-Y)))
(define coeffs (list-ref res 0))
(define kdom (list-ref res 3))

(map (lambda (nm)
       (map (lambda (nf)
              (let ((mode-wvl (/ (list-ref freqs nf)))
                    (mode-angle (rad->deg (acos (/ (vector3-x (list-ref kdom (+ (* nm nfreq) nf))) (list-ref freqs nf)))))
                    (mode-tran (/ (sqr (magnitude (array-ref coeffs nm nf 0))) (list-ref input-flux nf))))
                (if (> nm 0) (set! mode-tran (* 0.5 mode-tran)))
                (print "grating" nm ":, " mode-wvl ", " mode-angle ", " mode-tran "\n")))
              (arith-sequence 0 1 nfreq)))
     (arith-sequence 0 1 nmode))
