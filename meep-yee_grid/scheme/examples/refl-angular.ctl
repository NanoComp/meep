(set-param! resolution 200)          ; pixels/um

(define-param dpml 1)                ; PML thickness
(define-param sz 10)                 ; size of computational cell (without PMLs)
(set! sz (+ sz (* 2 dpml)))
(set! pml-layers (list (make pml (thickness dpml))))

(set! geometry-lattice (make lattice (size no-size no-size sz)))

(define-param wvl-min 0.4)           ; min wavelength
(define-param wvl-max 0.8)           ; max wavelength
(define fmin (/ wvl-max))            ; min frequency
(define fmax (/ wvl-min))            ; max frequency
(define fcen (* 0.5 (+ fmin fmax)))  ; center frequency
(define df (- fmax fmin))            ; frequency width
(define-param nfreq 50)              ; number of frequency bins

; rotation angle (in degrees) of source: CCW around Y axis, 0 degrees along +Z axis
(define-param theta 0)
(define theta-r (deg->rad theta))

; if normal incidence, force number of dimensions to be 1
(set! dimensions (if (= theta-r 0) 1 3))

; plane of incidence is xz
(set! k-point (vector3* fmin (vector3 (sin theta-r) 0 (cos theta-r))))

(set! sources (list (make source (src (make gaussian-src (frequency fcen) (fwidth df)))
			         (component Ex) (center 0 0 (+ (* -0.5 sz) dpml)))))

(define-param empty? true)

; add a block with n=3.5 for the air-dielectric interface
(if (not empty?)
    (set! geometry (list (make block (size infinity infinity (* 0.5 sz)) (center 0 0 (* 0.25 sz)) (material (make medium (index 3.5)))))))

(define refl (add-flux fcen df nfreq (make flux-region (center 0 0 (* -0.25 sz)))))

(if (not empty?) (load-minus-flux "refl-flux" refl))

(run-sources+ (stop-when-fields-decayed 50 Ex (vector3 0 0 (+ (* -0.5 sz) dpml)) 1e-9))

(if empty? (save-flux "refl-flux" refl))

(display-fluxes refl)
