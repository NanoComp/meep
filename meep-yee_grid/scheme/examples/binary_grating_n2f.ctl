(set-param! resolution 25)                   ; pixels/μm

(define-param dpml 1.0)                      ; PML thickness
(define-param dsub 3.0)                      ; substrate thickness
(define-param dpad 3.0)                      ; padding between grating and pml
(define-param gp 10.0)                       ; grating period
(define-param gh 0.5)                        ; grating height
(define-param gdc 0.5)                       ; grating duty cycle

(define-param nperiods 10)                   ; number of unit cells in finite periodic grating

(define-param ff-distance 1e8)               ; far-field distance from near-field monitor
(define-param ff-angle 20)                   ; far-field cone angle
(define-param ff-npts 500)                   ; number of far-field points

(define ff-length (* ff-distance (tan (deg->rad ff-angle))))
(define ff-res (/ ff-npts ff-length))

(define sx (+ dpml dsub gh dpad dpml))

(define-param wvl-min 0.4)                   ; min wavelength
(define-param wvl-max 0.6)                   ; max wavelength
(define fmin (/ wvl-max))                    ; min frequency
(define fmax (/ wvl-min))                    ; max frequency
(define fcen (* 0.5 (+ fmin fmax)))          ; center frequency
(define df (- fmax fmin))                    ; frequency width

(define glass (make medium (index 1.5)))

(set! geometry-lattice (make lattice (size sx no-size no-size)))

(set! pml-layers (list (make pml (thickness dpml) (direction X))))

(set! k-point (vector3 0))

(set! default-material glass)

(define src-pt (vector3 (+ (* -0.5 sx) dpml (* 0.5 dsub))))
(set! sources (list (make source
          (src (make gaussian-src (frequency fcen) (fwidth df)))
          (component Ez) (center src-pt))))

(define-param nfreq 21)
(define n2f-pt (vector3 (- (* 0.5 sx) dpml (* 0.5 dpad))))
(define n2f-obj (add-near2far fcen df nfreq (make near2far-region (center n2f-pt))))

(run-sources+ (stop-when-fields-decayed 50 Ez n2f-pt 1e-9))

(output-farfields n2f-obj "source" (volume (center ff-distance (* 0.5 ff-length)) (size 0 ff-length)) ff-res)

(reset-meep)

;;; unit cell with periodic boundaries

(define sy gp)
(set! geometry-lattice (make lattice (size sx sy no-size)))

(set! pml-layers (list (make pml (thickness dpml) (direction X))))

(set! sources (list (make source
          (src (make gaussian-src (frequency fcen) (fwidth df)))
          (component Ez) (center src-pt) (size 0 sy))))

(set! default-material air)

(set! geometry (list
                (make block (material glass) (size (+ dpml dsub) infinity infinity) (center (+ (* -0.5 sx) (* 0.5 (+ dpml dsub)))))
                (make block (material glass) (size gh (* gdc gp) infinity) (center (+ (* -0.5 sx) dpml dsub (* 0.5 gh))))))

(set! k-point (vector3 0))

(set! symmetries (list (make mirror-sym (direction Y))))

(set! n2f-obj (add-near2far fcen df nfreq (make near2far-region (center n2f-pt) (size 0 sy)) #:nperiods nperiods))

(run-sources+ (stop-when-fields-decayed 50 Ez n2f-pt 1e-9))

(output-farfields n2f-obj "unit-cell" (volume (center ff-distance (* 0.5 ff-length)) (size 0 ff-length)) ff-res)

(reset-meep)

;;; finite periodic grating with flat surface termination extending into PML

(define num-cells (+ (* 2 nperiods) 1))
(set! sy (+ dpml (* num-cells gp) dpml))
(set! geometry-lattice (make lattice (size sx sy no-size)))

(set! pml-layers (list (make pml (thickness dpml))))

(set! sources (list (make source
          (src (make gaussian-src (frequency fcen) (fwidth df) (is-integrated? true)))
          (component Ez) (center src-pt) (size 0 sy))))

(set! geometry (list (make block (material glass) (size (+ dpml dsub) infinity infinity) (center (+ (* -0.5 sx) (* 0.5 (+ dpml dsub)))))))

(set! geometry (append geometry
                       (map (lambda (n)
                              (make block (material glass) (size gh (* gdc gp) infinity)
                                    (center (+ (* -0.5 sx) dpml dsub (* 0.5 gh)) (+ (* -0.5 sy) dpml (* (+ n 0.5) gp)) 0)))
                            (arith-sequence 0 1 num-cells))))

(set! k-point (vector3 0))

(set! symmetries (list (make mirror-sym (direction Y))))

(set! n2f-obj (add-near2far fcen df nfreq (make near2far-region (center n2f-pt) (size 0 (- sy (* 2 dpml))))))

(run-sources+ (stop-when-fields-decayed 50 Ez n2f-pt 1e-9)
              (at-beginning output-epsilon))

(output-farfields n2f-obj "super-cell" (volume (center ff-distance (* 0.5 ff-length)) (size 0 ff-length)) ff-res)
