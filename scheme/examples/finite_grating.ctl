;; true:  plot the scattered fields in the air region adjacent to the grating
;; false: plot the diffraction spectra based on a 1d cross section of the scattered fields
(define-param field-profile? true)

(set-param! resolution 50)                   ; pixels/μm

(define-param dpml 1.0)                      ; PML thickness
(define-param dsub 2.0)                      ; substrate thickness
(define-param dpad 1.0)                      ; flat-surface padding
(define-param dair                           ; air region thickness adjacent to grating
  (if field-profile? 10 dpad))
(define-param gp 1.0)                        ; grating periodicity
(define-param gh 0.5)                        ; grating height
(define-param gdc 0.5)                       ; grating duty cycle
(define-param num-cells 5)                   ; number of grating unit cells

(define-param wvl 0.5)                       ; center wavelength
(define fcen (/ wvl))                        ; center frequency

(set! k-point (vector3 0))

(define glass (make medium (index 1.5)))

(set! pml-layers (list (make pml (thickness dpml))))

(set! symmetries (list (make mirror-sym (direction Y))))

(define sx (+ dpml dsub gh dair dpml))
(define sy (+ dpml dpad (* num-cells gp) dpad dpml))

(set! geometry-lattice (make lattice (size sx sy no-size)))

(define src-pt (vector3 (+ (* -0.5 sx) dpml (* 0.5 dsub))))
(define pw-source (list (make source
                          (src (make gaussian-src (frequency fcen) (fwidth (* 0.2 fcen)) (is-integrated? true)))
                          (component Ez)
                          (center src-pt)
                          (size 0 sy))))
(set! sources pw-source)

(set! geometry (list (make block
                       (material glass)
                       (size (+ dpml dsub) infinity infinity)
                       (center (+ (* -0.5 sx) (* 0.5 (+ dpml dsub)))))))

(define mon-pt (vector3 (- (* 0.5 sx) dpml (* 0.5 dair))))
(define flat-fields (add-dft-fields (list Ez) fcen fcen 1 (volume (center mon-pt) (size (if field-profile? dair 0) (- sy (* 2 dpml))))))

(run-sources+ 100)

(output-dft flat-fields "flat")

(reset-meep)

(set! pml-layers (list (make pml (thickness dpml))))

(set! symmetries (list (make mirror-sym (direction Y))))

(set! geometry-lattice (make lattice (size sx sy no-size)))

(set! k-point (vector3 0))

(set! sources pw-source)

(set! geometry (list (make block
                       (material glass)
                       (size (+ dpml dsub) infinity infinity)
                       (center (+ (* -0.5 sx) (* 0.5 (+ dpml dsub)))))))

(set! geometry (append geometry
                       (map (lambda (n)
                              (make block
                                (material glass)
                                (size gh (* gdc gp) infinity)
                                (center (+ (* -0.5 sx) dpml dsub (* 0.5 gh)) (+ (* -0.5 sy) dpml dpad (* (+ n 0.5) gp)) 0)))
                            (arith-sequence 0 1 num-cells))))

(define grating-fields (add-dft-fields (list Ez) fcen fcen 1 (volume (center mon-pt) (size (if field-profile? dair 0) (- sy (* 2 dpml))))))

(run-sources+ 100)

(output-dft grating-fields "grating")
