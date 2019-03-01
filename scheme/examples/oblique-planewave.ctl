(set-param! resolution 50) ; pixels/μm

(set! geometry-lattice (make lattice (size 14 10 no-size)))

(set! pml-layers (list (make pml (thickness 2) (direction X))))

; rotation angle (in degrees) of planewave, counter clockwise (CCW) around z-axis
(define-param rot-angle 0)
(set! rot-angle (deg->rad rot-angle))

(define-param fsrc 1.0) ; frequency of planewave (wavelength = 1/fsrc)

(define-param n 1.5) ; refractive index of homogeneous material
(set! default-material (make medium (index n)))

(define k (rotate-vector3 (vector3 0 0 1) rot-angle (vector3 (* fsrc n) 0 0)))
(set! k-point k)

(if (= rot-angle 0)
    (set! symmetries (list (make mirror-sym (direction Y)))))

(set! sources (list
               (make eigenmode-source
                 (src (make continuous-src (frequency fsrc)))
                 (center 0 0 0)
                 (size 0 10 0)
                 (direction (if (= rot-angle 0) AUTOMATIC NO-DIRECTION))
                 (eig-kpoint k)
                 (eig-band 1)
                 (eig-parity (if (= rot-angle 0) (+ EVEN-Y ODD-Z) ODD-Z))
                 (eig-match-freq? true))))

(run-until 100 (in-volume (volume (center 0 0 0) (size 10 10 0))
                          (at-end output-efield-z)))
