(set-param! resolution 50) ; pixels/μm

(set! geometry-lattice (make lattice (size 14 14 no-size)))

(set! pml-layers (list (make pml (thickness 2))))

(define-param rot-angle 0) ; rotation angle (in degrees) of waveguide, CCW around z-axis
(set! rot-angle (deg->rad rot-angle))

(set! geometry (list (make block
                       (center 0 0 0)
                       (size infinity 1 infinity)
                       (e1 (rotate-vector3 (vector3 0 0 1) rot-angle (vector3 1 0 0)))
                       (e2 (rotate-vector3 (vector3 0 0 1) rot-angle (vector3 0 1 0)))
                       (material (make medium (epsilon 12))))))

(define-param fsrc 0.15) ; frequency of eigenmode or continuous-wave (CW) source
(define-param kx 0.4)    ; initial guess for wavevector in x-direction of eigenmode 
(define-param bnum 1)    ; band index of eigenmode

(define-param eig-src? true)

(set! sources (list
               (if eig-src?
                   (make eigenmode-source
                     (src (make continuous-src (frequency fsrc)))
                     (center 0 0 0)
                     (size 0 14 0)
                     (direction (if (= rot-angle 0) AUTOMATIC NO-DIRECTION))
                     (eig-kpoint (rotate-vector3 (vector3 0 0 1) rot-angle (vector3 kx 0 0)))
                     (eig-band bnum)
                     (eig-parity (if (= rot-angle 0) (+ EVEN-Y ODD-Z) ODD-Z))
                     (eig-match-freq? true))
                   (make source
                     (src (make continuous-src (frequency fsrc)))
                     (center 0 0 0)
                     (size 0 2 0)
                     (component Ez)))))

(if (= rot-angle 0)
    (set! symmetries (list (make mirror-sym (direction Y)))))

(run-until 100 (in-volume (volume (center 0 0 0) (size 10 10 0))
                          (at-beginning output-epsilon)
                          (at-end output-efield-z)))
