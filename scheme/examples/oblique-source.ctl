(set-param! resolution 50) ; pixels/μm

(set! geometry-lattice (make lattice (size 14 14 no-size)))

(set! pml-layers (list (make pml (thickness 2))))

; rotation angle (in degrees) of waveguide, counter clockwise (CCW) around z-axis
(define-param rot-angle 20)
(set! rot-angle (deg->rad rot-angle))

; width of waveguide
(define-param w 1.0)

(set! geometry (list (make block
                       (center 0 0 0)
                       (size infinity w infinity)
                       (e1 (rotate-vector3 (vector3 0 0 1) rot-angle (vector3 1 0 0)))
                       (e2 (rotate-vector3 (vector3 0 0 1) rot-angle (vector3 0 1 0)))
                       (material (make medium (epsilon 12))))))

(define-param fsrc 0.15) ; frequency of eigenmode or constant-amplitude source
(define-param bnum 1)    ; band number of eigenmode

(define kpoint (rotate-vector3 (vector3 0 0 1) rot-angle (vector3 1 0 0)))

(define-param compute-flux? true) ; compute flux (true) or output the field profile (false)

(define-param eig-src? true)      ; eigenmode (true) or constant-amplitude (false) source

(set! sources (list
               (if eig-src?
                   (make eigenmode-source
                     (src (if compute-flux? (make gaussian-src (frequency fsrc) (fwidth (* 0.2 fsrc))) (make continuous-src (frequency fsrc))))
                     (center 0 0 0)
                     (size 0 (* 3 w) 0)
                     (direction (if (= rot-angle 0) AUTOMATIC NO-DIRECTION))
                     (eig-kpoint kpoint)
                     (eig-band bnum)
                     (eig-parity (if (= rot-angle 0) (+ EVEN-Y ODD-Z) ODD-Z))
                     (eig-match-freq? true))
                   (make source
                     (src (if compute-flux? (make gaussian-src (frequency fsrc) (fwidth (* 0.2 fsrc))) (make continuous-src (frequency fsrc))))
                     (center 0 0 0)
                     (size 0 (* 3 w) 0)
                     (component Ez)))))

(if (= rot-angle 0)
    (set! symmetries (list (make mirror-sym (direction Y)))))

(if compute-flux?
    (let ((tran (add-flux fsrc 0 1 (make flux-region (center 5 0 0) (size 0 14 0)))))
      (run-sources+ 50)
      (display-fluxes tran)
      (let ((res (get-eigenmode-coefficients tran
                                             (list 1)
                                             #:eig-parity (if (= rot-angle 0) (+ ODD-Z EVEN-Y) ODD-Z)
                                             #:direction NO-DIRECTION
                                             #:kpoint-func (lambda (f n) kpoint))))
        (print "mode-coeff-flux:, " (sqr (magnitude (array-ref (list-ref res 0) 0 0 0))) "\n")))
    (run-until 100 (in-volume (volume (center 0 0 0) (size 10 10 0))
                              (at-beginning output-epsilon)
                              (at-end output-efield-z))))
