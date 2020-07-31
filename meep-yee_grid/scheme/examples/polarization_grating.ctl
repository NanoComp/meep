;; polarization grating from C. Oh and M.J. Escuti, Optics Letters, Vol. 33, No. 20, pp. 2287-9, 2008
;; note: reference uses z as the propagation direction and y as the out-of-plane direction; this script uses x and z, respectively

(set-param! resolution 50)    ; pixels/μm

(define-param dpml 1.0)       ; PML thickness
(define-param dsub 1.0)       ; substrate thickness
(define-param dpad 1.0)       ; padding thickness

(set! k-point (vector3 0 0 0))

(define boundary-layers (list (make pml (thickness dpml) (direction X))))
(set! pml-layers boundary-layers)

(define n0 1.55)
(define delta-n 0.159)
(define eps-diag (matrix3x3 (vector3 (sqr n0) 0 0)
                            (vector3 0 (sqr n0) 0)
                            (vector3 0 0 (sqr (+ n0 delta-n)))))

(define-param wvl 0.54)      ; center wavelength
(define fcen (/ wvl))        ; center frequency

(define-param d 1.7)         ; chiral layer thickness
(define-param ph 70)         ; chiral layer twist angle
(define-param gp 6.5)        ; grating period
(define-param nmode 5)       ; number of mode coefficients to compute

(set! ph (deg->rad ph))

(define sx (+ dpml dsub d d dpad dpml))
(define sy gp)

(define cell (make lattice (size sx sy no-size)))
(set! geometry-lattice cell)

; linear-polarized planewave pulse source
(define src-pt (vector3 (+ (* -0.5 sx) dpml (* 0.3 dsub)) 0 0))
(define lp-src (list (make source
                       (src (make gaussian-src (frequency fcen) (fwidth (* 0.05 fcen))))
                       (component Ez)
                       (center src-pt)
                       (size 0 sy 0))
                     (make source
                       (src (make gaussian-src (frequency fcen) (fwidth (* 0.05 fcen))))
                       (component Ey)
                       (center src-pt)
                       (size 0 sy 0))))

(set! sources lp-src)

(set! default-material (make medium (index n0)))

(define tran-pt (vector3 (- (* 0.5 sx) dpml (* 0.5 dpad)) 0 0))
(define tran-flux1 (add-flux fcen 0 1 (make flux-region (center tran-pt) (size 0 sy 0))))

(run-sources+ 100)

(define input-flux (get-fluxes tran-flux1))

(reset-meep)

(set! geometry-lattice cell)

(set! pml-layers boundary-layers)

(set! sources lp-src)

(set! k-point (vector3 0 0 0))

(set! default-material air)

; twist angle of nematic director; from equation 1b
(define phi (lambda (p)
              (let ((xx (- (vector3-x p) (+ (* -0.5 sx) dpml dsub))))
                (if (and (>= xx 0) (<= xx d))
                    (+ (* pi (vector3-y p) (/ gp)) (* ph xx (/ d)))
                    (+ (* pi (vector3-y p) (/ gp)) (- (* ph xx (/ d))) (* 2 ph))))))

(define lc-epsilon-diag (vector3 0 0 0))
(define lc-epsilon-offdiag (vector3 0 0 0))
(define lc-epsilon (matrix3x3 (vector3 0 0 0) (vector3 0 0 0) (vector3 0 0 0)))

; return the anisotropic permittivity tensor for a uniaxial, twisted nematic liquid crystal
(define lc-mat (lambda (p)
                 (let
                     ; rotation matrix for rotation around x axis
                     ((Rx (matrix3x3 (vector3 1 0 0)
                                     (vector3 0 (cos (phi p)) (sin (phi p)))
                                     (vector3 0 (- (sin (phi p))) (cos (phi p))))))
                   (set! lc-epsilon (matrix3x3* Rx (matrix3x3* eps-diag (matrix3x3-transpose Rx))))
                   (set! lc-epsilon-diag (vector3 (vector3-x (vector3-x lc-epsilon))
                                                  (vector3-y (vector3-y lc-epsilon))
                                                  (vector3-z (vector3-z lc-epsilon))))
                   (set! lc-epsilon-offdiag (vector3 (vector3-x (vector3-y lc-epsilon))
                                                     (vector3-x (vector3-z lc-epsilon))
                                                     (vector3-y (vector3-z lc-epsilon))))
                   (make medium (epsilon-diag lc-epsilon-diag) (epsilon-offdiag lc-epsilon-offdiag)))))

(set! geometry (list
                (make block
                  (center (+ (* -0.5 sx) (* 0.5 (+ dpml dsub))) 0 0)
                  (size (+ dpml dsub) infinity infinity)
                  (material (make medium (index n0))))
                (make block
                  (center (+ (* -0.5 sx) dpml dsub d) 0 0)
                  (size (* 2 d) infinity infinity)
                  (material (make material-function (material-func lc-mat))))))

(define tran-flux2 (add-flux fcen 0 1 (make flux-region (center tran-pt) (size 0 sy 0))))

(run-sources+ 300)

(define res1 (get-eigenmode-coefficients tran-flux2 (arith-sequence 1 1 nmode) #:eig-parity (+ ODD-Z EVEN-Y)))
(define res2 (get-eigenmode-coefficients tran-flux2 (arith-sequence 1 1 nmode) #:eig-parity (+ EVEN-Z ODD-Y)))

(define t-coeffs1 (list-ref res1 0))
(define t-coeffs2 (list-ref res2 0))
(define kdom (list-ref res1 3))

(map (lambda (nm)
       (let ((mode-angle (acos (/ (vector3-x (list-ref kdom nm)) fcen)))
             (mode-tran (/ (+ (sqr (magnitude (array-ref t-coeffs1 nm 0 0))) (sqr (magnitude (array-ref t-coeffs2 nm 0 0)))) (list-ref input-flux 0))))
         (print "tran:, " nm ", " (rad->deg mode-angle) ", " mode-tran "\n")))
     (arith-sequence 0 1 nmode))
