;; far fields at focal-length distance of a metalens with binary-phase grating unit cell from the Meep tutorial

(set-param! resolution 50)  ; pixels/μm

(define-param dpml 1.0)     ; PML thickness
(define-param dsub 2.0)     ; substrate thickness
(define-param dpad 2.0)     ; padding between grating and PML
(define-param gp 0.3)       ; grating period
(define-param gh 1.8)       ; grating height

(define-param focal-length 200) ; focal length of metalens
(define-param spot-length 100)  ; far field line length
(define-param ff-res 10)        ; far field resolution (points/μm)

; list of grating duty cycles
(define-param gdc-list (list '()))

; # of cells
(define num-cells (length gdc-list))

; return gdc of nth cell
(define gdc-cell (lambda (n) (list-ref gdc-list n)))

(define sx (+ dpml dsub gh dpad dpml))
(define sy (* num-cells gp))

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

(set! geometry (list (make block
                       (material glass)
                       (size (+ dpml dsub) infinity infinity)
                       (center (+ (* -0.5 sx) (* 0.5 (+ dpml dsub))) 0 0))))

(set! geometry (append geometry
                       (map (lambda (n)
                              (make block
                                (material glass)
                                (size gh (* (gdc-cell n) gp) infinity)
                                (center (+ (* -0.5 sx) dpml dsub (* 0.5 gh)) (+ (* -0.5 sy) (* (+ n 0.5) gp)) 0)))
                            (arith-sequence 0 1 num-cells))))

(define symm (list (make mirror-sym (direction Y))))
(set! symmetries symm)

(define mon-pt (vector3 (- (* 0.5 sx) dpml (* 0.5 dpad)) 0 0))
(define n2f-obj (add-near2far fcen 0 1 (make near2far-region (center mon-pt) (size 0 sy 0))))

(run-sources+ 500)

(output-farfields n2f-obj (string-append "numcells-" (number->string num-cells)) (volume (center (+ (* -0.5 sx) dpml dsub gh focal-length) 0 0) (size spot-length 0 0)) ff-res)
