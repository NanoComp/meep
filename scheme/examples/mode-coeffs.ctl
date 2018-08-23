(set-param! resolution 15)

(define-param w 1)  ; width of waveguide
(define-param L 10)  ; length of waveguide

(define Si (make medium (epsilon 12.0)))

(define-param dair 3.0)
(define-param dpml 3.0)

(define sx (+ dpml L dpml))
(define sy (+ dpml dair w dair dpml))
(set! geometry-lattice (make lattice (size sx sy no-size)))

(define prism_x (+ sx 1))
(define prism_y (/ w 2))
(define verts
  (list
   (vector3 (* prism_x -1) prism_y)
   (vector3 prism_x prism_y)
   (vector3 prism_x (* prism_y -1))
   (vector3 (* prism_x -1) (* prism_y -1))))

(set! geometry
      (list
       (make prism
         (center (vector3 0 0))
         (vertices verts)
         (height infinity)
         (material (make medium (epsilon 12.0))))))

(set! pml-layers (list (make pml (thickness dpml))))

; mode frequency
(define-param fcen 0.20)  ; > 0.5/sqrt(11) to have at least 2 modes

(define-param mode-num 1)

(set! sources (list
               (make eigenmode-source
                 (src (make gaussian-src (frequency fcen) (fwidth (* 0.5 fcen))))
                 (center (vector3 (+ (* -0.5 sx) dpml) 0))
                 (component ALL-COMPONENTS)
                 (size (vector3 0 (- sy (* 2 dpml))))
                 (eig-match-freq? true)
                 (eig-band mode-num)
                 (eig-resolution 32))))

(set! symmetries
      (list
       (make mirror-sym (direction Y) (phase (if (odd? mode-num) 1 -1)))))

(define xm (- (* 0.5 sx) dpml))  ; x-coordinate of monitor
(define mflux
  (add-mode-monitor fcen 0 1
                    (make mode-region (center (vector3 xm 0)) (size (vector3 0 (- sy (* 2 dpml)))))))
(define mode-flux
  (add-flux fcen 0 1
            (make flux-region (center (vector3 xm 0)) (size (vector3 0 (- sy (* 2 dpml)))))))
(run-sources+ 100)

(define modes-to-check '(1 2))  ; indices of modes for which to compute expansion coefficients
(define alpha-vgrp-kpoints (get-eigenmode-coefficients mflux modes-to-check))
(define alpha (first alpha-vgrp-kpoints))
(define vgrp (second alpha-vgrp-kpoints))
(define kpoints (third alpha-vgrp-kpoints))

;; self.assertTrue(kpoints[0].close(mp.Vector3(0.604301, 0, 0)))
;; self.assertTrue(kpoints[1].close(mp.Vector3(0.494353, 0, 0), tol=1e-2))

(define mode-power (car (get-fluxes mode-flux)))

(define test-passed #t)
(define tolerance 5.0e-3)
(define c0 (array-ref alpha (- mode-num 1) 0 0))  ; coefficient of forward-traveling wave for mode # mode_num

(map (lambda (nm)
       (if (not (= mode-num nm))
           (let ((cfrel (/ (magnitude (array-ref alpha (- nm 1) 0 0)) (magnitude c0)))
                 (cbrel (/ (magnitude (array-ref alpha (- nm 1) 0 1)) (magnitude c0))))
             (if (or (> cfrel tolerance) (> cbrel tolerance))
                 (set! test-passed #f)))))
     modes-to-check)

; test 1: coefficient of excited mode >> coeffs of all other modes
(if (not test-passed)
    (begin
      (print "Test failed")
      (exit 1)))

; test 2: |mode coeff|^2 = power
(if (> (abs (- 1 (/ mode-power (* (magnitude c0) (magnitude c0))))) 0.1)
    (begin
      (print "Test failed")
      (exit 1)))
