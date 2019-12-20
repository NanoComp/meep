(set-param! resolution 30)  ; pixels/um

(define-param perpendicular? true)
(define src-cmpt (if perpendicular? Hz Ez))
(define fcen (if perpendicular? 0.21 0.17))  ; pulse center frequency

(define-param n 3.4)  ; index of waveguide
(define-param w 1)    ; width of waveguide
(define-param r 1)    ; inner radius of ring
(define-param pad 4)  ; padding between waveguide and edge of PML
(define-param dpml 2) ; thickness of PML
(set-param! m 5)      ; angular dependence

(set! pml-layers (list (make pml (thickness dpml))))

(define sxy (* 2 (+ r w pad dpml)))
(set! geometry-lattice (make lattice (size sxy sxy no-size)))

(set! symmetries (list (make mirror-sym (direction X) (phase (if perpendicular? +1 -1)))
                       (make mirror-sym (direction Y) (phase (if perpendicular? -1 +1)))))

(set! geometry (list
		(make cylinder
                  (material (make dielectric (index n)))
                  (radius (+ r w))
                  (height infinity)
                  (center 0 0))
		(make cylinder
                  (material vacuum)
                  (radius r)
                  (height infinity)
                  (center 0 0))))

;; find resonant frequency of unperturbed geometry using broadband source

(define df (* 0.2 fcen))  ; pulse width (in frequency)

(set! sources (list
               (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (component src-cmpt)
                 (center (+ r 0.1)))
               (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (component src-cmpt)
                 (center (- (+ r 0.1)))
                 (amplitude -1))))

(run-sources+ 100 (after-sources (harminv src-cmpt (vector3 (+ r 0.1)) fcen df)))

(define frq-unperturbed (harminv-freq-re (car harminv-results)))

(reset-meep)

;; unperturbed geometry with narrowband source

(set! pml-layers (list (make pml (thickness dpml))))

(set! geometry-lattice (make lattice (size sxy sxy no-size)))

(set! geometry (list
		(make cylinder
                  (material (make dielectric (index n)))
                  (radius (+ r w))
                  (height infinity)
                  (center 0 0))
		(make cylinder
                  (material vacuum)
                  (radius r)
                  (height infinity)
                  (center 0 0))))

(set! fcen frq-unperturbed)
(set! df (* 0.05 fcen))

(set! sources (list
               (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (component src-cmpt)
                 (center (+ r 0.1)))
               (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (component src-cmpt)
                 (center (- (+ r 0.1)))
                 (amplitude -1))))

(run-sources+ 100)

(define deps (- 1 (* n n)))
(define deps-inv (- 1 (/ (* n n))))

(define numerator-integral 0)

(if perpendicular?
    (let ((para-integral (* deps 2 pi (- (* r (sqr (magnitude (get-field-point Ey (vector3 r))))) (* (+ r w) (sqr (magnitude (get-field-point Ey (vector3 (+ r w)))))))))
          (perp-integral (* deps-inv 2 pi (- (* (+ r w) (sqr (magnitude (get-field-point Dy (vector3 0 (+ r w)))))) (* r (sqr (magnitude (get-field-point Dy (vector3 0 r)))))))))
      (set! numerator-integral (+ para-integral perp-integral)))
    (set! numerator-integral (* deps 2 pi (- (* r (sqr (magnitude (get-field-point Ez (vector3 r))))) (* (+ r w) (sqr (magnitude (get-field-point Ez (vector3 (+ r w))))))))))

(define denominator-integral (electric-energy-in-box (volume (center 0 0) (size (- sxy (* 2 dpml)) (- sxy (* 2 dpml))))))
(define perturb-theory-dw-dR (/ (* -1 frq-unperturbed numerator-integral) (* 8 denominator-integral)))

(reset-meep)

;; perturbed geometry with narrowband source

(define-param dr 0.04)

(set! pml-layers (list (make pml (thickness dpml))))

(set! geometry-lattice (make lattice (size sxy sxy no-size)))

(set! geometry (list
		(make cylinder
                  (material (make dielectric (index n)))
                  (radius (+ r dr w))
                  (height infinity)
                  (center 0 0))
		(make cylinder
                  (material vacuum)
                  (radius (+ r dr))
                  (height infinity)
                  (center 0 0))))

(set! sources (list
               (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (component src-cmpt)
                 (center (+ r dr 0.1)))
               (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (component src-cmpt)
                 (center (- (+ r dr 0.1)))
                 (amplitude -1))))

(run-sources+ 100 (after-sources (harminv src-cmpt (vector3 (+ r 0.1)) fcen df)))

(define frq-perturbed (harminv-freq-re (car harminv-results)))

(define finite-diff-dw-dR (/ (- frq-perturbed frq-unperturbed) dr))

(print "dwdR:, " perturb-theory-dw-dR " (pert. theory), " finite-diff-dw-dR " (finite diff.)\n")
