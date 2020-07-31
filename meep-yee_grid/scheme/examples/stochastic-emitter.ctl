(set-param! resolution 50)     ;; resolution (pixels/um)

(define-param nr 20)           ;; number of random trials (method 1)
(define-param nd 10)           ;; number of dipoles
(define-param nf 500)          ;; number of frequencies
(define-param textured? false) ;; flat (default) or textured surface
(define-param method 1)        ;; type of method (1 or 2)

(define dpml 1.0)
(define dair 1.0)
(define hrod 0.7)
(define wrod 0.5)
(define dsub 5.0)
(define dAg 0.5)

(define sx 1.1)
(define sy (+ dpml dair hrod dsub dAg))
(set! geometry-lattice (make lattice (size sx sy no-size)))

(set! pml-layers (list (make pml (direction Y) (thickness dpml) (side High))))

(define fcen 1.0)
(define df 0.2)
(define run-time (* 2 (/ nf df)))

(set! geometry (list (make block
                       (material (make medium (index 3.45)))
                       (center 0 (- (* 0.5 sy) dpml dair hrod (* 0.5 dsub)))
                       (size infinity dsub infinity))
                     (make block
                       (material Ag)
                       (center 0 (+ (* -0.5 sy) (* 0.5 dAg)))
                       (size infinity dAg infinity))))

(if textured?
    (set! geometry (append geometry (list (make block
                                            (material (make medium (index 3.45)))
                                            (center 0 (- (* 0.5 sy) dpml dair (* 0.5 hrod)))
                                            (size wrod hrod infinity))))))

(set! k-point (vector3 0 0 0))

(define (compute-flux . args)
  (let ((m (get-keyword-value args #:m 1))
        (n (get-keyword-value args #:n 0)))
    (reset-meep)
    (if (= m 1)
        ;; method 1
        (map (lambda (nn)
               (set! sources (append sources (list (make source
                                                     (src (make custom-src (src-func (lambda (t) (random:normal)))))
                                                     (component Ez)
                                                     (center (* sx (+ -0.5 (/ nn nd))) (+ (* -0.5 sy) dAg (* 0.5 dsub))))))))
             (arith-sequence 0 1 nd))
        ;; method 2
        (set! sources (list (make source
                              (src (make gaussian-src (frequency fcen) (fwidth df)))
                              (component Ez)
                              (center (* sx (+ -0.5 (/ n nd))) (+ (* -0.5 sy) dAg (* 0.5 dsub)))))))
    (set! geometry-lattice geometry-lattice)
    (set! pml-layers pml-layers)
    (set! geometry geometry)
    (let ((flux-mon (add-flux fcen df nf (make flux-region (center 0 (- (* 0.5 sy) dpml)) (size sx 0)))))
      (run-until run-time)
      (display-fluxes flux-mon))))

(if (= method 1)
    (map (lambda (t)
           (compute-flux #:m 1))
         (arith-sequence 0 1 nr))
    (map (lambda (d)
           (compute-flux #:m 2 #:n d))
         (arith-sequence 0 1 nd)))
