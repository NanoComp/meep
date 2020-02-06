(set-param! resolution 25)       ;; pixels/um

(define-param dpml 1.0)          ;; PML thickness
(define-param dsub 2.0)          ;; substrate thickness
(define-param dpad 2.0)          ;; padding between zone plate and PML
(define-param zh 0.5)            ;; zone-plate height
(define-param zN 25)             ;; number of zones (odd zones: pi phase shift, even zones: none)
(define-param focal-length 200)  ;; focal length of zone plate
(define-param spot-length 100)   ;; far-field line length
(define-param ff-res 10)         ;; far-field resolution (points/um)

(set! pml-layers (list (make pml (thickness dpml))))

(define-param wvl-cen 0.5)
(define frq-cen (/ wvl-cen))
(define dfrq (* 0.2 frq-cen))

(define r (map (lambda (n)
                 (sqrt (* n wvl-cen (+ focal-length (/ (* n wvl-cen) 4)))))
                 (arith-sequence 1 1 zN)))

(define sr (+ (list-ref r (- zN 1)) dpad dpml))
(define sz (+ dpml dsub zh dpad dpml))
(define cell (make lattice (size sr 0 sz)))
(set! geometry-lattice cell)
(set! dimensions CYLINDRICAL)
(set-param! m -1)

(set! sources (list (make source
                      (src (make gaussian-src (frequency frq-cen) (fwidth dfrq) (is-integrated? true)))
                      (center (* 0.5 sr) 0 (+ (* -0.5 sz) dpml))
                      (size sr 0 0)
                      (component Er))
                    (make source
                      (src (make gaussian-src (frequency frq-cen) (fwidth dfrq) (is-integrated? true)))
                      (center (* 0.5 sr) 0 (+ (* -0.5 sz) dpml))
                      (size sr 0 0)
                      (component Ep)
                      (amplitude 0-1i))))

(define glass (make medium (index 1.5)))

(set! geometry (list (make block
                       (material glass)
                       (size sr 0 (+ dpml dsub))
                       (center (* 0.5 sr) 0 (+ (* -0.5 sz) (* 0.5 (+ dpml dsub)))))))

(set! geometry (append geometry
                       (map (lambda (n)
                              (make block
                                (material (if (= (modulo n 2) 0) glass vacuum))
                                (size (list-ref r n) 0 zh)
                                (center (* 0.5 (list-ref r n)) 0 (+ (* -0.5 sz) dpml dsub (* 0.5 zh)))))
                            (reverse (arith-sequence 0 1 zN)))))

(define n2f-obj (add-near2far frq-cen 0 1
                              (make near2far-region (center (* 0.5 (- sr dpml)) 0 (- (* 0.5 sz) dpml)) (size (- sr dpml) 0 0))
                              (make near2far-region (center (- sr dpml) 0 (- (* 0.5 sz) (* 0.5 (+ dsub zh dpad)))) (size 0 0 (+ dsub zh dpad)))))

(run-sources+ 100)

(output-farfields n2f-obj "r" (volume (center (* 0.5 (- sr dpml)) 0 (+ (* -0.5 sz) dpml dsub zh focal-length)) (size (- sr dpml) 0 0)) ff-res)
(output-farfields n2f-obj "z" (volume (center 0 0 (+ (* -0.5 sz) dpml dsub zh focal-length)) (size 0 0 spot-length)) ff-res)
