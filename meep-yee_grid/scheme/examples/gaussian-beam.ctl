;; generate a titled Gaussian beam profile by defining the amplitude function of the source

(set-param! resolution 40) ; pixels/μm

(set! geometry-lattice (make lattice (size 20 10 no-size)))

(set! pml-layers (list (make pml (thickness 1.0) (direction Y))))

(define-param fcen 1.0) ; center frequency of CW source (wavelength is 1 μm)

(define-param tilt-angle -10)
(set! tilt-angle (deg->rad tilt-angle))

(define k (vector3-scale fcen (rotate-vector3 (vector3 0 0 1) tilt-angle (vector3 0 1 0))))

(define-param beam-sigma 1.5) ; beam width

(define (gaussian-beam sigma k x0)
  (lambda (x)
   (exp (- (* 0+2i pi (vector3-dot k (vector3- x x0)))
           (/ (vector3-dot (vector3- x x0) (vector3- x x0)) (* 2 sigma sigma))))))

(define src-pt (vector3 0 4 0))
(set! sources (list (make source
                      (src (make continuous-src (frequency fcen) (fwidth (* 0.2 fcen))))
                      (component Ez)
                      (center src-pt)
                      (size 20 0 0)
                      (amp-func (gaussian-beam beam-sigma k src-pt)))))

(run-until 50 (in-volume (volume (center 0 0 0) (size 20 8 0))
                         (at-end output-efield-z)))
