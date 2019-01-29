;; From the Meep tutorial: plotting permittivity and fields of a straight waveguide

(set! geometry-lattice (make lattice (size 16 8 no-size)))

(set! geometry (list
                (make block (center 0 0) (size infinity 1 infinity)
                      (material (make medium (epsilon 12))))))

(set! sources (list
               (make source
                 (src (make continuous-src (frequency 0.15)))
                 (component Ez)
                 (center -7 0))))

(set! pml-layers (list (make pml (thickness 1.0))))

(set! resolution 10)

(run-until 200
           (at-beginning output-epsilon)
           (at-end output-efield-z))
