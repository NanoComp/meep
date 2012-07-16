; Example file illustrating an eigenmode source, generating a waveguide mode
; (requires recent MPB version to be installed before Meep is compiled)

(set! geometry-lattice (make lattice (size 16 8 no-size)))

; an asymmetrical dielectric waveguide:
(set! geometry (list
                (make block (center 0 0) (size infinity 1 infinity)
                      (material (make dielectric (epsilon 12))))
		(make block (center 0 0.3) (size infinity 0.1 infinity)
		      (material air))))

; create a transparent source that excites a right-going waveguide mode
(set! sources (list
               (make eigenmode-source
                 (src (make continuous-src (frequency 0.15)))
		 (size 0 6 0)
                 (center -5 0)
		 (component ALL-COMPONENTS)
		 (eig-parity TM)
		 )))
(set! pml-layers (list (make pml (thickness 1.0))))

(set-param! force-complex-fields? true) ; so we can get time-average flux

(set-param! resolution 10)

(run-until 200
           (at-beginning output-epsilon)
	   (at-end (output-png+h5 Ez "-a yarg -A $EPS -S3 -Zc dkbluered")))

(print "left-going flux = " ; (averaged over y region of width 1.8)
       (/ (flux-in-box X (volume (center -6 0) (size 1.8 6 0))) -1.8)
       "\n")

(print "right-going flux = " ; (averaged over y region of width 1.8)
       (/ (flux-in-box X (volume (center +6 0) (size 1.8 6 0))) +1.8)
       "\n")
