;given m1 m2, make a composit index n; the inverse is casimir-source-info below
; m1 = r-c, m2 = c  => r = m1+m2;  n - s = m2 => n = m2 + s; s = (sum_(k=1)^r k)
; => n = m2 + (sum_(k=1)^(m1+m2) k) = m2 + 1/2 (m1+m2) * (m1+m2+1)
(define (make-casimir-src-index m1 m2) (+ m2 (* (/ 2) (+ m1 m2) (+ m1 m2 1))))

; return a list (source-vol mx my mz)
; given the volume integration-vol and n, pick out the appropriate side and mx my mz to use
; sides are ordered by decreasing weight
; weights are w(side, m) = area(side)/total area * 1/(m+1)^4, a rough estimate of the
; contribution to the stress tensor from that side and that m
(define (casimir-source-info integration-vol n)
  (define (modround x n) (modulo (inexact->exact (round x)) n))
  (define (get-src-index n) ;given n, extract out the two values of m for 3-d
    (let* ((s 0) ;sum of diagonals
	 (r 0) ;row intersection
	 (c 0));column intersection
      (while (< (+ s r) n)
	     (set! r (+ r 1))
	     (set! s (+ s r)))
      (set! c (- n s))
      (list (- r c) c)))

  (let* ((min-corner (meep-volume-get-min-corner integration-vol))
	 (max-corner (meep-volume-get-max-corner integration-vol))
	 (size-vec (vector3- max-corner min-corner))
	 (center-vec (vector3+ (vector3-scale 0.5 size-vec) min-corner))
	 (sx (vector3-x size-vec))
	 (sy (vector3-y size-vec))
	 (sz (vector3-z size-vec))
	 (xshift (vector3 (/ sx 2)    0       0))
	 (yshift (vector3    0    (/ sy 2)    0))
	 (zshift (vector3    0        0   (/ sz 2))))
    (if (and (> sy 1e-15) (> sz 1e-15))  ;3d cartesian: n = 6*f(m1,m2) + s
	(let* ((s (modround n 6))
	       (nr (/ (- n s) 6))
	       (ms (get-src-index nr)) ;get (m1 m2)
	       (m1 (first ms))
	       (m2 (second ms))
	       (x-const (vector3 0 sy sz))
	       (y-const (vector3 sx 0 sz))
	       (z-const (vector3 sx sy 0))
	       (center-list
		(list (vector3- center-vec xshift) (vector3+ center-vec xshift)
		      (vector3- center-vec yshift) (vector3+ center-vec yshift)
		      (vector3- center-vec zshift) (vector3+ center-vec zshift)))
	       (m-list
		(list (vector3 0 m1 m2) (vector3 0 m1 m2) (vector3 m1 0 m2)
		      (vector3 m1 0 m2) (vector3 m1 m2 0) (vector3 m1 m2 0)))
	       (orientation-list (list -1 1 -1 1 -1 1))
	       (size-list
		(list x-const x-const y-const y-const z-const z-const))
	       (surface-vol
		(volume (center (list-ref center-list s)) (size (list-ref size-list s))))
	       (surface-m (list-ref m-list s)))
	  (print "Computing in 3d\n")
	  (list surface-vol
		(vector3-x surface-m) (vector3-y surface-m) (vector3-z surface-m)
		(list-ref orientation-list s) 1))
	(if (= dimensions -2) ;cylindricals - must make sure that the volume has only r >= 0
	    (let* ((3-sides? (if (<= (vector3-x min-corner) 0) true false)) ;volume passes through the origin
		   (s (if 3-sides? (modround n 3) (modround n 4)))
		   (nr (if 3-sides? (/ (- n s) 3) (/ (- n s) 4))) ;reduced index
		   (ms (get-src-index nr)) ;extract out both m-phi and m-dct
		   (m-phi (first ms))
		   (m-dct (second ms))
		   (DR (if 3-sides? (/ 1 resolution) 0)) ;cannot include r = 0!!
		   (sr (- (vector3-x max-corner) (+ (max 0 (vector3-x min-corner)) DR)))
		   (r-cen (+ (* 0.5 sr) (max 0 (vector3-x min-corner)) DR))
		   (new-center-vec (vector3 r-cen 0 (vector3-z center-vec)))
		   (r-shift (vector3 (/ sr 2) 0 0))
		   (z-shift (vector3 0 0 (/ sz 2)))
		   (r-const-size (vector3 0 0 sz))
		   (z-const-size (vector3 sr 0 0))
		   (center-list ;if 3-sides? = true, only first 3 list elements are used
		    (list (vector3- new-center-vec z-shift) (vector3+ new-center-vec z-shift)
			  (vector3+ new-center-vec r-shift) (vector3- new-center-vec r-shift)))
		   (m-list
		    (list (vector3 m-dct m-phi 0) (vector3 m-dct m-phi 0)
			  (vector3 0 m-phi m-dct) (vector3 0 m-phi m-dct)))
		   (orientation-list (list -1 1 1 -1))
		   (size-list
		    (list z-const-size z-const-size r-const-size r-const-size))
		   (surface-vol
		    (volume (center (list-ref center-list s)) (size (list-ref size-list s))))
		   (surface-m (list-ref m-list s)))
	      (print "Computing in Cylindrical coordinates: m-phi = "m-phi", m-dct = "m-dct", 3-sides? = "3-sides?"\n")
	      (list surface-vol
		    (vector3-x surface-m) (vector3-y surface-m) (vector3-z surface-m)
		    (list-ref orientation-list s) (if (or (= s 0) (= s 1)) 1 0)))
	    (let* ((s (modround n 4)) ;2d or quasi-3d cartesian: n = 4m + s, no ambiguity in m
		   (m (/ (- n s) 4))
		   (x-const-size (vector3 0 sy )) ;sz may be non-zero for quasi-3d systems
		   (y-const-size (vector3 sx 0 ))
		   (center-list
		    (list (vector3- center-vec xshift) (vector3+ center-vec xshift)
			  (vector3- center-vec yshift) (vector3+ center-vec yshift)))
		   (m-list
		    (list (vector3 0 m 0) (vector3 0 m 0) (vector3 m 0 0) (vector3 m 0 0)))
		   (orientation-list (list -1 1 -1 1))
		   (size-list
		    (list x-const-size x-const-size y-const-size y-const-size))
		   (surface-vol
		    (volume (center (list-ref center-list s)) (size (list-ref size-list s))))
		   (surface-m (list-ref m-list s)))
	      (print "Casimir.scm: working in 2 dimensions\n")
	      (print "  Surface center: "(list-ref center-list s)"\n")
	      (print "  Surface size: "(list-ref size-list s)"\n")
	      (list surface-vol
		    (vector3-x surface-m) (vector3-y surface-m) (vector3-z surface-m)
		    (list-ref orientation-list s) 1))))))

;compute the casimir force for a single n and single polarization
;n contains both the side number and the harmonic expansion index
(define (casimir-force-contrib force-direction integration-vol n Sigma T source-component gt . step-funcs)
  (define (cos-func X mx my mz source-vol)
    (let*
	((min-corner (meep-volume-get-min-corner source-vol))
	 (max-corner (meep-volume-get-max-corner source-vol))
	 (size-vec (vector3- max-corner min-corner))
	 (X-start (vector3+ X (vector3-scale 0.5 size-vec)))
	 (sx (vector3-x size-vec))
	 (sy (vector3-y size-vec))
	 (sz (vector3-z size-vec))
	 (x (vector3-x X-start))
	 (y (vector3-y X-start))
	 (z (vector3-z X-start))
	 (kx (if (> sx 1e-15) (/ (* mx pi) sx) 0))
	 (ky (if (> sy 1e-15) (/ (* my pi) sy) 0))
	 (kz (if (> sz 1e-15) (/ (* mz pi) sz) 0))
	 (Nx (if (> sx 1e-15) (/ (if (= mx 0) 1 2) sx) 1))
	 (Ny (if (> sy 1e-15) (/ (if (= my 0) 1 2) sy) 1))
	 (Nz (if (> sz 1e-15) (/ (if (= mz 0) 1 2) sz) 1)))
      (* (sqrt (* Nx Ny Nz)) (cos (* kx x)) (cos (* ky y)) (cos (* kz z)))))

  (let* ((ft (meep-type source-component))
	 (source-info (casimir-source-info integration-vol n))
	 (source-vol (first source-info))
	 (mx (second source-info))
	 (my (third source-info)) ;m-phi in cylindrical coordinates
	 (mz (fourth source-info))
	 (source-orientation (fifth source-info))
	 (dt (/ Courant resolution)))
    (if (= ft E-stuff)
	(begin
	  (set! global-D-conductivity Sigma)
	  (set! global-B-conductivity 0))
	(begin
	  (set! global-B-conductivity Sigma)
	  (set! global-D-conductivity 0)))
    (if (eq? dimensions -2)
	(begin (print "Cylindricals: m = "my" and (nr nz) = ("mx", "mz")\n")
	       (print "  surface center = "(meep-volume-center source-vol)"\n")
	       (print "  source size = "(vector3-
					  (meep-volume-get-max-corner source-vol)
					  (meep-volume-get-min-corner source-vol)))
	       (set! m my))) ;set exp(i m phi) field dependence
    (set! sources
	  (list (make source
		  (src (make custom-src ; delta function pulse
			 (src-func (lambda (t) (/ 1 dt)))
			 (start-time (* -0.25 dt))
			 (end-time (* 0.75 dt))
			 (is-integrated? false)))
		  (center (meep-volume-center source-vol))
		  (size (vector3- (meep-volume-get-max-corner source-vol)
				  (meep-volume-get-min-corner source-vol)))
		  (component source-component)
		  (amp-func (lambda (p) (cos-func p mx my mz source-vol))))))
    (reset-meep)
    (init-fields)
    (let* ((counter 0)
	   (force-integral 0))
      (define (integrate-function)
	(let* ((f-temp (meep-fields-casimir-stress-dct-integral
			fields
			force-direction
			(meep-component-direction source-component)
			mx (if (eq? dimensions -2) 0 my) mz
			ft source-vol)))
	  (set! force-integral
		(+ force-integral
		   (imag-part (* (list-ref gt counter) dt
				 source-orientation
				 (if (eq? dimensions -2)
				     (* (if (eq? my 0) 1 2)
					(real-part f-temp))
				     f-temp)))))
	  (set! counter (+ counter 1))))
      (apply run-until (cons (- T 1) (cons integrate-function step-funcs)))
      force-integral)))

;%%%%%%%%%%%%%%%%%%%%% BLOCH PBCS %%%%%%%%%%%%%%%%%%%%%%
;here the source is specified in
;the form exp(i k x), k = pi/L (m + k_red),
;m = (mx,my,mz) are integers (reciprocal lattice vectors
;k_red = (kx,ky,kz) is in the 1st BZ, m an integer
;source-vol is assumed to occupy one entire plane intersecting
;the computational cell, so we don't need to extract out
;its information - there is only one side to it

;pass the vector (mx my mz) and (kx ky kz) ready-made, since
;this surface consists of only one face
(define (casimir-force-contrib-bloch force-direction source-vol k-vec Sigma T source-component gt . step-funcs)
  ;sources of the form exp(i g x); surface integration in
  ;casimir.cpp integrates against exp(-i g x)
  (define (casimir-bloch-func X gx gy gz source-vol)
    (let*
	((min-corner (meep-volume-get-min-corner source-vol))
	 (max-corner (meep-volume-get-max-corner source-vol))
	 (size-vec (vector3- max-corner min-corner)) ;crossection of computational cell
	 (sx (vector3-x size-vec))
	 (sy (vector3-y size-vec))
	 (sz (vector3-z size-vec))
	 (x (vector3-x X))
	 (y (vector3-y X))
	 (z (vector3-z X))
	 (Kx (if (> sx 1e-15) (/ (* gx pi) 1) 0)) ;phase winding is independent of unit cell size
	 (Ky (if (> sy 1e-15) (/ (* gy pi) 1) 0))
	 (Kz (if (> sz 1e-15) (/ (* gz pi) 1) 0))
	 (Nx (if (> sx 1e-15) (/ sx) 1))
	 (Ny (if (> sy 1e-15) (/ sy) 1))
	 (Nz (if (> sz 1e-15) (/ sz) 1)))
      (* (sqrt (* Nx Ny Nz)) (exp (* (sqrt -1) (+ (* Kx x) (* Ky y) (* Kz z)))))))

  (let* ((ft (meep-type source-component))
	 (dt (/ Courant resolution))
	 ;Bloch phases - exp( i * (2*pi*m + pi*k) x/L)
	 (gx (vector3-x k-vec))
	 (gy (vector3-y k-vec))
	 (gz (vector3-z k-vec)))
    (set! force-complex-fields? true)
    (set! k-point (vector3-scale 0.5 k-vec))
    (if (= ft E-stuff)
	(begin
	  (set! global-D-conductivity Sigma)
	  (set! global-B-conductivity 0))
	(begin
	  (set! global-D-conductivity 0)
	  (set! global-B-conductivity Sigma)))
    (set! sources
	  (list (make source
		  (src (make custom-src
			 (src-func (lambda (t) (/ 1 dt)))
			 (start-time (* -0.25 dt))
			 (end-time (* 0.75 dt))
			 (is-integrated? false)))
		  (center (meep-volume-center source-vol))
		  (size (vector3- (meep-volume-get-max-corner source-vol)
				  (meep-volume-get-min-corner source-vol)))
		  (component source-component)
		  (amp-func (lambda (p) (casimir-bloch-func p gx gy gz source-vol))))))
    (reset-meep)
    (init-fields)
    (let* ((counter 0)
	   (force-integral 0))
      (define (integrate-function)
	(set! force-integral
	      (+ force-integral
		 (imag-part (* (list-ref gt counter) dt
			       (meep-fields-casimir-stress-dct-integral
				fields
				force-direction
				(meep-component-direction source-component)
				gx gy gz
				ft source-vol true)))))
	(set! counter (+ counter 1)))
      (apply run-until (cons (- T 1) (cons integrate-function step-funcs)))
      force-integral)))
