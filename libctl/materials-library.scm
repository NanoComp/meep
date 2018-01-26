; materials library for metallic optoelectronic materials
; from A.D. Rakic et al., Applied Optics, Vol. 37, No. 22, pp. 5271-83 (1998)

; conversion factor for eV to 1/um [=1/hc]
(define eV-um-scale (/ 1.23984193))

;------------------------------------------------------------------

(define Ag-plasma-frq (* 9.01 eV-um-scale))

(define Ag-f0 0.845)
(define Ag-frq0 1e-10)
(define Ag-gam0 (* 0.048 eV-um-scale))
(define Ag-sig0 (/ (* Ag-f0 (sqr Ag-plasma-frq)) (sqr Ag-frq0)))

(define Ag-f1 0.065)
(define Ag-frq1 (* 0.816 eV-um-scale)) ; 1.519 um
(define Ag-gam1 (* 3.886 eV-um-scale))
(define Ag-sig1 (/ (* Ag-f1 (sqr Ag-plasma-frq)) (sqr Ag-frq1)))

(define Ag-f2 0.124)
(define Ag-frq2 (* 4.481 eV-um-scale)) ; 0.273 um
(define Ag-gam2 (* 0.452 eV-um-scale))
(define Ag-sig2 (/ (* Ag-f2 (sqr Ag-plasma-frq)) (sqr Ag-frq2)))

(define Ag-f3 0.011)
(define Ag-frq3 (* 8.185 eV-um-scale)) ; 0.152 um
(define Ag-gam3 (* 0.065 eV-um-scale))
(define Ag-sig3 (/ (* Ag-f3 (sqr Ag-plasma-frq)) (sqr Ag-frq3)))

(define Ag-f4 0.840)
(define Ag-frq4 (* 9.083 eV-um-scale)) ; 0.137 um
(define Ag-gam4 (* 0.916 eV-um-scale))
(define Ag-sig4 (/ (* Ag-f4 (sqr Ag-plasma-frq)) (sqr Ag-frq4)))

(define Ag (make medium (epsilon 1.0)
  (E-susceptibilities
     (make drude-susceptibility
       (frequency Ag-frq0) (gamma Ag-gam0) (sigma Ag-sig0))
     (make lorentzian-susceptibility
       (frequency Ag-frq1) (gamma Ag-gam1) (sigma Ag-sig1))
     (make lorentzian-susceptibility
       (frequency Ag-frq2) (gamma Ag-gam2) (sigma Ag-sig2))
     (make lorentzian-susceptibility
       (frequency Ag-frq3) (gamma Ag-gam3) (sigma Ag-sig3))
     (make lorentzian-susceptibility
       (frequency Ag-frq4) (gamma Ag-gam4) (sigma Ag-sig4)))))

;------------------------------------------------------------------

(define Au-plasma-frq (* 9.03 eV-um-scale))

(define Au-f0 0.760)
(define Au-frq0 1e-10)
(define Au-gam0 (* 0.053 eV-um-scale))
(define Au-sig0 (/ (* Au-f0 (sqr Au-plasma-frq)) (sqr Au-frq0)))

(define Au-f1 0.024)
(define Au-frq1 (* 0.415 eV-um-scale)) ; 2.988 um
(define Au-gam1 (* 0.241 eV-um-scale))
(define Au-sig1 (/ (* Au-f1 (sqr Au-plasma-frq)) (sqr Au-frq1)))

(define Au-f2 0.010)
(define Au-frq2 (* 0.830 eV-um-scale)) ; 1.494 um
(define Au-gam2 (* 0.345 eV-um-scale))
(define Au-sig2 (/ (* Au-f2 (sqr Au-plasma-frq)) (sqr Au-frq2)))

(define Au-f3 0.071)
(define Au-frq3 (* 2.969 eV-um-scale)) ; 0.418 um
(define Au-gam3 (* 0.870 eV-um-scale))
(define Au-sig3 (/ (* Au-f3 (sqr Au-plasma-frq)) (sqr Au-frq3)))

(define Au-f4 0.601)
(define Au-frq4 (* 4.304 eV-um-scale)) ; 0.288 um
(define Au-gam4 (* 2.494 eV-um-scale))
(define Au-sig4 (/ (* Au-f4 (sqr Au-plasma-frq)) (sqr Au-frq4)))

(define Au (make medium (epsilon 1.0)
  (E-susceptibilities
     (make drude-susceptibility
       (frequency Au-frq0) (gamma Au-gam0) (sigma Au-sig0))
     (make lorentzian-susceptibility
       (frequency Au-frq1) (gamma Au-gam1) (sigma Au-sig1))
     (make lorentzian-susceptibility
       (frequency Au-frq2) (gamma Au-gam2) (sigma Au-sig2))
     (make lorentzian-susceptibility
       (frequency Au-frq3) (gamma Au-gam3) (sigma Au-sig3))
     (make lorentzian-susceptibility
       (frequency Au-frq4) (gamma Au-gam4) (sigma Au-sig4)))))

;------------------------------------------------------------------

(define Cu-plasma-frq (* 10.83 eV-um-scale))

(define Cu-f0 0.575)
(define Cu-frq0 1e-10)
(define Cu-gam0 (* 0.030 eV-um-scale))
(define Cu-sig0 (/ (* Cu-f0 (sqr Cu-plasma-frq)) (sqr Cu-frq0)))

(define Cu-f1 0.061)
(define Cu-frq1 (* 0.291 eV-um-scale)) ; 4.261 um
(define Cu-gam1 (* 0.378 eV-um-scale))
(define Cu-sig1 (/ (* Cu-f1 (sqr Cu-plasma-frq)) (sqr Cu-frq1)))

(define Cu-f2 0.104)
(define Cu-frq2 (* 2.957 eV-um-scale)) ; 0.419 um
(define Cu-gam2 (* 1.056 eV-um-scale))
(define Cu-sig2 (/ (* Cu-f2 (sqr Cu-plasma-frq)) (sqr Cu-frq2)))

(define Cu-f3 0.723)
(define Cu-frq3 (* 5.300 eV-um-scale)) ; 0.234 um
(define Cu-gam3 (* 3.213 eV-um-scale))
(define Cu-sig3 (/ (* Cu-f3 (sqr Cu-plasma-frq)) (sqr Cu-frq3)))

(define Cu-f4 0.638)
(define Cu-frq4 (* 11.18 eV-um-scale)) ; 0.111 um
(define Cu-gam4 (* 4.305 eV-um-scale))
(define Cu-sig4 (/ (* Cu-f4 (sqr Cu-plasma-frq)) (sqr Cu-frq4)))

(define Cu (make medium (epsilon 1.0)
  (E-susceptibilities
     (make drude-susceptibility
       (frequency Cu-frq0) (gamma Cu-gam0) (sigma Cu-sig0))
     (make lorentzian-susceptibility
       (frequency Cu-frq1) (gamma Cu-gam1) (sigma Cu-sig1))
     (make lorentzian-susceptibility
       (frequency Cu-frq2) (gamma Cu-gam2) (sigma Cu-sig2))
     (make lorentzian-susceptibility
       (frequency Cu-frq3) (gamma Cu-gam3) (sigma Cu-sig3))
     (make lorentzian-susceptibility
       (frequency Cu-frq4) (gamma Cu-gam4) (sigma Cu-sig4)))))

;------------------------------------------------------------------

(define Al-plasma-frq (* 14.98 eV-um-scale))

(define Al-f0 0.523)
(define Al-frq0 1e-10)
(define Al-gam0 (* 0.047 eV-um-scale))
(define Al-sig0 (/ (* Al-f0 (sqr Al-plasma-frq)) (sqr Al-frq0)))

(define Al-f1 0.227)
(define Al-frq1 (* 0.162 eV-um-scale)) ; 7.654 um
(define Al-gam1 (* 0.333 eV-um-scale))
(define Al-sig1 (/ (* Al-f1 (sqr Al-plasma-frq)) (sqr Al-frq1)))

(define Al-f2 0.050)
(define Al-frq2 (* 1.544 eV-um-scale)) ; 0.803 um
(define Al-gam2 (* 0.312 eV-um-scale))
(define Al-sig2 (/ (* Al-f2 (sqr Al-plasma-frq)) (sqr Al-frq2)))

(define Al-f3 0.166)
(define Al-frq3 (* 1.808 eV-um-scale)) ; 0.686 um
(define Al-gam3 (* 1.351 eV-um-scale))
(define Al-sig3 (/ (* Al-f3 (sqr Al-plasma-frq)) (sqr Al-frq3)))

(define Al-f4 0.030)
(define Al-frq4 (* 3.473 eV-um-scale)) ; 0.357 um
(define Al-gam4 (* 3.382 eV-um-scale))
(define Al-sig4 (/ (* Al-f4 (sqr Al-plasma-frq)) (sqr Al-frq4)))

(define Al (make medium (epsilon 1.0)
  (E-susceptibilities
     (make drude-susceptibility
       (frequency Al-frq0) (gamma Al-gam0) (sigma Al-sig0))
     (make lorentzian-susceptibility
       (frequency Al-frq1) (gamma Al-gam1) (sigma Al-sig1))
     (make lorentzian-susceptibility
       (frequency Al-frq2) (gamma Al-gam2) (sigma Al-sig2))
     (make lorentzian-susceptibility
       (frequency Al-frq3) (gamma Al-gam3) (sigma Al-sig3))
     (make lorentzian-susceptibility
       (frequency Al-frq4) (gamma Al-gam4) (sigma Al-sig4)))))

;------------------------------------------------------------------

(define Be-plasma-frq (* 18.51 eV-um-scale))

(define Be-f0 0.084)
(define Be-frq0 1e-10)
(define Be-gam0 (* 0.035 eV-um-scale))
(define Be-sig0 (/ (* Be-f0 (sqr Be-plasma-frq)) (sqr Be-frq0)))

(define Be-f1 0.031)
(define Be-frq1 (* 0.100 eV-um-scale)) ; 12.398 um
(define Be-gam1 (* 1.664 eV-um-scale))
(define Be-sig1 (/ (* Be-f1 (sqr Be-plasma-frq)) (sqr Be-frq1)))

(define Be-f2 0.140)
(define Be-frq2 (* 1.032 eV-um-scale)) ; 1.201 um
(define Be-gam2 (* 3.395 eV-um-scale))
(define Be-sig2 (/ (* Be-f2 (sqr Be-plasma-frq)) (sqr Be-frq2)))

(define Be-f3 0.530)
(define Be-frq3 (* 3.183 eV-um-scale)) ; 0.390 um
(define Be-gam3 (* 4.454 eV-um-scale))
(define Be-sig3 (/ (* Be-f3 (sqr Be-plasma-frq)) (sqr Be-frq3)))

(define Be-f4 0.130)
(define Be-frq4 (* 4.604 eV-um-scale)) ; 0.269 um
(define Be-gam4 (* 1.802 eV-um-scale))
(define Be-sig4 (/ (* Be-f4 (sqr Be-plasma-frq)) (sqr Be-frq4)))

(define Be (make medium (epsilon 1.0)
  (E-susceptibilities
     (make drude-susceptibility
       (frequency Be-frq0) (gamma Be-gam0) (sigma Be-sig0))
     (make lorentzian-susceptibility
       (frequency Be-frq1) (gamma Be-gam1) (sigma Be-sig1))
     (make lorentzian-susceptibility
       (frequency Be-frq2) (gamma Be-gam2) (sigma Be-sig2))
     (make lorentzian-susceptibility
       (frequency Be-frq3) (gamma Be-gam3) (sigma Be-sig3))
     (make lorentzian-susceptibility
       (frequency Be-frq4) (gamma Be-gam4) (sigma Be-sig4)))))

;------------------------------------------------------------------

(define Cr-plasma-frq (* 10.75 eV-um-scale))

(define Cr-f0 0.168)
(define Cr-frq0 1e-10)
(define Cr-gam0 (* 0.047 eV-um-scale))
(define Cr-sig0 (/ (* Cr-f0 (sqr Cr-plasma-frq)) (sqr Cr-frq0)))

(define Cr-f1 0.151)
(define Cr-frq1 (* 0.121 eV-um-scale)) ; 10.247 um
(define Cr-gam1 (* 3.175 eV-um-scale))
(define Cr-sig1 (/ (* Cr-f1 (sqr Cr-plasma-frq)) (sqr Cr-frq1)))

(define Cr-f2 0.150)
(define Cr-frq2 (* 0.543 eV-um-scale)) ; 2.283 um
(define Cr-gam2 (* 1.305 eV-um-scale))
(define Cr-sig2 (/ (* Cr-f2 (sqr Cr-plasma-frq)) (sqr Cr-frq2)))

(define Cr-f3 1.149)
(define Cr-frq3 (* 1.970 eV-um-scale)) ; 0.629 um
(define Cr-gam3 (* 2.676 eV-um-scale))
(define Cr-sig3 (/ (* Cr-f3 (sqr Cr-plasma-frq)) (sqr Cr-frq3)))

(define Cr-f4 0.825)
(define Cr-frq4 (* 8.775 eV-um-scale)) ; 0.141 um
(define Cr-gam4 (* 1.335 eV-um-scale))
(define Cr-sig4 (/ (* Cr-f4 (sqr Cr-plasma-frq)) (sqr Cr-frq4)))

(define Cr (make medium (epsilon 1.0)
  (E-susceptibilities
     (make drude-susceptibility
       (frequency Cr-frq0) (gamma Cr-gam0) (sigma Cr-sig0))
     (make lorentzian-susceptibility
       (frequency Cr-frq1) (gamma Cr-gam1) (sigma Cr-sig1))
     (make lorentzian-susceptibility
       (frequency Cr-frq2) (gamma Cr-gam2) (sigma Cr-sig2))
     (make lorentzian-susceptibility
       (frequency Cr-frq3) (gamma Cr-gam3) (sigma Cr-sig3))
     (make lorentzian-susceptibility
       (frequency Cr-frq4) (gamma Cr-gam4) (sigma Cr-sig4)))))

;------------------------------------------------------------------

(define Ni-plasma-frq (* 15.92 eV-um-scale))

(define Ni-f0 0.096)
(define Ni-frq0 1e-10)
(define Ni-gam0 (* 0.048 eV-um-scale))
(define Ni-sig0 (/ (* Ni-f0 (sqr Ni-plasma-frq)) (sqr Ni-frq0)))

(define Ni-f1 0.100)
(define Ni-frq1 (* 0.174 eV-um-scale)) ; 7.126 um
(define Ni-gam1 (* 4.511 eV-um-scale))
(define Ni-sig1 (/ (* Ni-f1 (sqr Ni-plasma-frq)) (sqr Ni-frq1)))

(define Ni-f2 0.135)
(define Ni-frq2 (* 0.582 eV-um-scale)) ; 2.130 um
(define Ni-gam2 (* 1.334 eV-um-scale))
(define Ni-sig2 (/ (* Ni-f2 (sqr Ni-plasma-frq)) (sqr Ni-frq2)))

(define Ni-f3 0.106)
(define Ni-frq3 (* 1.597 eV-um-scale)) ; 0.776 um
(define Ni-gam3 (* 2.178 eV-um-scale))
(define Ni-sig3 (/ (* Ni-f3 (sqr Ni-plasma-frq)) (sqr Ni-frq3)))

(define Ni-f4 0.729)
(define Ni-frq4 (* 6.089 eV-um-scale)) ; 0.204 um
(define Ni-gam4 (* 6.292 eV-um-scale))
(define Ni-sig4 (/ (* Ni-f4 (sqr Ni-plasma-frq)) (sqr Ni-frq4)))

(define Ni (make medium (epsilon 1.0)
  (E-susceptibilities
     (make drude-susceptibility
       (frequency Ni-frq0) (gamma Ni-gam0) (sigma Ni-sig0))
     (make lorentzian-susceptibility
       (frequency Ni-frq1) (gamma Ni-gam1) (sigma Ni-sig1))
     (make lorentzian-susceptibility
       (frequency Ni-frq2) (gamma Ni-gam2) (sigma Ni-sig2))
     (make lorentzian-susceptibility
       (frequency Ni-frq3) (gamma Ni-gam3) (sigma Ni-sig3))
     (make lorentzian-susceptibility
       (frequency Ni-frq4) (gamma Ni-gam4) (sigma Ni-sig4)))))

;------------------------------------------------------------------

(define Pd-plasma-frq (* 9.72 eV-um-scale))

(define Pd-f0 0.330)
(define Pd-frq0 1e-10)
(define Pd-gam0 (* 0.008 eV-um-scale))
(define Pd-sig0 (/ (* Pd-f0 (sqr Pd-plasma-frq)) (sqr Pd-frq0)))

(define Pd-f1 0.649)
(define Pd-frq1 (* 0.336 eV-um-scale)) ; 3.690 um
(define Pd-gam1 (* 2.950 eV-um-scale))
(define Pd-sig1 (/ (* Pd-f1 (sqr Pd-plasma-frq)) (sqr Pd-frq1)))

(define Pd-f2 0.121)
(define Pd-frq2 (* 0.501 eV-um-scale)) ; 2.475 um
(define Pd-gam2 (* 0.555 eV-um-scale))
(define Pd-sig2 (/ (* Pd-f2 (sqr Pd-plasma-frq)) (sqr Pd-frq2)))

(define Pd-f3 0.638)
(define Pd-frq3 (* 1.659 eV-um-scale)) ; 0.747 um
(define Pd-gam3 (* 4.621 eV-um-scale))
(define Pd-sig3 (/ (* Pd-f3 (sqr Pd-plasma-frq)) (sqr Pd-frq3)))

(define Pd-f4 0.453)
(define Pd-frq4 (* 5.715 eV-um-scale)) ; 0.217 um
(define Pd-gam4 (* 3.236 eV-um-scale))
(define Pd-sig4 (/ (* Pd-f4 (sqr Pd-plasma-frq)) (sqr Pd-frq4)))

(define Pd (make medium (epsilon 1.0)
  (E-susceptibilities
     (make drude-susceptibility
       (frequency Pd-frq0) (gamma Pd-gam0) (sigma Pd-sig0))
     (make lorentzian-susceptibility
       (frequency Pd-frq1) (gamma Pd-gam1) (sigma Pd-sig1))
     (make lorentzian-susceptibility
       (frequency Pd-frq2) (gamma Pd-gam2) (sigma Pd-sig2))
     (make lorentzian-susceptibility
       (frequency Pd-frq3) (gamma Pd-gam3) (sigma Pd-sig3))
     (make lorentzian-susceptibility
       (frequency Pd-frq4) (gamma Pd-gam4) (sigma Pd-sig4)))))

;------------------------------------------------------------------

(define Pt-plasma-frq (* 9.59 eV-um-scale))

(define Pt-f0 0.333)
(define Pt-frq0 1e-10)
(define Pt-gam0 (* 0.080 eV-um-scale))
(define Pt-sig0 (/ (* Pt-f0 (sqr Pt-plasma-frq)) (sqr Pt-frq0)))

(define Pt-f1 0.191)
(define Pt-frq1 (* 0.780 eV-um-scale)) ; 1.590 um
(define Pt-gam1 (* 0.517 eV-um-scale))
(define Pt-sig1 (/ (* Pt-f1 (sqr Pt-plasma-frq)) (sqr Pt-frq1)))

(define Pt-f2 0.659)
(define Pt-frq2 (* 1.314 eV-um-scale)) ; 0.944 um
(define Pt-gam2 (* 1.838 eV-um-scale))
(define Pt-sig2 (/ (* Pt-f2 (sqr Pt-plasma-frq)) (sqr Pt-frq2)))

(define Pt-f3 0.547)
(define Pt-frq3 (* 3.141 eV-um-scale)) ; 0.395 um
(define Pt-gam3 (* 3.668 eV-um-scale))
(define Pt-sig3 (/ (* Pt-f3 (sqr Pt-plasma-frq)) (sqr Pt-frq3)))

(define Pt-f4 3.576)
(define Pt-frq4 (* 9.249 eV-um-scale)) ; 0.134 um
(define Pt-gam4 (* 8.517 eV-um-scale))
(define Pt-sig4 (/ (* Pt-f4 (sqr Pt-plasma-frq)) (sqr Pt-frq4)))

(define Pt (make medium (epsilon 1.0)
  (E-susceptibilities
     (make drude-susceptibility
       (frequency Pt-frq0) (gamma Pt-gam0) (sigma Pt-sig0))
     (make lorentzian-susceptibility
       (frequency Pt-frq1) (gamma Pt-gam1) (sigma Pt-sig1))
     (make lorentzian-susceptibility
       (frequency Pt-frq2) (gamma Pt-gam2) (sigma Pt-sig2))
     (make lorentzian-susceptibility
       (frequency Pt-frq3) (gamma Pt-gam3) (sigma Pt-sig3))
     (make lorentzian-susceptibility
       (frequency Pt-frq4) (gamma Pt-gam4) (sigma Pt-sig4)))))

;------------------------------------------------------------------

(define Ti-plasma-frq (* 7.29 eV-um-scale))

(define Ti-f0 0.148)
(define Ti-frq0 1e-10)
(define Ti-gam0 (* 0.082 eV-um-scale))
(define Ti-sig0 (/ (* Ti-f0 (sqr Ti-plasma-frq)) (sqr Ti-frq0)))

(define Ti-f1 0.899)
(define Ti-frq1 (* 0.777 eV-um-scale)) ; 1.596 um
(define Ti-gam1 (* 2.276 eV-um-scale))
(define Ti-sig1 (/ (* Ti-f1 (sqr Ti-plasma-frq)) (sqr Ti-frq1)))

(define Ti-f2 0.393)
(define Ti-frq2 (* 1.545 eV-um-scale)) ; 0.802 um
(define Ti-gam2 (* 2.518 eV-um-scale))
(define Ti-sig2 (/ (* Ti-f2 (sqr Ti-plasma-frq)) (sqr Ti-frq2)))

(define Ti-f3 0.187)
(define Ti-frq3 (* 2.509 eV-um-scale)) ; 0.494 um
(define Ti-gam3 (* 1.663 eV-um-scale))
(define Ti-sig3 (/ (* Ti-f3 (sqr Ti-plasma-frq)) (sqr Ti-frq3)))

(define Ti-f4 0.001)
(define Ti-frq4 (* 19.43 eV-um-scale)) ; 0.064 um
(define Ti-gam4 (* 1.762 eV-um-scale))
(define Ti-sig4 (/ (* Ti-f4 (sqr Ti-plasma-frq)) (sqr Ti-frq4)))

(define Ti (make medium (epsilon 1.0)
  (E-susceptibilities
     (make drude-susceptibility
       (frequency Ti-frq0) (gamma Ti-gam0) (sigma Ti-sig0))
     (make lorentzian-susceptibility
       (frequency Ti-frq1) (gamma Ti-gam1) (sigma Ti-sig1))
     (make lorentzian-susceptibility
       (frequency Ti-frq2) (gamma Ti-gam2) (sigma Ti-sig2))
     (make lorentzian-susceptibility
       (frequency Ti-frq3) (gamma Ti-gam3) (sigma Ti-sig3))
     (make lorentzian-susceptibility
       (frequency Ti-frq4) (gamma Ti-gam4) (sigma Ti-sig4)))))

;------------------------------------------------------------------

(define W-plasma-frq (* 13.22 eV-um-scale))

(define W-f0 0.206)
(define W-frq0 1e-10)
(define W-gam0 (* 0.064 eV-um-scale))
(define W-sig0 (/ (* W-f0 (sqr W-plasma-frq)) (sqr W-frq0)))

(define W-f1 0.054)
(define W-frq1 (* 1.004 eV-um-scale)) ; 1.235 um
(define W-gam1 (* 0.530 eV-um-scale))
(define W-sig1 (/ (* W-f1 (sqr W-plasma-frq)) (sqr W-frq1)))

(define W-f2 0.166)
(define W-frq2 (* 1.917 eV-um-scale)) ; 0.647 um
(define W-gam2 (* 1.281 eV-um-scale))
(define W-sig2 (/ (* W-f2 (sqr W-plasma-frq)) (sqr W-frq2)))

(define W-f3 0.706)
(define W-frq3 (* 3.580 eV-um-scale)) ; 0.346 um
(define W-gam3 (* 3.332 eV-um-scale))
(define W-sig3 (/ (* W-f3 (sqr W-plasma-frq)) (sqr W-frq3)))

(define W-f4 2.590)
(define W-frq4 (* 7.498 eV-um-scale)) ; 0.165 um
(define W-gam4 (* 5.836 eV-um-scale))
(define W-sig4 (/ (* W-f4 (sqr W-plasma-frq)) (sqr W-frq4)))

(define W (make medium (epsilon 1.0)
  (E-susceptibilities
     (make drude-susceptibility
       (frequency W-frq0) (gamma W-gam0) (sigma W-sig0))
     (make lorentzian-susceptibility
       (frequency W-frq1) (gamma W-gam1) (sigma W-sig1))
     (make lorentzian-susceptibility
       (frequency W-frq2) (gamma W-gam2) (sigma W-sig2))
     (make lorentzian-susceptibility
       (frequency W-frq3) (gamma W-gam3) (sigma W-sig3))
     (make lorentzian-susceptibility
       (frequency W-frq4) (gamma W-gam4) (sigma W-sig4)))))
