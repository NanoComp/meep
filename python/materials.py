# -*- coding: utf-8 -*-
# Materials Library

import meep as mp

# default unit length is 1 um
um_scale = 1.0

# conversion factor for eV to 1/um [=1/hc]
eV_um_scale = um_scale/1.23984193

#------------------------------------------------------------------

# crystalline silicon (cSi) from A. Deinega et al., J. Optical Society of America A, Vol. 28, No. 5, pp. 770-77, 2011
# based on experimental data for intrinsic silicon at T=300K from M.A. Green and M. Keevers, Progress in Photovoltaics, Vol. 3, pp. 189-92, 1995
# wavelength range: 0.4 - 1.0 um

cSi_range = mp.FreqRange(min=um_scale, max=um_scale/0.4)

cSi_frq1 = 3.64/um_scale
cSi_gam1 = 0
cSi_sig1 = 8
cSi_frq2 = 2.76/um_scale
cSi_gam2 = 2*0.063/um_scale
cSi_sig2 = 2.85
cSi_frq3 = 1.73/um_scale
cSi_gam3 = 2*2.5/um_scale
cSi_sig3 = -0.107

cSi_susc = [mp.LorentzianSusceptibility(frequency=cSi_frq1, gamma=cSi_gam1, sigma=cSi_sig1),
            mp.LorentzianSusceptibility(frequency=cSi_frq2, gamma=cSi_gam2, sigma=cSi_sig2),
            mp.LorentzianSusceptibility(frequency=cSi_frq3, gamma=cSi_gam3, sigma=cSi_sig3)]

cSi = mp.Medium(epsilon=1.0, E_susceptibilities=cSi_susc, valid_freq_range=cSi_range)

#------------------------------------------------------------------

# amorphous silicon (a-Si) from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.21 - 0.83 um

aSi_range = mp.FreqRange(min=um_scale/0.83, max=um_scale/0.21)

aSi_frq1 = 1/(0.315481407124682*um_scale)
aSi_gam1 = 1/(0.645751005208333*um_scale)
aSi_sig1 = 14.571

aSi_susc = [mp.LorentzianSusceptibility(frequency=aSi_frq1, gamma=aSi_gam1, sigma=aSi_sig1)]

aSi = mp.Medium(epsilon=3.109, E_susceptibilities=aSi_susc, valid_freq_range=aSi_range)

#------------------------------------------------------------------

# hydrogenated amorphous silicon (a-Si:H) from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.21 - 0.83 um

aSi_H_range = mp.FreqRange(min=um_scale/0.83, max=um_scale/0.21)

aSi_H_frq1 = 1/(0.334189199460916*um_scale)
aSi_H_gam1 = 1/(0.579365387850467*um_scale)
aSi_H_sig1 = 12.31

aSi_H_susc = [mp.LorentzianSusceptibility(frequency=aSi_H_frq1, gamma=aSi_H_gam1, sigma=aSi_H_sig1)]

aSi_H = mp.Medium(epsilon=3.22, E_susceptibilities=aSi_H_susc, valid_freq_range=aSi_H_range)

#------------------------------------------------------------------

# indium tin oxide (ITO) from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.21 - 0.83 um

ITO_range = mp.FreqRange(min=um_scale/0.83, max=um_scale/0.21)

ITO_frq1 = 1/(0.182329695588235*um_scale)
ITO_gam1 = 1/(1.94637665620094*um_scale)
ITO_sig1 = 2.5

ITO_susc = [mp.LorentzianSusceptibility(frequency=ITO_frq1, gamma=ITO_gam1, sigma=ITO_sig1)]

ITO = mp.Medium(epsilon=1.0, E_susceptibilities=ITO_susc, valid_freq_range=ITO_range)

#------------------------------------------------------------------

# alumina (Al2O3) from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.21 - 2.07 um

Al2O3_range = mp.FreqRange(min=um_scale/2.07, max=um_scale/0.21)

Al2O3_frq1 = 1/(0.101476668030774*um_scale)
Al2O3_gam1 = 0
Al2O3_sig1 = 1.52

Al2O3_susc = [mp.LorentzianSusceptibility(frequency=Al2O3_frq1, gamma=Al2O3_gam1, sigma=Al2O3_sig1)]

Al2O3 = mp.Medium(epsilon=1.0, E_susceptibilities=Al2O3_susc, valid_freq_range=Al2O3_range)

#------------------------------------------------------------------

# aluminum nitride (AlN) from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.26 - 1.65 um

AlN_range = mp.FreqRange(min=um_scale/1.65, max=um_scale/0.26)

AlN_frq1 = 1/(0.139058089950651*um_scale)
AlN_gam1 = 0
AlN_sig1 = 3.306

AlN_susc = [mp.LorentzianSusceptibility(frequency=AlN_frq1, gamma=AlN_gam1, sigma=AlN_sig1)]

AlN = mp.Medium(epsilon=1.0, E_susceptibilities=AlN_susc, valid_freq_range=AlN_range)

#------------------------------------------------------------------

# aluminum arsenide (AlAs) from R.E. Fern and A. Onton, J. Applied Physics, Vol. 42, pp. 3499-500, 1971
# ref: https://refractiveindex.info/?shelf=main&book=AlAs&page=Fern
# wavelength range: 0.56 - 2.2 um

AlAs_range = mp.FreqRange(min=um_scale/2.2, max=um_scale/0.56)

AlAs_frq1 = 1/(0.2822*um_scale)
AlAs_gam1 = 0
AlAs_sig1 = 6.0840
AlAs_frq2 = 1/(27.62*um_scale)
AlAs_gam2 = 0
AlAs_sig2 = 1.900

AlAs_susc = [mp.LorentzianSusceptibility(frequency=AlAs_frq1, gamma=AlAs_gam1, sigma=AlAs_sig1),
             mp.LorentzianSusceptibility(frequency=AlAs_frq2, gamma=AlAs_gam2, sigma=AlAs_sig2)]

AlAs = mp.Medium(epsilon=2.0792, E_susceptibilities=AlAs_susc, valid_freq_range=AlAs_range)

#------------------------------------------------------------------

# borosilicate glass (BK7) from SCHOTT Zemax catalog 2017-01-20b
# ref: https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT
# wavelength range: 0.3 - 2.5 um

BK7_range = mp.FreqRange(min=um_scale/2.5, max=um_scale/0.3)

BK7_frq1 = 1/(0.07746417668832478*um_scale)
BK7_gam1 = 0
BK7_sig1 = 1.03961212
BK7_frq2 = 1/(0.14148467902921502*um_scale)
BK7_gam2 = 0
BK7_sig2 = 0.231792344
BK7_frq3 = 1/(10.176475470417055*um_scale)
BK7_gam3 = 0
BK7_sig3 = 1.01046945

BK7_susc = [mp.LorentzianSusceptibility(frequency=BK7_frq1, gamma=BK7_gam1, sigma=BK7_sig1),
            mp.LorentzianSusceptibility(frequency=BK7_frq2, gamma=BK7_gam2, sigma=BK7_sig2),
            mp.LorentzianSusceptibility(frequency=BK7_frq3, gamma=BK7_gam3, sigma=BK7_sig3)]

BK7 = mp.Medium(epsilon=1.0, E_susceptibilities=BK7_susc, valid_freq_range=BK7_range)

#------------------------------------------------------------------

# fused quartz (silica) from I.H. Malitson, J. Optical Society of America, Vol. 55, pp. 1205-9, 1965
# ref: https://refractiveindex.info/?shelf=glass&book=fused_silica&page=Malitson
# wavelength range: 0.21 - 6.7 um

fused_quartz_range = mp.FreqRange(min=um_scale/6.7, max=um_scale/0.21)

fused_quartz_frq1 = 1/(0.0684043*um_scale)
fused_quartz_gam1 = 0
fused_quartz_sig1 = 0.696166300
fused_quartz_frq2 = 1/(0.1162414*um_scale)
fused_quartz_gam2 = 0
fused_quartz_sig2 = 0.407942600
fused_quartz_frq3 = 1/(9.896161*um_scale)
fused_quartz_gam3 = 0
fused_quartz_sig3 = 0.897479400

fused_quartz_susc = [mp.LorentzianSusceptibility(frequency=fused_quartz_frq1, gamma=fused_quartz_gam1, sigma=fused_quartz_sig1),
                     mp.LorentzianSusceptibility(frequency=fused_quartz_frq2, gamma=fused_quartz_gam2, sigma=fused_quartz_sig2),
                     mp.LorentzianSusceptibility(frequency=fused_quartz_frq3, gamma=fused_quartz_gam3, sigma=fused_quartz_sig3)]

fused_quartz = mp.Medium(epsilon=1.0, E_susceptibilities=fused_quartz_susc, valid_freq_range=fused_quartz_range)

#------------------------------------------------------------------

# gallium arsenide (GaAs) from T. Skauli et al., J. Applied Physics, Vol. 94, pp. 6447-55, 2003
# ref: https://refractiveindex.info/?shelf=main&book=GaAs&page=Skauli
# wavelength range: 0.97 - 17 um

GaAs_range = mp.FreqRange(min=um_scale/17, max=um_scale/0.97)

GaAs_frq1 = 1/(0.4431307*um_scale)
GaAs_gam1 = 0
GaAs_sig1 = 5.466742
GaAs_frq2 = 1/(0.8746453*um_scale)
GaAs_gam2 = 0
GaAs_sig2 = 0.02429960
GaAs_frq3 = 1/(36.9166*um_scale)
GaAs_gam3 = 0
GaAs_sig3 = 1.957522

GaAs_susc = [mp.LorentzianSusceptibility(frequency=GaAs_frq1, gamma=GaAs_gam1, sigma=GaAs_sig1),
             mp.LorentzianSusceptibility(frequency=GaAs_frq2, gamma=GaAs_gam2, sigma=GaAs_sig2),
             mp.LorentzianSusceptibility(frequency=GaAs_frq3, gamma=GaAs_gam3, sigma=GaAs_sig3)]

GaAs = mp.Medium(epsilon=5.372514, E_susceptibilities=GaAs_susc, valid_freq_range=GaAs_range)

#------------------------------------------------------------------
# stoichiometric silicon nitride (Si3N4) from H. R. Philipp, Optical properties of silicon nitride, J. Electrochem. Soc. 120, 295-300, 1973
# fit from https://refractiveindex.info/?shelf=main&book=Si3N4&page=Philipp
# wavelength range: 0.207 - 1.24 um

Si3N4_VISNIR_range = mp.FreqRange(min=um_scale/1.24, max=um_scale/0.207)

Si3N4_VISNIR_frq1 = 1/(0.13967*um_scale)
Si3N4_VISNIR_gam1 = 0
Si3N4_VISNIR_sig1 = 2.8939

Si3N4_VISNIR_susc = [mp.LorentzianSusceptibility(frequency=Si3N4_VISNIR_frq1, gamma=Si3N4_VISNIR_gam1, sigma=Si3N4_VISNIR_sig1)]

Si3N4_VISNIR = mp.Medium(epsilon=1.0, E_susceptibilities=Si3N4_VISNIR_susc, valid_freq_range=Si3N4_VISNIR_range)

#------------------------------------------------------------------
# stoichiometric silicon nitride (Si3N4) from K. Luke, et. al., Broadband mid-infrared frequency comb generation in a Si3N4 microresonator, Opt. Lett. 40, 4823-4826, 2015
# fit from https://refractiveindex.info/?shelf=main&book=Si3N4&page=Luke
# wavelength range: 0.310 - 5.504 um

Si3N4_NIR_range = mp.FreqRange(min=um_scale/5.504, max=um_scale/0.310)

Si3N4_NIR_frq1 = 0.1353406
Si3N4_NIR_gam1 = 0
Si3N4_NIR_sig1 = 3.0249
Si3N4_NIR_frq2 = 1239.842
Si3N4_NIR_gam2 = 0
Si3N4_NIR_sig2 = 40314

Si3N4_NIR_susc = [mp.LorentzianSusceptibility(frequency=Si3N4_NIR_frq1, gamma=Si3N4_NIR_gam1, sigma=Si3N4_NIR_sig1),
                  mp.LorentzianSusceptibility(frequency=Si3N4_NIR_frq2, gamma=Si3N4_NIR_gam2, sigma=Si3N4_NIR_sig2)]

Si3N4_NIR = mp.Medium(epsilon=1.0, E_susceptibilities=Si3N4_NIR_susc, valid_freq_range=Si3N4_NIR_range)

#------------------------------------------------------------------

# elemental metals from A.D. Rakic et al., Applied Optics, Vol. 37, No. 22, pp. 5271-83, 1998
# wavelength range: 0.2 - 12.4 um

metal_range = mp.FreqRange(min=um_scale/12.4, max=um_scale/0.2)

# silver (Ag)

Ag_plasma_frq = 9.01*eV_um_scale
Ag_f0 = 0.845
Ag_frq0 = 1e-10
Ag_gam0 = 0.048*eV_um_scale
Ag_sig0 = Ag_f0*Ag_plasma_frq**2/Ag_frq0**2
Ag_f1 = 0.065
Ag_frq1 = 0.816*eV_um_scale      # 1.519 um
Ag_gam1 = 3.886*eV_um_scale
Ag_sig1 = Ag_f1*Ag_plasma_frq**2/Ag_frq1**2
Ag_f2 = 0.124
Ag_frq2 = 4.481*eV_um_scale      # 0.273 um
Ag_gam2 = 0.452*eV_um_scale
Ag_sig2 = Ag_f2*Ag_plasma_frq**2/Ag_frq2**2
Ag_f3 = 0.011
Ag_frq3 = 8.185*eV_um_scale      # 0.152 um
Ag_gam3 = 0.065*eV_um_scale
Ag_sig3 = Ag_f3*Ag_plasma_frq**2/Ag_frq3**2
Ag_f4 = 0.840
Ag_frq4 = 9.083*eV_um_scale      # 0.137 um
Ag_gam4 = 0.916*eV_um_scale
Ag_sig4 = Ag_f4*Ag_plasma_frq**2/Ag_frq4**2

Ag_susc = [mp.DrudeSusceptibility(frequency=Ag_frq0, gamma=Ag_gam0, sigma=Ag_sig0),
           mp.LorentzianSusceptibility(frequency=Ag_frq1, gamma=Ag_gam1, sigma=Ag_sig1),
           mp.LorentzianSusceptibility(frequency=Ag_frq2, gamma=Ag_gam2, sigma=Ag_sig2),
           mp.LorentzianSusceptibility(frequency=Ag_frq3, gamma=Ag_gam3, sigma=Ag_sig3),
           mp.LorentzianSusceptibility(frequency=Ag_frq4, gamma=Ag_gam4, sigma=Ag_sig4)]

Ag = mp.Medium(epsilon=1.0, E_susceptibilities=Ag_susc, valid_freq_range=metal_range)

#------------------------------------------------------------------

# gold (Au)

Au_plasma_frq = 9.03*eV_um_scale
Au_f0 = 0.760
Au_frq0 = 1e-10
Au_gam0 = 0.053*eV_um_scale
Au_sig0 = Au_f0*Au_plasma_frq**2/Au_frq0**2
Au_f1 = 0.024
Au_frq1 = 0.415*eV_um_scale      # 2.988 um
Au_gam1 = 0.241*eV_um_scale
Au_sig1 = Au_f1*Au_plasma_frq**2/Au_frq1**2
Au_f2 = 0.010
Au_frq2 = 0.830*eV_um_scale      # 1.494 um
Au_gam2 = 0.345*eV_um_scale
Au_sig2 = Au_f2*Au_plasma_frq**2/Au_frq2**2
Au_f3 = 0.071
Au_frq3 = 2.969*eV_um_scale      # 0.418 um
Au_gam3 = 0.870*eV_um_scale
Au_sig3 = Au_f3*Au_plasma_frq**2/Au_frq3**2
Au_f4 = 0.601
Au_frq4 = 4.304*eV_um_scale      # 0.288 um
Au_gam4 = 2.494*eV_um_scale
Au_sig4 = Au_f4*Au_plasma_frq**2/Au_frq4**2

Au_susc = [mp.DrudeSusceptibility(frequency=Au_frq0, gamma=Au_gam0, sigma=Au_sig0),
           mp.LorentzianSusceptibility(frequency=Au_frq1, gamma=Au_gam1, sigma=Au_sig1),
           mp.LorentzianSusceptibility(frequency=Au_frq2, gamma=Au_gam2, sigma=Au_sig2),
           mp.LorentzianSusceptibility(frequency=Au_frq3, gamma=Au_gam3, sigma=Au_sig3),
           mp.LorentzianSusceptibility(frequency=Au_frq4, gamma=Au_gam4, sigma=Au_sig4)]

Au = mp.Medium(epsilon=1.0, E_susceptibilities=Au_susc, valid_freq_range=metal_range)

#------------------------------------------------------------------

# copper (Cu)

Cu_plasma_frq = 10.83*eV_um_scale
Cu_f0 = 0.575
Cu_frq0 = 1e-10
Cu_gam0 = 0.030*eV_um_scale
Cu_sig0 = Cu_f0*Cu_plasma_frq**2/Cu_frq0**2
Cu_f1 = 0.061
Cu_frq1 = 0.291*eV_um_scale      # 4.261 um
Cu_gam1 = 0.378*eV_um_scale
Cu_sig1 = Cu_f1*Cu_plasma_frq**2/Cu_frq1**2
Cu_f2 = 0.104
Cu_frq2 = 2.957*eV_um_scale      # 0.419 um
Cu_gam2 = 1.056*eV_um_scale
Cu_sig2 = Cu_f2*Cu_plasma_frq**2/Cu_frq2**2
Cu_f3 = 0.723
Cu_frq3 = 5.300*eV_um_scale      # 0.234 um
Cu_gam3 = 3.213*eV_um_scale
Cu_sig3 = Cu_f3*Cu_plasma_frq**2/Cu_frq3**2
Cu_f4 = 0.638
Cu_frq4 = 11.18*eV_um_scale      # 0.111 um
Cu_gam4 = 4.305*eV_um_scale
Cu_sig4 = Cu_f4*Cu_plasma_frq**2/Cu_frq4**2

Cu_susc = [mp.DrudeSusceptibility(frequency=Cu_frq0, gamma=Cu_gam0, sigma=Cu_sig0),
           mp.LorentzianSusceptibility(frequency=Cu_frq1, gamma=Cu_gam1, sigma=Cu_sig1),
           mp.LorentzianSusceptibility(frequency=Cu_frq2, gamma=Cu_gam2, sigma=Cu_sig2),
           mp.LorentzianSusceptibility(frequency=Cu_frq3, gamma=Cu_gam3, sigma=Cu_sig3),
           mp.LorentzianSusceptibility(frequency=Cu_frq4, gamma=Cu_gam4, sigma=Cu_sig4)]

Cu = mp.Medium(epsilon=1.0, E_susceptibilities=Cu_susc, valid_freq_range=metal_range)

#------------------------------------------------------------------

# aluminum (Al)

Al_plasma_frq = 14.98*eV_um_scale
Al_f0 = 0.523
Al_frq0 = 1e-10
Al_gam0 = 0.047*eV_um_scale
Al_sig0 = Al_f0*Al_plasma_frq**2/Al_frq0**2
Al_f1 = 0.227
Al_frq1 = 0.162*eV_um_scale      # 7.654 um
Al_gam1 = 0.333*eV_um_scale
Al_sig1 = Al_f1*Al_plasma_frq**2/Al_frq1**2
Al_f2 = 0.050
Al_frq2 = 1.544*eV_um_scale      # 0.803 um
Al_gam2 = 0.312*eV_um_scale
Al_sig2 = Al_f2*Al_plasma_frq**2/Al_frq2**2
Al_f3 = 0.166
Al_frq3 = 1.808*eV_um_scale      # 0.686 um
Al_gam3 = 1.351*eV_um_scale
Al_sig3 = Al_f3*Al_plasma_frq**2/Al_frq3**2
Al_f4 = 0.030
Al_frq4 = 3.473*eV_um_scale      # 0.357 um
Al_gam4 = 3.382*eV_um_scale
Al_sig4 = Al_f4*Al_plasma_frq**2/Al_frq4**2

Al_susc = [mp.DrudeSusceptibility(frequency=Al_frq0, gamma=Al_gam0, sigma=Al_sig0),
           mp.LorentzianSusceptibility(frequency=Al_frq1, gamma=Al_gam1, sigma=Al_sig1),
           mp.LorentzianSusceptibility(frequency=Al_frq2, gamma=Al_gam2, sigma=Al_sig2),
           mp.LorentzianSusceptibility(frequency=Al_frq3, gamma=Al_gam3, sigma=Al_sig3),
           mp.LorentzianSusceptibility(frequency=Al_frq4, gamma=Al_gam4, sigma=Al_sig4)]

Al = mp.Medium(epsilon=1.0, E_susceptibilities=Al_susc, valid_freq_range=metal_range)

#------------------------------------------------------------------

# beryllium (Be)

Be_plasma_frq = 18.51*eV_um_scale
Be_f0 = 0.084
Be_frq0 = 1e-10
Be_gam0 = 0.035*eV_um_scale
Be_sig0 = Be_f0*Be_plasma_frq**2/Be_frq0**2
Be_f1 = 0.031
Be_frq1 = 0.100*eV_um_scale     # 12.398 um
Be_gam1 = 1.664*eV_um_scale
Be_sig1 = Be_f1*Be_plasma_frq**2/Be_frq1**2
Be_f2 = 0.140
Be_frq2 = 1.032*eV_um_scale      # 1.201 um
Be_gam2 = 3.395*eV_um_scale
Be_sig2 = Be_f2*Be_plasma_frq**2/Be_frq2**2
Be_f3 = 0.530
Be_frq3 = 3.183*eV_um_scale      # 0.390 um
Be_gam3 = 4.454*eV_um_scale
Be_sig3 = Be_f3*Be_plasma_frq**2/Be_frq3**2
Be_f4 = 0.130
Be_frq4 = 4.604*eV_um_scale      # 0.269 um
Be_gam4 = 1.802*eV_um_scale
Be_sig4 = Be_f4*Be_plasma_frq**2/Be_frq4**2

Be_susc = [mp.DrudeSusceptibility(frequency=Be_frq0, gamma=Be_gam0, sigma=Be_sig0),
           mp.LorentzianSusceptibility(frequency=Be_frq1, gamma=Be_gam1, sigma=Be_sig1),
           mp.LorentzianSusceptibility(frequency=Be_frq2, gamma=Be_gam2, sigma=Be_sig2),
           mp.LorentzianSusceptibility(frequency=Be_frq3, gamma=Be_gam3, sigma=Be_sig3),
           mp.LorentzianSusceptibility(frequency=Be_frq4, gamma=Be_gam4, sigma=Be_sig4)]

Be = mp.Medium(epsilon=1.0, E_susceptibilities=Be_susc, valid_freq_range=metal_range)

#------------------------------------------------------------------

# chromium (Cr)

Cr_plasma_frq = 10.75*eV_um_scale
Cr_f0 = 0.168
Cr_frq0 = 1e-10
Cr_gam0 = 0.047*eV_um_scale
Cr_sig0 = Cr_f0*Cr_plasma_frq**2/Cr_frq0**2
Cr_f1 = 0.151
Cr_frq1 = 0.121*eV_um_scale     # 10.247 um
Cr_gam1 = 3.175*eV_um_scale
Cr_sig1 = Cr_f1*Cr_plasma_frq**2/Cr_frq1**2
Cr_f2 = 0.150
Cr_frq2 = 0.543*eV_um_scale      # 2.283 um
Cr_gam2 = 1.305*eV_um_scale
Cr_sig2 = Cr_f2*Cr_plasma_frq**2/Cr_frq2**2
Cr_f3 = 1.149
Cr_frq3 = 1.970*eV_um_scale      # 0.629 um
Cr_gam3 = 2.676*eV_um_scale
Cr_sig3 = Cr_f3*Cr_plasma_frq**2/Cr_frq3**2
Cr_f4 = 0.825
Cr_frq4 = 8.775*eV_um_scale      # 0.141 um
Cr_gam4 = 1.335*eV_um_scale
Cr_sig4 = Cr_f4*Cr_plasma_frq**2/Cr_frq4**2

Cr_susc = [mp.DrudeSusceptibility(frequency=Cr_frq0, gamma=Cr_gam0, sigma=Cr_sig0),
           mp.LorentzianSusceptibility(frequency=Cr_frq1, gamma=Cr_gam1, sigma=Cr_sig1),
           mp.LorentzianSusceptibility(frequency=Cr_frq2, gamma=Cr_gam2, sigma=Cr_sig2),
           mp.LorentzianSusceptibility(frequency=Cr_frq3, gamma=Cr_gam3, sigma=Cr_sig3),
           mp.LorentzianSusceptibility(frequency=Cr_frq4, gamma=Cr_gam4, sigma=Cr_sig4)]

Cr = mp.Medium(epsilon=1.0, E_susceptibilities=Cr_susc, valid_freq_range=metal_range)

#------------------------------------------------------------------

# nickel (Ni)

Ni_plasma_frq = 15.92*eV_um_scale
Ni_f0 = 0.096
Ni_frq0 = 1e-10
Ni_gam0 = 0.048*eV_um_scale
Ni_sig0 = Ni_f0*Ni_plasma_frq**2/Ni_frq0**2
Ni_f1 = 0.100
Ni_frq1 = 0.174*eV_um_scale      # 7.126 um
Ni_gam1 = 4.511*eV_um_scale
Ni_sig1 = Ni_f1*Ni_plasma_frq**2/Ni_frq1**2
Ni_f2 = 0.135
Ni_frq2 = 0.582*eV_um_scale      # 2.130 um
Ni_gam2 = 1.334*eV_um_scale
Ni_sig2 = Ni_f2*Ni_plasma_frq**2/Ni_frq2**2
Ni_f3 = 0.106
Ni_frq3 = 1.597*eV_um_scale      # 0.776 um
Ni_gam3 = 2.178*eV_um_scale
Ni_sig3 = Ni_f3*Ni_plasma_frq**2/Ni_frq3**2
Ni_f4 = 0.729
Ni_frq4 = 6.089*eV_um_scale      # 0.204 um
Ni_gam4 = 6.292*eV_um_scale
Ni_sig4 = Ni_f4*Ni_plasma_frq**2/Ni_frq4**2

Ni_susc = [mp.DrudeSusceptibility(frequency=Ni_frq0, gamma=Ni_gam0, sigma=Ni_sig0),
           mp.LorentzianSusceptibility(frequency=Ni_frq1, gamma=Ni_gam1, sigma=Ni_sig1),
           mp.LorentzianSusceptibility(frequency=Ni_frq2, gamma=Ni_gam2, sigma=Ni_sig2),
           mp.LorentzianSusceptibility(frequency=Ni_frq3, gamma=Ni_gam3, sigma=Ni_sig3),
           mp.LorentzianSusceptibility(frequency=Ni_frq4, gamma=Ni_gam4, sigma=Ni_sig4)]

Ni = mp.Medium(epsilon=1.0, E_susceptibilities=Ni_susc, valid_freq_range=metal_range)

#------------------------------------------------------------------

# palladium (Pd)

Pd_plasma_frq = 9.72*eV_um_scale
Pd_f0 = 0.330
Pd_frq0 = 1e-10
Pd_gam0 = 0.008*eV_um_scale
Pd_sig0 = Pd_f0*Pd_plasma_frq**2/Pd_frq0**2
Pd_f1 = 0.649
Pd_frq1 = 0.336*eV_um_scale      # 3.690 um
Pd_gam1 = 2.950*eV_um_scale
Pd_sig1 = Pd_f1*Pd_plasma_frq**2/Pd_frq1**2
Pd_f2 = 0.121
Pd_frq2 = 0.501*eV_um_scale      # 2.475 um
Pd_gam2 = 0.555*eV_um_scale
Pd_sig2 = Pd_f2*Pd_plasma_frq**2/Pd_frq2**2
Pd_f3 = 0.638
Pd_frq3 = 1.659*eV_um_scale      # 0.747 um
Pd_gam3 = 4.621*eV_um_scale
Pd_sig3 = Pd_f3*Pd_plasma_frq**2/Pd_frq3**2
Pd_f4 = 0.453
Pd_frq4 = 5.715*eV_um_scale      # 0.217 um
Pd_gam4 = 3.236*eV_um_scale
Pd_sig4 = Pd_f4*Pd_plasma_frq**2/Pd_frq4**2

Pd_susc = [mp.DrudeSusceptibility(frequency=Pd_frq0, gamma=Pd_gam0, sigma=Pd_sig0),
           mp.LorentzianSusceptibility(frequency=Pd_frq1, gamma=Pd_gam1, sigma=Pd_sig1),
           mp.LorentzianSusceptibility(frequency=Pd_frq2, gamma=Pd_gam2, sigma=Pd_sig2),
           mp.LorentzianSusceptibility(frequency=Pd_frq3, gamma=Pd_gam3, sigma=Pd_sig3),
           mp.LorentzianSusceptibility(frequency=Pd_frq4, gamma=Pd_gam4, sigma=Pd_sig4)]

Pd = mp.Medium(epsilon=1.0, E_susceptibilities=Pd_susc, valid_freq_range=metal_range)

#------------------------------------------------------------------

# platinum (Pt)

Pt_plasma_frq = 9.59*eV_um_scale
Pt_f0 = 0.333
Pt_frq0 = 1e-10
Pt_gam0 = 0.080*eV_um_scale
Pt_sig0 = Pt_f0*Pt_plasma_frq**2/Pt_frq0**2
Pt_f1 = 0.191
Pt_frq1 = 0.780*eV_um_scale      # 1.590 um
Pt_gam1 = 0.517*eV_um_scale
Pt_sig1 = Pt_f1*Pt_plasma_frq**2/Pt_frq1**2
Pt_f2 = 0.659
Pt_frq2 = 1.314*eV_um_scale      # 0.944 um
Pt_gam2 = 1.838*eV_um_scale
Pt_sig2 = Pt_f2*Pt_plasma_frq**2/Pt_frq2**2
Pt_f3 = 0.547
Pt_frq3 = 3.141*eV_um_scale      # 0.395 um
Pt_gam3 = 3.668*eV_um_scale
Pt_sig3 = Pt_f3*Pt_plasma_frq**2/Pt_frq3**2
Pt_f4 = 3.576
Pt_frq4 = 9.249*eV_um_scale      # 0.134 um
Pt_gam4 = 8.517*eV_um_scale
Pt_sig4 = Pt_f4*Pt_plasma_frq**2/Pt_frq4**2

Pt_susc = [mp.DrudeSusceptibility(frequency=Pt_frq0, gamma=Pt_gam0, sigma=Pt_sig0),
           mp.LorentzianSusceptibility(frequency=Pt_frq1, gamma=Pt_gam1, sigma=Pt_sig1),
           mp.LorentzianSusceptibility(frequency=Pt_frq2, gamma=Pt_gam2, sigma=Pt_sig2),
           mp.LorentzianSusceptibility(frequency=Pt_frq3, gamma=Pt_gam3, sigma=Pt_sig3),
           mp.LorentzianSusceptibility(frequency=Pt_frq4, gamma=Pt_gam4, sigma=Pt_sig4)]

Pt = mp.Medium(epsilon=1.0, E_susceptibilities=Pt_susc, valid_freq_range=metal_range)

#------------------------------------------------------------------

# titanium (Ti)

Ti_plasma_frq = 7.29*eV_um_scale
Ti_f0 = 0.148
Ti_frq0 = 1e-10
Ti_gam0 = 0.082*eV_um_scale
Ti_sig0 = Ti_f0*Ti_plasma_frq**2/Ti_frq0**2
Ti_f1 = 0.899
Ti_frq1 = 0.777*eV_um_scale      # 1.596 um
Ti_gam1 = 2.276*eV_um_scale
Ti_sig1 = Ti_f1*Ti_plasma_frq**2/Ti_frq1**2
Ti_f2 = 0.393
Ti_frq2 = 1.545*eV_um_scale      # 0.802 um
Ti_gam2 = 2.518*eV_um_scale
Ti_sig2 = Ti_f2*Ti_plasma_frq**2/Ti_frq2**2
Ti_f3 = 0.187
Ti_frq3 = 2.509*eV_um_scale      # 0.494 um
Ti_gam3 = 1.663*eV_um_scale
Ti_sig3 = Ti_f3*Ti_plasma_frq**2/Ti_frq3**2
Ti_f4 = 0.001
Ti_frq4 = 19.43*eV_um_scale      # 0.064 um
Ti_gam4 = 1.762*eV_um_scale
Ti_sig4 = Ti_f4*Ti_plasma_frq**2/Ti_frq4**2

Ti_susc = [mp.DrudeSusceptibility(frequency=Ti_frq0, gamma=Ti_gam0, sigma=Ti_sig0),
           mp.LorentzianSusceptibility(frequency=Ti_frq1, gamma=Ti_gam1, sigma=Ti_sig1),
           mp.LorentzianSusceptibility(frequency=Ti_frq2, gamma=Ti_gam2, sigma=Ti_sig2),
           mp.LorentzianSusceptibility(frequency=Ti_frq3, gamma=Ti_gam3, sigma=Ti_sig3),
           mp.LorentzianSusceptibility(frequency=Ti_frq4, gamma=Ti_gam4, sigma=Ti_sig4)]

Ti = mp.Medium(epsilon=1.0, E_susceptibilities=Ti_susc, valid_freq_range=metal_range)

#------------------------------------------------------------------

# tungsten (W)

W_plasma_frq = 13.22*eV_um_scale
W_f0 = 0.206
W_frq0 = 1e-10
W_gam0 = 0.064*eV_um_scale
W_sig0 = W_f0*W_plasma_frq**2/W_frq0**2
W_f1 = 0.054
W_frq1 = 1.004*eV_um_scale      # 1.235 um
W_gam1 = 0.530*eV_um_scale
W_sig1 = W_f1*W_plasma_frq**2/W_frq1**2
W_f2 = 0.166
W_frq2 = 1.917*eV_um_scale      # 0.647 um
W_gam2 = 1.281*eV_um_scale
W_sig2 = W_f2*W_plasma_frq**2/W_frq2**2
W_f3 = 0.706
W_frq3 = 3.580*eV_um_scale      # 0.346 um
W_gam3 = 3.332*eV_um_scale
W_sig3 = W_f3*W_plasma_frq**2/W_frq3**2
W_f4 = 2.590
W_frq4 = 7.498*eV_um_scale      # 0.165 um
W_gam4 = 5.836*eV_um_scale
W_sig4 = W_f4*W_plasma_frq**2/W_frq4**2

W_susc = [mp.DrudeSusceptibility(frequency=W_frq0, gamma=W_gam0, sigma=W_sig0),
          mp.LorentzianSusceptibility(frequency=W_frq1, gamma=W_gam1, sigma=W_sig1),
          mp.LorentzianSusceptibility(frequency=W_frq2, gamma=W_gam2, sigma=W_sig2),
          mp.LorentzianSusceptibility(frequency=W_frq3, gamma=W_gam3, sigma=W_sig3),
          mp.LorentzianSusceptibility(frequency=W_frq4, gamma=W_gam4, sigma=W_sig4)]

W = mp.Medium(epsilon=1.0, E_susceptibilities=W_susc, valid_freq_range=metal_range)

#------------------------------------------------------------------

# metals from D. Barchiesi and T. Grosges, J. Nanophotonics, Vol. 8, 08996, 2015
# wavelength range: 0.4 - 0.8 um

metal_visible_range = mp.FreqRange(min=um_scale/0.8, max=um_scale/0.4)

# gold (Au)
# fit to P.B. Johnson and R.W. Christy, Physical Review B, Vol. 6, pp. 4370-9, 1972
Au_JC_visible_frq0 = 1/(0.139779231751333*um_scale)
Au_JC_visible_gam0 = 1/(26.1269913352870*um_scale)
Au_JC_visible_sig0 = 1

Au_JC_visible_frq1 = 1/(0.404064525036786*um_scale)
Au_JC_visible_gam1 = 1/(1.12834046202759*um_scale)
Au_JC_visible_sig1 = 2.07118534879440

Au_JC_visible_susc = [mp.DrudeSusceptibility(frequency=Au_JC_visible_frq0, gamma=Au_JC_visible_gam0, sigma=Au_JC_visible_sig0),
                      mp.LorentzianSusceptibility(frequency=Au_JC_visible_frq1, gamma=Au_JC_visible_gam1, sigma=Au_JC_visible_sig1)]

Au_JC_visible = mp.Medium(epsilon=6.1599, E_susceptibilities=Au_JC_visible_susc)

#------------------------------------------------------------------

# gold (Au)
# fit to E.D. Palik, Handbook of Optical Constants, Academic Press, 1985 
Au_visible_frq0 = 1/(0.0473629248511456*um_scale)
Au_visible_gam0 = 1/(0.255476199605166*um_scale)
Au_visible_sig0 = 1

Au_visible_frq1 = 1/(0.800619321082804*um_scale)
Au_visible_gam1 = 1/(0.381870287531951*um_scale)
Au_visible_sig1 = -169.060953137985

Au_visible_susc = [mp.DrudeSusceptibility(frequency=Au_visible_frq0, gamma=Au_visible_gam0, sigma=Au_visible_sig0),
                   mp.LorentzianSusceptibility(frequency=Au_visible_frq1, gamma=Au_visible_gam1, sigma=Au_visible_sig1)]

Au_visible = mp.Medium(epsilon=0.6888, E_susceptibilities=Au_visible_susc, valid_freq_range=metal_visible_range)

#------------------------------------------------------------------

## WARNING: unstable; field divergence may occur

# silver (Au)
# fit to E.D. Palik, Handbook of Optical Constants, Academic Press, 1985 
Ag_visible_frq0 = 1/(0.142050162130618*um_scale)
Ag_visible_gam0 = 1/(18.0357292925015*um_scale)
Ag_visible_sig0 = 1

Ag_visible_frq1 = 1/(0.115692151792108*um_scale)
Ag_visible_gam1 = 1/(0.257794324096575*um_scale)
Ag_visible_sig1 = 3.74465275944019

Ag_visible_susc = [mp.DrudeSusceptibility(frequency=Ag_visible_frq0, gamma=Ag_visible_gam0, sigma=Ag_visible_sig0),
                   mp.LorentzianSusceptibility(frequency=Ag_visible_frq1, gamma=Ag_visible_gam1, sigma=Ag_visible_sig1)]

Ag_visible = mp.Medium(epsilon=0.0067526, E_susceptibilities=Ag_visible_susc, valid_freq_range=metal_visible_range)

#------------------------------------------------------------------

## WARNING: unstable; field divergence may occur

# aluminum (Al)
# fit to E.D. Palik, Handbook of Optical Constants, Academic Press, 1985 
Al_visible_frq0 = 1/(0.0625841659042985*um_scale)
Al_visible_gam0 = 1/(0.606007002962666*um_scale)
Al_visible_sig0 = 1

Al_visible_frq1 = 1/(0.528191199577075*um_scale)
Al_visible_gam1 = 1/(0.291862527666814*um_scale)
Al_visible_sig1 = -44.4456675577921

Al_visible_susc = [mp.DrudeSusceptibility(frequency=Al_visible_frq0, gamma=Al_visible_gam0, sigma=Al_visible_sig0),
                   mp.LorentzianSusceptibility(frequency=Al_visible_frq1, gamma=Al_visible_gam1, sigma=Al_visible_sig1)]

Al_visible = mp.Medium(epsilon=0.13313, E_susceptibilities=Al_visible_susc, valid_freq_range=metal_visible_range)

#------------------------------------------------------------------

# chroimium (Cr)
# fit to E.D. Palik, Handbook of Optical Constants, Academic Press, 1985 
Cr_visible_frq0 = 1/(0.118410119507342*um_scale)
Cr_visible_gam0 = 1/(0.628596264869804*um_scale)
Cr_visible_sig0 = 1

Cr_visible_frq1 = 1/(0.565709598452496*um_scale)
Cr_visible_gam1 = 1/(0.731117670900812*um_scale)
Cr_visible_sig1 = 13.2912419951294

Cr_visible_susc = [mp.DrudeSusceptibility(frequency=Cr_visible_frq0, gamma=Cr_visible_gam0, sigma=Cr_visible_sig0),
                   mp.LorentzianSusceptibility(frequency=Cr_visible_frq1, gamma=Cr_visible_gam1, sigma=Cr_visible_sig1)]

Cr_visible = mp.Medium(epsilon=2.7767, E_susceptibilities=Cr_visible_susc, valid_freq_range=metal_visible_range)

#------------------------------------------------------------------

## WARNING: unstable; field divergence may occur

# titanium (Ti)
# fit to E.D. Palik, Handbook of Optical Constants, Academic Press, 1985 
Ti_visible_frq0 = 1/(0.101331651921602*um_scale)
Ti_visible_gam0 = 1/(0.365743382258719*um_scale)
Ti_visible_sig0 = 1

Ti_visible_frq1 = 1/(4.56839173979216e-09*um_scale)
Ti_visible_gam1 = 1/(5.86441957443603e-10*um_scale)
Ti_visible_sig1 = 54742662.1963414

Ti_visible_susc = [mp.DrudeSusceptibility(frequency=Ti_visible_frq0, gamma=Ti_visible_gam0, sigma=Ti_visible_sig0),
                   mp.LorentzianSusceptibility(frequency=Ti_visible_frq1, gamma=Ti_visible_gam1, sigma=Ti_visible_sig1)]

Ti_visible = mp.Medium(epsilon=-5.4742e7, E_susceptibilities=Ti_visible_susc, valid_freq_range=metal_visible_range)

#------------------------------------------------------------------

# aluminum (Al) from Horiba Technical Note 09: Drude Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Drude_Dispersion_Model.pdf
# wavelength range: 0.19 - 0.83 um

Al_drude_range = mp.FreqRange(min=um_scale/0.83, max=um_scale/0.19)

Al_drude_frq = 1/(0.0789607648707171*um_scale)
Al_drude_gam = 1/(1.78138208333333*um_scale)
Al_drude_sig = 1

Al_drude_susc = [mp.DrudeSusceptibility(frequency=Al_drude_frq, gamma=Al_drude_gam, sigma=Al_drude_sig)]

Al_drude = mp.Medium(epsilon=1.0, E_susceptibilities=Al_drude_susc, valid_freq_range=Al_drude_range)

#------------------------------------------------------------------

# cobalt (Co) from Horiba Technical Note 09: Drude Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Drude_Dispersion_Model.pdf
# wavelength range: 0.26 - 1.65 um

Co_range = mp.FreqRange(min=um_scale/1.65, max=um_scale/0.26)

Co_frq = 1/(0.0789607648707171*um_scale)
Co_gam = 1/(0.213802712536644*um_scale)
Co_sig = 1

Co_susc = [mp.DrudeSusceptibility(frequency=Co_frq, gamma=Co_gam, sigma=Co_sig)]

Co = mp.Medium(epsilon=3.694, E_susceptibilities=Co_susc, valid_freq_range=Co_range)

#------------------------------------------------------------------

## WARNING: unstable; field divergence may occur

# molybdenum (Mo) from Horiba Technical Note 09: Drude Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Drude_Dispersion_Model.pdf
# wavelength range: 0.25 - 0.83 um

Mo_range = mp.FreqRange(min=um_scale/0.83, max=um_scale/0.25)

Mo_frq = 1/(0.0620790071099539*um_scale)
Mo_gam = 1/(0.148359690080172*um_scale)
Mo_sig = 1

Mo_susc = [mp.DrudeSusceptibility(frequency=Mo_frq, gamma=Mo_gam, sigma=Mo_sig)]

Mo = mp.Medium(epsilon=-1.366, E_susceptibilities=Mo_susc, valid_freq_range=Mo_range)

#------------------------------------------------------------------

# nickel chrome (NiCr) from Horiba Technical Note 09: Drude Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Drude_Dispersion_Model.pdf
# wavelength range: 0.25 - 0.83 um

NiCr_range = mp.FreqRange(min=um_scale/0.83, max=um_scale/0.25)

NiCr_frq = 1/(0.0868845080588648*um_scale)
NiCr_gam = 1/(0.308418390547264*um_scale)
NiCr_sig = 1

NiCr_susc = [ mp.DrudeSusceptibility(frequency=NiCr_frq, gamma=NiCr_gam, sigma=NiCr_sig) ]

NiCr = mp.Medium(epsilon=1.0, E_susceptibilities=NiCr_susc, valid_freq_range=NiCr_range)

#------------------------------------------------------------------

# nickel iron (NiFe) from Horiba Technical Note 09: Drude Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Drude_Dispersion_Model.pdf
# wavelength range: 0.25 - 0.83 um

NiFe_range = mp.FreqRange(min=um_scale/0.83, max=um_scale/0.25)

NiFe_frq = 1/(0.0838297450980392*um_scale)
NiFe_gam = 1/(0.259381156903766*um_scale)
NiFe_sig = 1

NiFe_susc = [mp.DrudeSusceptibility(frequency=NiFe_frq, gamma=NiFe_gam, sigma=NiFe_sig)]

NiFe = mp.Medium(epsilon=1.0, E_susceptibilities=NiFe_susc, valid_freq_range=NiFe_range)

#------------------------------------------------------------------

# titanium (Ti) from Horiba Technical Note 09: Drude Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Drude_Dispersion_Model.pdf
# wavelength range: 0.21 - 1.24 um

Ti_drude_range = mp.FreqRange(min=um_scale/1.24, max=um_scale/0.21)

Ti_drude_frq = 1/(0.113746966055046*um_scale)
Ti_drude_gam = 1/(0.490056098814229*um_scale)
Ti_drude_sig = 1

Ti_drude_susc = [mp.DrudeSusceptibility(frequency=Ti_drude_frq, gamma=Ti_drude_gam, sigma=Ti_drude_sig)]

Ti_drude = mp.Medium(epsilon=1.0, E_susceptibilities=Ti_drude_susc, valid_freq_range=Ti_drude_range)

#------------------------------------------------------------------

# silicon nitride (SiN), non-stoichiometric, from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.21 - 2.07 um

SiN_range = mp.FreqRange(min=um_scale/2.07, max=um_scale/0.21)

SiN_frq1 = 1/(0.190891752117013*um_scale)
SiN_gam1 = 1/(3.11518072864322*um_scale)
SiN_sig1 = 1.2650

SiN_susc = [mp.LorentzianSusceptibility(frequency=SiN_frq1, gamma=SiN_gam1, sigma=SiN_sig1)]

SiN = mp.Medium(epsilon=2.320, E_susceptibilities=SiN_susc, valid_freq_range=SiN_range)

#------------------------------------------------------------------

# silicon nitride (Si3N4), stoichiometric, from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.23 - 0.83 um

Si3N4_range = mp.FreqRange(min=um_scale/0.83, max=um_scale/0.23)

Si3N4_frq1 = 1/(0.389153148148148*um_scale)
Si3N4_gam1 = 1/(0.693811936205932*um_scale)
Si3N4_sig1 = 4.377

Si3N4_susc = [mp.LorentzianSusceptibility(frequency=Si3N4_frq1, gamma=Si3N4_gam1, sigma=Si3N4_sig1)]

Si3N4 = mp.Medium(epsilon=1.0, E_susceptibilities=Si3N4_susc, valid_freq_range=Si3N4_range)

#------------------------------------------------------------------

# silicon dioxide (SiO2) from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.25 - 1.77 um

SiO2_range = mp.FreqRange(min=um_scale/1.77, max=um_scale/0.25)

SiO2_frq1 = 1/(0.103320160833333*um_scale)
SiO2_gam1 = 1/(12.3984193000000*um_scale)
SiO2_sig1 = 1.12

SiO2_susc = [mp.LorentzianSusceptibility(frequency=SiO2_frq1, gamma=SiO2_gam1, sigma=SiO2_sig1)]

SiO2 = mp.Medium(epsilon=1.0, E_susceptibilities=SiO2_susc, valid_freq_range=SiO2_range)

#------------------------------------------------------------------

# indium phosphide (InP) from Handbook of Optics, 2nd edition, Vol. 2., McGraw-Hill, 1994
# ref: https://refractiveindex.info/?shelf=main&book=InP&page=Pettit
# wavelength range: 0.95 - 10 um

InP_range = mp.FreqRange(min=um_scale/10, max=um_scale/0.95)

InP_frq1 = 1/(0.6263*um_scale)
InP_gam1 = 0
InP_sig1 = 2.316
InP_frq2 = 1/(32.935*um_scale)
InP_gam2 = 0
InP_sig2 = 2.765

InP_susc = [mp.LorentzianSusceptibility(frequency=InP_frq1, gamma=InP_gam1, sigma=InP_sig1),
            mp.LorentzianSusceptibility(frequency=InP_frq2, gamma=InP_gam2, sigma=InP_sig2)]

InP = mp.Medium(epsilon=7.255, E_susceptibilities=InP_susc, valid_freq_range=InP_range)

#------------------------------------------------------------------

# germanium (Ge) from N. P. Barnes and M. S. Piltch, J. Optical Society America, Vol. 69, pp. 178-180, 1979
# ref: https://refractiveindex.info/?shelf=main&book=Ge&page=Icenogle
# wavelength range: 2.5 - 12 um

Ge_range = mp.FreqRange(min=um_scale/12, max=um_scale/2.5)

Ge_frq1 = 1/(0.6641159*um_scale)
Ge_gam1 = 0
Ge_sig1 = 6.7288
Ge_frq2 = 1/(62.210127*um_scale)
Ge_gam2 = 0
Ge_sig2 = 0.21307

Ge_susc = [mp.LorentzianSusceptibility(frequency=Ge_frq1, gamma=Ge_gam1, sigma=Ge_sig1),
           mp.LorentzianSusceptibility(frequency=Ge_frq2, gamma=Ge_gam2, sigma=Ge_sig2)]

Ge = mp.Medium(epsilon=9.28156, E_susceptibilities=Ge_susc, valid_freq_range=Ge_range)

#------------------------------------------------------------------

# silicon (Si) from C. D. Salzberg and J. J. Villa, , J. Optical Society America, Vol. 47, pp. 244-246, 1957
# ref: https://refractiveindex.info/?shelf=main&book=Si&page=Salzberg
# wavelength range: 1.36 - 11 um

Si_range = mp.FreqRange(min=um_scale/11, max=um_scale/1.36)

Si_frq1 = 1/(0.301516485*um_scale)
Si_gam1 = 0
Si_sig1 = 10.6684293
Si_frq2 = 1/(1.13475115*um_scale)
Si_gam2 = 0
Si_sig2 = 0.0030434748
Si_frq3 = 1/(1104*um_scale)
Si_gam3 = 0
Si_sig3 = 1.54133408

Si_susc = [mp.LorentzianSusceptibility(frequency=Si_frq1, gamma=Si_gam1, sigma=Si_sig1),
           mp.LorentzianSusceptibility(frequency=Si_frq2, gamma=Si_gam2, sigma=Si_sig2),
           mp.LorentzianSusceptibility(frequency=Si_frq3, gamma=Si_gam3, sigma=Si_sig3)]

Si = mp.Medium(epsilon=1.0, E_susceptibilities=Si_susc, valid_freq_range=Si_range)

#------------------------------------------------------------------

# poly(methyl methacrylate) (PMMA) from N. Sultanova et al., Acta Physica Polonica A, Vol. 116, pp. 585-7, 2009
# ref: https://refractiveindex.info/?shelf=organic&book=poly%28methyl_methacrylate%29&page=Sultanova
# wavelength range: 0.437-1.052 um

PMMA_range = mp.FreqRange(min=um_scale/1.052, max=um_scale/0.437)

PMMA_frq1 = 1/(0.106362587407415*um_scale)
PMMA_gam1 = 0
PMMA_sig1 = 1.1819

PMMA_susc = [mp.LorentzianSusceptibility(frequency=PMMA_frq1, gamma=PMMA_gam1, sigma=PMMA_sig1)]

PMMA = mp.Medium(epsilon=1.0, E_susceptibilities=PMMA_susc, valid_freq_range=PMMA_range)

#------------------------------------------------------------------

# polycarbonate (PC) from N. Sultanova et al., Acta Physica Polonica A, Vol. 116, pp. 585-7, 2009
# ref: https://refractiveindex.info/?shelf=organic&book=polycarbonate&page=Sultanova
# wavelength range: 0.437-1.052 um

PC_range = mp.FreqRange(min=um_scale/1.052, max=um_scale/0.437)

PC_frq1 = 1/(0.145958898324152*um_scale)
PC_gam1 = 0
PC_sig1 = 1.4182

PC_susc = [mp.LorentzianSusceptibility(frequency=PC_frq1, gamma=PC_gam1, sigma=PC_sig1)]

PC = mp.Medium(epsilon=1.0, E_susceptibilities=PC_susc, valid_freq_range=PC_range)

#------------------------------------------------------------------

# polystyrene (PS) from N. Sultanova et al., Acta Physica Polonica A, Vol. 116, pp. 585-7, 2009
# ref: https://refractiveindex.info/?shelf=organic&book=polystyren&page=Sultanova
# wavelength range: 0.437-1.052 um

PS_range = mp.FreqRange(min=um_scale/1.052, max=um_scale/0.437)

PS_frq1 = 1/(0.142182980697410*um_scale)
PS_gam1 = 0
PS_sig1 = 1.4435

PS_susc = [mp.LorentzianSusceptibility(frequency=PS_frq1, gamma=PS_gam1, sigma=PS_sig1)]

PS = mp.Medium(epsilon=1.0, E_susceptibilities=PS_susc, valid_freq_range=PS_range)

#------------------------------------------------------------------

# cellulose (CLS) from N. Sultanova et al., Acta Physica Polonica A, Vol. 116, pp. 585-7, 2009
# ref: https://refractiveindex.info/?shelf=organic&book=cellulose&page=Sultanova
# wavelength range: 0.437-1.052 um

CLS_range = mp.FreqRange(min=um_scale/1.052, max=um_scale/0.437)

CLS_frq1 = 1/(0.105294824184287*um_scale)
CLS_gam1 = 0
CLS_sig1 = 1.124

CLS_susc = [mp.LorentzianSusceptibility(frequency=CLS_frq1, gamma=CLS_gam1, sigma=CLS_sig1)]

CLS = mp.Medium(epsilon=1.0, E_susceptibilities=CLS_susc, valid_freq_range=CLS_range)

