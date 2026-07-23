# Materials Library

import meep as mp

# default unit length is 1 μm
um_scale = 1.0

# conversion factor for eV to 1/μm [=1/hc]
eV_um_scale = um_scale / 1.23984193


# ------------------------------------------------------------------
# Helper builders
#
# Each helper assembles the susceptibility list and the mp.Medium so that the
# per-material definitions below stay concise. Arithmetic for the pole
# parameters is kept inline at the call sites, so the resulting Medium objects
# are numerically identical to spelling each susceptibility out by hand.


def _lorentzian(epsilon, freq_range, poles):
    """Isotropic medium from a list of (frequency, gamma, sigma) Lorentz poles."""
    susc = [
        mp.LorentzianSusceptibility(frequency=frq, gamma=gam, sigma=sig)
        for frq, gam, sig in poles
    ]
    return mp.Medium(
        epsilon=epsilon, E_susceptibilities=susc, valid_freq_range=freq_range
    )


def _drude_lorentz(epsilon, freq_range, drude, lorentzians=()):
    """Isotropic medium from one Drude pole plus optional Lorentz poles.

    `drude` and each entry of `lorentzians` are (frequency, gamma, sigma) tuples.
    """
    susc = [mp.DrudeSusceptibility(frequency=drude[0], gamma=drude[1], sigma=drude[2])]
    susc += [
        mp.LorentzianSusceptibility(frequency=frq, gamma=gam, sigma=sig)
        for frq, gam, sig in lorentzians
    ]
    return mp.Medium(
        epsilon=epsilon, E_susceptibilities=susc, valid_freq_range=freq_range
    )


def _metal_LD(plasma_eV, freq_range, oscillators):
    """Lorentz-Drude metal in the Rakic et al. parameterization.

    `oscillators` is a list of (f, frequency_eV, gamma_eV). The first entry is
    the Drude (intraband) term, evaluated at a near-zero frequency; the rest are
    Lorentzian (interband) terms. Each sigma is derived from the plasma frequency
    as sigma = f * plasma_frq**2 / frq**2.
    """
    plasma_frq = plasma_eV * eV_um_scale
    susc = []
    for i, (f, frq_eV, gam_eV) in enumerate(oscillators):
        if i == 0:
            frq = 1e-10
            sig = f * plasma_frq**2 / frq**2
            susc.append(
                mp.DrudeSusceptibility(
                    frequency=frq, gamma=gam_eV * eV_um_scale, sigma=sig
                )
            )
        else:
            frq = frq_eV * eV_um_scale
            sig = f * plasma_frq**2 / frq**2
            susc.append(
                mp.LorentzianSusceptibility(
                    frequency=frq, gamma=gam_eV * eV_um_scale, sigma=sig
                )
            )
    return mp.Medium(epsilon=1.0, E_susceptibilities=susc, valid_freq_range=freq_range)


def _uniaxial(freq_range, o_poles, e_poles, epsilon=None, epsilon_diag=None):
    """Uniaxial (birefringent) medium with ordinary (o) axes in X and Y and the
    extraordinary (e) axis in Z.

    `o_poles` and `e_poles` are lists of (frequency, gamma, sigma) tuples.
    """
    susc = [
        mp.LorentzianSusceptibility(
            frequency=frq, gamma=gam, sigma_diag=sig * mp.Vector3(1, 1, 0)
        )
        for frq, gam, sig in o_poles
    ]
    susc += [
        mp.LorentzianSusceptibility(
            frequency=frq, gamma=gam, sigma_diag=sig * mp.Vector3(0, 0, 1)
        )
        for frq, gam, sig in e_poles
    ]
    if epsilon_diag is not None:
        kwargs = {"epsilon_diag": epsilon_diag}
    else:
        kwargs = {"epsilon": epsilon}
    return mp.Medium(E_susceptibilities=susc, valid_freq_range=freq_range, **kwargs)


# ------------------------------------------------------------------
# crystalline silicon (c-Si) from A. Deinega et al., J. Optical Society of America A, Vol. 28, No. 5, pp. 770-77 (2011)
# based on experimental data for intrinsic silicon at T=300K from M.A. Green and M. Keevers, Progress in Photovoltaics, Vol. 3, pp. 189-92 (1995)
# wavelength range: 0.4 - 1.0 μm

cSi = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale, max=um_scale / 0.4),
    [
        (3.64 / um_scale, 0, 8),
        (2.76 / um_scale, 2 * 0.063 / um_scale, 2.85),
        (1.73 / um_scale, 2 * 2.5 / um_scale, -0.107),
    ],
)

# ------------------------------------------------------------------
# amorphous silicon (a-Si) from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.21 - 0.83 μm

aSi = _lorentzian(
    3.109,
    mp.FreqRange(min=um_scale / 0.83, max=um_scale / 0.21),
    [(1 / (0.315481407124682 * um_scale), 1 / (0.645751005208333 * um_scale), 14.571)],
)

# ------------------------------------------------------------------
# hydrogenated amorphous silicon (a-Si:H) from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.21 - 0.83 μm

aSi_H = _lorentzian(
    3.22,
    mp.FreqRange(min=um_scale / 0.83, max=um_scale / 0.21),
    [(1 / (0.334189199460916 * um_scale), 1 / (0.579365387850467 * um_scale), 12.31)],
)

# ------------------------------------------------------------------
# indium tin oxide (ITO) from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.21 - 0.83 μm

ITO = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale / 0.83, max=um_scale / 0.21),
    [(1 / (0.182329695588235 * um_scale), 1 / (1.94637665620094 * um_scale), 2.5)],
)

# ------------------------------------------------------------------
# alumina (Al2O3) from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.21 - 2.07 μm

Al2O3 = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale / 2.07, max=um_scale / 0.21),
    [(1 / (0.101476668030774 * um_scale), 0, 1.52)],
)

# ------------------------------------------------------------------
# aluminum nitride (AlN) from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.26 - 1.65 μm

AlN = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale / 1.65, max=um_scale / 0.26),
    [(1 / (0.139058089950651 * um_scale), 0, 3.306)],
)

# ------------------------------------------------------------------
# aluminum arsenide (AlAs) from R.E. Fern and A. Onton, J. Applied Physics, Vol. 42, pp. 3499-500 (1971)
# ref: https://refractiveindex.info/?shelf=main&book=AlAs&page=Fern
# wavelength range: 0.56 - 2.2 μm

AlAs = _lorentzian(
    2.0792,
    mp.FreqRange(min=um_scale / 2.2, max=um_scale / 0.56),
    [
        (1 / (0.2822 * um_scale), 0, 6.0840),
        (1 / (27.62 * um_scale), 0, 1.900),
    ],
)

# ------------------------------------------------------------------
# borosilicate glass (BK7) from SCHOTT Zemax catalog 2017-01-20b
# ref: https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT
# wavelength range: 0.3 - 2.5 μm

BK7 = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale / 2.5, max=um_scale / 0.3),
    [
        (1 / (0.07746417668832478 * um_scale), 0, 1.03961212),
        (1 / (0.14148467902921502 * um_scale), 0, 0.231792344),
        (1 / (10.176475470417055 * um_scale), 0, 1.01046945),
    ],
)

# ------------------------------------------------------------------
# fused quartz (silica) from I.H. Malitson, J. Optical Society of America, Vol. 55, pp. 1205-9 (1965)
# ref: https://refractiveindex.info/?shelf=glass&book=fused_silica&page=Malitson
# wavelength range: 0.21 - 6.7 μm

fused_quartz = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale / 6.7, max=um_scale / 0.21),
    [
        (1 / (0.0684043 * um_scale), 0, 0.696166300),
        (1 / (0.1162414 * um_scale), 0, 0.407942600),
        (1 / (9.896161 * um_scale), 0, 0.897479400),
    ],
)

# ------------------------------------------------------------------
# gallium arsenide (GaAs) from T. Skauli et al., J. Applied Physics, Vol. 94, pp. 6447-55 (2003)
# ref: https://refractiveindex.info/?shelf=main&book=GaAs&page=Skauli
# wavelength range: 0.97 - 17 μm

GaAs = _lorentzian(
    5.372514,
    mp.FreqRange(min=um_scale / 17, max=um_scale / 0.97),
    [
        (1 / (0.4431307 * um_scale), 0, 5.466742),
        (1 / (0.8746453 * um_scale), 0, 0.02429960),
        (1 / (36.9166 * um_scale), 0, 1.957522),
    ],
)

# ------------------------------------------------------------------
# silicon nitride (Si3N4) from H. R. Philipp, J. Electrochemical Society 120, 295-300 (1973)
# ref: https://refractiveindex.info/?shelf=main&book=Si3N4&page=Philipp
# wavelength range: 0.207 - 1.24 μm

Si3N4_VISNIR = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale / 1.24, max=um_scale / 0.207),
    [(1 / (0.13967 * um_scale), 0, 2.8939)],
)

# ------------------------------------------------------------------
# silicon nitride (Si3N4) from K. Luke, et. al., Optics Letters, Vol. 40, pp. 4823-26 (2015)
# ref: https://refractiveindex.info/?shelf=main&book=Si3N4&page=Luke
# wavelength range: 0.310 - 5.504 μm

Si3N4_NIR = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale / 5.504, max=um_scale / 0.310),
    [
        (1 / (0.1353406 * um_scale), 0, 3.0249),
        (1 / (1239.842 * um_scale), 0, 40314),
    ],
)

# ------------------------------------------------------------------
# elemental metals from A.D. Rakic et al., Applied Optics, Vol. 37, No. 22, pp. 5271-83 (1998)
# wavelength range: 0.2 - 12.4 μm

metal_range_Ag = mp.FreqRange(min=um_scale / 12.398, max=um_scale / 0.24797)
metal_range_Au = mp.FreqRange(min=um_scale / 6.1992, max=um_scale / 0.24797)
metal_range = mp.FreqRange(min=um_scale / 12.398, max=um_scale / 0.20664)

# silver (Ag)
Ag = _metal_LD(
    9.01,
    metal_range_Ag,
    [
        (0.845, 0, 0.048),  # Drude
        (0.065, 0.816, 3.886),  # 1.519 μm
        (0.124, 4.481, 0.452),  # 0.273 μm
        (0.011, 8.185, 0.065),  # 0.152 μm
        (0.840, 9.083, 0.916),  # 0.137 μm
        (5.646, 20.29, 2.419),  # 0.061 μm
    ],
)

# gold (Au)
Au = _metal_LD(
    9.03,
    metal_range_Au,
    [
        (0.760, 0, 0.053),  # Drude
        (0.024, 0.415, 0.241),  # 2.988 μm
        (0.010, 0.830, 0.345),  # 1.494 μm
        (0.071, 2.969, 0.870),  # 0.418 μm
        (0.601, 4.304, 2.494),  # 0.288 μm
        (4.384, 13.32, 2.214),  # 0.093 μm
    ],
)

# copper (Cu)
Cu = _metal_LD(
    10.83,
    metal_range,
    [
        (0.575, 0, 0.03),  # Drude
        (0.061, 0.291, 0.378),  # 4.261 μm
        (0.104, 2.957, 1.056),  # 0.419 μm
        (0.723, 5.300, 3.213),  # 0.234 μm
        (0.638, 11.18, 4.305),  # 0.111 μm
    ],
)

# aluminum (Al)
Al = _metal_LD(
    14.98,
    metal_range,
    [
        (0.523, 0, 0.047),  # Drude
        (0.227, 0.162, 0.333),  # 7.654 μm
        (0.050, 1.544, 0.312),  # 0.803 μm
        (0.166, 1.808, 1.351),  # 0.686 μm
        (0.030, 3.473, 3.382),  # 0.357 μm
    ],
)

# beryllium (Be)
Be = _metal_LD(
    18.51,
    metal_range,
    [
        (0.084, 0, 0.035),  # Drude
        (0.031, 0.100, 1.664),  # 12.398 μm
        (0.140, 1.032, 3.395),  # 1.201 μm
        (0.530, 3.183, 4.454),  # 0.390 μm
        (0.130, 4.604, 1.802),  # 0.269 μm
    ],
)

# chromium (Cr)
Cr = _metal_LD(
    10.75,
    metal_range,
    [
        (0.168, 0, 0.047),  # Drude
        (0.151, 0.121, 3.175),  # 10.247 μm
        (0.150, 0.543, 1.305),  # 2.283 μm
        (1.149, 1.970, 2.676),  # 0.629 μm
        (0.825, 8.775, 1.335),  # 0.141 μm
    ],
)

# nickel (Ni)
Ni = _metal_LD(
    15.92,
    metal_range,
    [
        (0.096, 0, 0.048),  # Drude
        (0.100, 0.174, 4.511),  # 7.126 μm
        (0.135, 0.582, 1.334),  # 2.130 μm
        (0.106, 1.597, 2.178),  # 0.776 μm
        (0.729, 6.089, 6.292),  # 0.204 μm
    ],
)

# palladium (Pd)
Pd = _metal_LD(
    9.72,
    metal_range,
    [
        (0.330, 0, 0.008),  # Drude
        (0.649, 0.336, 2.950),  # 3.690 μm
        (0.121, 0.501, 0.555),  # 2.475 μm
        (0.638, 1.659, 4.621),  # 0.747 μm
        (0.453, 5.715, 3.236),  # 0.217 μm
    ],
)

# platinum (Pt)
Pt = _metal_LD(
    9.59,
    metal_range,
    [
        (0.333, 0, 0.080),  # Drude
        (0.191, 0.780, 0.517),  # 1.590 μm
        (0.659, 1.314, 1.838),  # 0.944 μm
        (0.547, 3.141, 3.668),  # 0.395 μm
        (3.576, 9.249, 8.517),  # 0.134 μm
    ],
)

# titanium (Ti)
Ti = _metal_LD(
    7.29,
    metal_range,
    [
        (0.148, 0, 0.082),  # Drude
        (0.899, 0.777, 2.276),  # 1.596 μm
        (0.393, 1.545, 2.518),  # 0.802 μm
        (0.187, 2.509, 1.663),  # 0.494 μm
        (0.001, 19.43, 1.762),  # 0.064 μm
    ],
)

# tungsten (W)
W = _metal_LD(
    13.22,
    metal_range,
    [
        (0.206, 0, 0.064),  # Drude
        (0.054, 1.004, 0.530),  # 1.235 μm
        (0.166, 1.917, 1.281),  # 0.647 μm
        (0.706, 3.580, 3.332),  # 0.346 μm
        (2.590, 7.498, 5.836),  # 0.165 μm
    ],
)

# ------------------------------------------------------------------
# metals from D. Barchiesi and T. Grosges, J. Nanophotonics, Vol. 8, 083097 (2014)
# including Errata from J. Nanophotonics, Vol. 8, 089996 (2014)
# wavelength range: 0.4 - 0.8 μm

metal_visible_range = mp.FreqRange(min=um_scale / 0.8, max=um_scale / 0.4)

# gold (Au)
# fit to P.B. Johnson and R.W. Christy, Physical Review B, Vol. 6, pp. 4370-9 (1972)

Au_JC_visible = _drude_lorentz(
    6.1599,
    metal_visible_range,
    (1 / (0.139779231751333 * um_scale), 1 / (26.1269913352870 * um_scale), 1),
    [
        (
            1 / (0.404064525036786 * um_scale),
            1 / (1.12834046202759 * um_scale),
            2.07118534879440,
        )
    ],
)

# ------------------------------------------------------------------
# gold (Au)
# fit to E.D. Palik, Handbook of Optical Constants, Academic Press (1985)

Au_visible = _drude_lorentz(
    0.6888,
    metal_visible_range,
    (1 / (0.0473629248511456 * um_scale), 1 / (0.255476199605166 * um_scale), 1),
    [
        (
            1 / (0.800619321082804 * um_scale),
            1 / (0.381870287531951 * um_scale),
            -169.060953137985,
        )
    ],
)

# ------------------------------------------------------------------
## WARNING: unstable; field divergence may occur

# silver (Ag)
# fit to E.D. Palik, Handbook of Optical Constants, Academic Press (1985)

Ag_visible = _drude_lorentz(
    0.0067526,
    metal_visible_range,
    (1 / (0.142050162130618 * um_scale), 1 / (18.0357292925015 * um_scale), 1),
    [
        (
            1 / (0.115692151792108 * um_scale),
            1 / (0.257794324096575 * um_scale),
            3.74465275944019,
        )
    ],
)

# ------------------------------------------------------------------
## WARNING: unstable; field divergence may occur

# aluminum (Al)
# fit to E.D. Palik, Handbook of Optical Constants, Academic Press (1985)

Al_visible = _drude_lorentz(
    0.13313,
    metal_visible_range,
    (1 / (0.0625841659042985 * um_scale), 1 / (0.606007002962666 * um_scale), 1),
    [
        (
            1 / (0.528191199577075 * um_scale),
            1 / (0.291862527666814 * um_scale),
            -44.4456675577921,
        )
    ],
)

# ------------------------------------------------------------------
# chromium (Cr)
# fit to E.D. Palik, Handbook of Optical Constants, Academic Press (1985)

Cr_visible = _drude_lorentz(
    2.7767,
    metal_visible_range,
    (1 / (0.118410119507342 * um_scale), 1 / (0.628596264869804 * um_scale), 1),
    [
        (
            1 / (0.565709598452496 * um_scale),
            1 / (0.731117670900812 * um_scale),
            13.2912419951294,
        )
    ],
)

# ------------------------------------------------------------------
# titanium (Ti)
# fit to E.D. Palik, Handbook of Optical Constants, Academic Press (1985)

Ti_visible = _drude_lorentz(
    2.17069,
    metal_visible_range,
    (1 / (0.2213799986964903 * um_scale), 1 / (12.815176733218491 * um_scale), 1),
    [
        (
            1 / (0.6425774603632576 * um_scale),
            1 / (0.1738808794709548 * um_scale),
            74.4496,
        )
    ],
)

# ------------------------------------------------------------------
# aluminum (Al) from Horiba Technical Note 09: Drude Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Drude_Dispersion_Model.pdf
# wavelength range: 0.19 - 0.83 μm

Al_drude = _drude_lorentz(
    1.0,
    mp.FreqRange(min=um_scale / 0.83, max=um_scale / 0.19),
    (1 / (0.0789607648707171 * um_scale), 1 / (1.78138208333333 * um_scale), 1),
)

# ------------------------------------------------------------------
# cobalt (Co) from Horiba Technical Note 09: Drude Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Drude_Dispersion_Model.pdf
# wavelength range: 0.26 - 1.65 μm

Co = _drude_lorentz(
    3.694,
    mp.FreqRange(min=um_scale / 1.65, max=um_scale / 0.26),
    (1 / (0.0789607648707171 * um_scale), 1 / (0.213802712536644 * um_scale), 1),
)

# ------------------------------------------------------------------
## WARNING: unstable; field divergence may occur

# molybdenum (Mo) from Horiba Technical Note 09: Drude Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Drude_Dispersion_Model.pdf
# wavelength range: 0.25 - 0.83 μm

Mo = _drude_lorentz(
    -1.366,
    mp.FreqRange(min=um_scale / 0.83, max=um_scale / 0.25),
    (1 / (0.0620790071099539 * um_scale), 1 / (0.148359690080172 * um_scale), 1),
)

# ------------------------------------------------------------------
# nickel chrome (NiCr) from Horiba Technical Note 09: Drude Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Drude_Dispersion_Model.pdf
# wavelength range: 0.25 - 0.83 μm

NiCr = _drude_lorentz(
    1.0,
    mp.FreqRange(min=um_scale / 0.83, max=um_scale / 0.25),
    (1 / (0.0868845080588648 * um_scale), 1 / (0.308418390547264 * um_scale), 1),
)

# ------------------------------------------------------------------
# nickel iron (NiFe) from Horiba Technical Note 09: Drude Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Drude_Dispersion_Model.pdf
# wavelength range: 0.25 - 0.83 μm

NiFe = _drude_lorentz(
    1.0,
    mp.FreqRange(min=um_scale / 0.83, max=um_scale / 0.25),
    (1 / (0.0838297450980392 * um_scale), 1 / (0.259381156903766 * um_scale), 1),
)

# ------------------------------------------------------------------
# titanium (Ti) from Horiba Technical Note 09: Drude Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Drude_Dispersion_Model.pdf
# wavelength range: 0.21 - 1.24 μm

Ti_drude = _drude_lorentz(
    1.0,
    mp.FreqRange(min=um_scale / 1.24, max=um_scale / 0.21),
    (1 / (0.113746966055046 * um_scale), 1 / (0.490056098814229 * um_scale), 1),
)

# ------------------------------------------------------------------
# silicon nitride (SiN), non-stoichiometric, from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.21 - 2.07 μm

SiN = _lorentzian(
    2.320,
    mp.FreqRange(min=um_scale / 2.07, max=um_scale / 0.21),
    [(1 / (0.190891752117013 * um_scale), 1 / (3.11518072864322 * um_scale), 1.2650)],
)

# ------------------------------------------------------------------
# silicon nitride (Si3N4), stoichiometric, from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.23 - 0.83 μm

Si3N4 = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale / 0.83, max=um_scale / 0.23),
    [(1 / (0.389153148148148 * um_scale), 1 / (0.693811936205932 * um_scale), 4.377)],
)

# ------------------------------------------------------------------
# silicon dioxide (SiO2) from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.25 - 1.77 μm

SiO2 = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale / 1.77, max=um_scale / 0.25),
    [(1 / (0.103320160833333 * um_scale), 1 / (12.3984193000000 * um_scale), 1.12)],
)

# ------------------------------------------------------------------
# indium phosphide (InP) from Handbook of Optics, 2nd edition, Vol. 2, McGraw-Hill (1994)
# ref: https://refractiveindex.info/?shelf=main&book=InP&page=Pettit
# wavelength range: 0.95 - 10 μm

InP = _lorentzian(
    7.255,
    mp.FreqRange(min=um_scale / 10, max=um_scale / 0.95),
    [
        (1 / (0.6263 * um_scale), 0, 2.316),
        (1 / (32.935 * um_scale), 0, 2.765),
    ],
)

# ------------------------------------------------------------------
# germanium (Ge) from N. P. Barnes and M. S. Piltch, J. Optical Society America, Vol. 69, pp. 178-180 (1979)
# ref: https://refractiveindex.info/?shelf=main&book=Ge&page=Icenogle
# wavelength range: 2.5 - 12 μm

Ge = _lorentzian(
    9.28156,
    mp.FreqRange(min=um_scale / 12, max=um_scale / 2.5),
    [
        (1 / (0.6641159 * um_scale), 0, 6.7288),
        (1 / (62.210127 * um_scale), 0, 0.21307),
    ],
)

# ------------------------------------------------------------------
# silicon (Si) from C. D. Salzberg and J. J. Villa, , J. Optical Society America, Vol. 47, pp. 244-246 (1957)
# ref: https://refractiveindex.info/?shelf=main&book=Si&page=Salzberg
# wavelength range: 1.36 - 11 μm

Si = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale / 11, max=um_scale / 1.36),
    [
        (1 / (0.301516485 * um_scale), 0, 10.6684293),
        (1 / (1.13475115 * um_scale), 0, 0.0030434748),
        (1 / (1104 * um_scale), 0, 1.54133408),
    ],
)

# ------------------------------------------------------------------
# poly(methyl methacrylate) (PMMA) from N. Sultanova et al., Acta Physica Polonica A, Vol. 116, pp. 585-7 (2009)
# ref: https://refractiveindex.info/?shelf=organic&book=poly%28methyl_methacrylate%29&page=Sultanova
# wavelength range: 0.437 - 1.052 μm

PMMA = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale / 1.052, max=um_scale / 0.437),
    [(1 / (0.106362587407415 * um_scale), 0, 1.1819)],
)

# ------------------------------------------------------------------
# polycarbonate (PC) from N. Sultanova et al., Acta Physica Polonica A, Vol. 116, pp. 585-7 (2009)
# ref: https://refractiveindex.info/?shelf=organic&book=polycarbonate&page=Sultanova
# wavelength range: 0.437 - 1.052 μm

PC = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale / 1.052, max=um_scale / 0.437),
    [(1 / (0.145958898324152 * um_scale), 0, 1.4182)],
)

# ------------------------------------------------------------------
# polystyrene (PS) from N. Sultanova et al., Acta Physica Polonica A, Vol. 116, pp. 585-7 (2009)
# ref: https://refractiveindex.info/?shelf=organic&book=polystyren&page=Sultanova
# wavelength range: 0.437 - 1.052 μm

PS = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale / 1.052, max=um_scale / 0.437),
    [(1 / (0.142182980697410 * um_scale), 0, 1.4435)],
)

# ------------------------------------------------------------------
# cellulose (CLS) from N. Sultanova et al., Acta Physica Polonica A, Vol. 116, pp. 585-7 (2009)
# ref: https://refractiveindex.info/?shelf=organic&book=cellulose&page=Sultanova
# wavelength range: 0.437 - 1.052 μm

CLS = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale / 1.052, max=um_scale / 0.437),
    [(1 / (0.105294824184287 * um_scale), 0, 1.124)],
)

# ------------------------------------------------------------------
# barium borate (BaB2O4), beta phase, from G. Tamosauskas et al., Optical Materials Express, Vol. 8, pp. 1410-18 (2018)
# ref: https://refractiveindex.info/?shelf=main&book=BaB2O4&page=Tamosauskas-o
# ref: https://refractiveindex.info/?shelf=main&book=BaB2O4&page=Tamosauskas-e
# wavelength range: 0.188 - 5.2 μm

## NOTE: ordinary (o) axes in X and Y, extraordinary (e) axis in Z

BaB2O4 = _uniaxial(
    mp.FreqRange(min=um_scale / 5.2, max=um_scale / 0.188),
    o_poles=[
        (1 / (0.06265780079128216 * um_scale), 0, 0.90291),
        (1 / (0.13706202975295528 * um_scale), 0, 0.83155),
        (1 / (7.746612162745725 * um_scale), 0, 0.76536),
    ],
    e_poles=[
        (1 / (0.0845103543951864 * um_scale), 0, 1.151075),
        (1 / (0.15029970059850417 * um_scale), 0, 0.21803),
        (1 / (16.217274740226856 * um_scale), 0, 0.656),
    ],
    epsilon=1.0,
)

# ------------------------------------------------------------------
# lithium niobate (LiNbO3) from D.E. Zelmon et al., J. Optical Society of America B, Vol. 14, pp. 3319-22 (1997)
# ref: https://refractiveindex.info/?shelf=main&book=LiNbO3&page=Zelmon-o
# ref: https://refractiveindex.info/?shelf=main&book=LiNbO3&page=Zelmon-e
# wavelength range: 0.4 - 5.0 μm

## NOTE: ordinary (o) axes in X and Y, extraordinary (e) axis in Z

LiNbO3 = _uniaxial(
    mp.FreqRange(min=um_scale / 5.0, max=um_scale / 0.4),
    o_poles=[
        (1 / (0.13281566172707193 * um_scale), 0, 2.6734),
        (1 / (0.24318717071424636 * um_scale), 0, 1.2290),
        (1 / (21.78531615561271 * um_scale), 0, 12.614),
    ],
    e_poles=[
        (1 / (0.14307340773183533 * um_scale), 0, 2.9804),
        (1 / (0.2580697580112788 * um_scale), 0, 0.5981),
        (1 / (20.39803912144498 * um_scale), 0, 8.9543),
    ],
    epsilon=1.0,
)

# ------------------------------------------------------------------
# calcium tungstate (CaWO4) from W.L. Bond, J. Applied Physics, Vol. 36, pp. 1674-77 (1965)
# ref: https://refractiveindex.info/?shelf=main&book=CaWO4&page=Bond-o
# ref: https://refractiveindex.info/?shelf=main&book=CaWO4&page=Bond-e
# wavelength range: 0.45 - 4.0 μm

## NOTE: ordinary (o) axes in X and Y, extraordinary (e) axis in Z

CaWO4 = _uniaxial(
    mp.FreqRange(min=um_scale / 4.0, max=um_scale / 0.45),
    o_poles=[
        (1 / (0.1347 * um_scale), 0, 2.5493),
        (1 / (10.815 * um_scale), 0, 0.9200),
    ],
    e_poles=[
        (1 / (0.1379 * um_scale), 0, 2.6041),
        (1 / (21.371 * um_scale), 0, 4.1237),
    ],
    epsilon=1.0,
)

# ------------------------------------------------------------------
# calcium carbonate (CaCO3) from G. Ghosh, Optics Communication, Vol. 163, pp. 95-102 (1999)
# ref: https://refractiveindex.info/?shelf=main&book=CaCO3&page=Ghosh-o
# ref: https://refractiveindex.info/?shelf=main&book=CaCO3&page=Ghosh-e
# wavelength range: 0.204 - 2.172 μm

## NOTE: ordinary (o) axes in X and Y, extraordinary (e) axis in Z

CaCO3 = _uniaxial(
    mp.FreqRange(min=um_scale / 2.172, max=um_scale / 0.204),
    o_poles=[
        (1 / (0.13940057496294625 * um_scale), 0, 0.96464345),
        (1 / (10.954451150103322 * um_scale), 0, 1.82831454),
    ],
    e_poles=[
        (1 / (0.1032906302623815 * um_scale), 0, 0.82427830),
        (1 / (10.954451150103322 * um_scale), 0, 0.14429128),
    ],
    epsilon_diag=mp.Vector3(1.73358749, 1.73358749, 1.35859695),
)

# ------------------------------------------------------------------
# silicon dioxide (SiO2) from G. Ghosh, Optics Communication, Vol. 163, pp. 95-102 (1999)
# ref: https://refractiveindex.info/?shelf=main&book=SiO2&page=Ghosh-o
# ref: https://refractiveindex.info/?shelf=main&book=SiO2&page=Ghosh-e
# wavelength range: 0.198 - 2.0531 μm

## NOTE: ordinary (o) axes in X and Y, extraordinary (e) axis in Z

SiO2_aniso = _uniaxial(
    mp.FreqRange(min=um_scale / 2.0531, max=um_scale / 0.198),
    o_poles=[
        (1 / (0.10029257051247614 * um_scale), 0, 1.07044083),
        (1 / (10 * um_scale), 0, 1.10202242),
    ],
    e_poles=[
        (1 / (0.10104546699382412 * um_scale), 0, 1.09509924),
        (1 / (10 * um_scale), 0, 1.15662475),
    ],
    epsilon_diag=mp.Vector3(1.28604141, 1.28604141, 1.28851804),
)

# ------------------------------------------------------------------
# gallium nitride (GaN), alpha phase (wurtzite), from A.S. Barker Jr. and M. Ilegems, Physical Review B, Vol. 7, pp. 743-50 (1973)
# ref: https://refractiveindex.info/?shelf=main&book=GaN&page=Barker-o
# ref: https://refractiveindex.info/?shelf=main&book=GaN&page=Barker-e
# wavelength range: 0.35 - 10 μm

## NOTE: ordinary (o) axes in X and Y, extraordinary (e) axis in Z

GaN = _uniaxial(
    mp.FreqRange(min=um_scale / 10.0, max=um_scale / 0.35),
    o_poles=[
        (1 / (0.256 * um_scale), 0, 1.75),
        (1 / (17.86 * um_scale), 0, 4.1),
    ],
    e_poles=[(1 / (18.76 * um_scale), 0, 5.08)],
    epsilon_diag=mp.Vector3(3.6, 3.6, 5.35),
)

# ------------------------------------------------------------------
# aluminum nitride (AlN) from J. Pastrnak and L. Roskovcova, Physica Status Solidi, Vol. 14, K5-8 (1966)
# ref: https://refractiveindex.info/?shelf=main&book=AlN&page=Pastrnak-o
# ref: https://refractiveindex.info/?shelf=main&book=AlN&page=Pastrnak-e
# wavelength range: 0.22 - 5 μm

## NOTE: ordinary (o) axes in X and Y, extraordinary (e) axis in Z

AlN_aniso = _uniaxial(
    mp.FreqRange(min=um_scale / 5.0, max=um_scale / 0.22),
    o_poles=[
        (1 / (0.1715 * um_scale), 0, 1.3786),
        (1 / (15.03 * um_scale), 0, 3.861),
    ],
    e_poles=[
        (1 / (0.1746 * um_scale), 0, 1.6173),
        (1 / (15.03 * um_scale), 0, 4.139),
    ],
    epsilon_diag=mp.Vector3(3.1399, 3.1399, 3.0729),
)

# ------------------------------------------------------------------
# alumina/sapphire (Al2O3) from I.H. Malitson and M.J. Dodge, J. Optical Society of America, Vol. 62, pp. 1405 (1972)
# ref: https://refractiveindex.info/?shelf=main&book=Al2O3&page=Malitson-o
# ref: https://refractiveindex.info/?shelf=main&book=Al2O3&page=Malitson-e
# wavelength range: 0.2 - 5 μm

## NOTE: ordinary (o) axes in X and Y, extraordinary (e) axis in Z

Al2O3_aniso = _uniaxial(
    mp.FreqRange(min=um_scale / 5.0, max=um_scale / 0.2),
    o_poles=[
        (1 / (0.0726631 * um_scale), 0, 1.4313493),
        (1 / (0.1193242 * um_scale), 0, 0.65054713),
        (1 / (18.02825 * um_scale), 0, 5.3414021),
    ],
    e_poles=[
        (1 / (0.0740288 * um_scale), 0, 1.5039759),
        (1 / (0.1216529 * um_scale), 0, 0.55069141),
        (1 / (20.072248 * um_scale), 0, 6.5927379),
    ],
    epsilon=1,
)

# ------------------------------------------------------------------
# yttrium oxide (Y2O3) from Y. Nigara, Japanese J. of Applied Physics, Vol. 7, pp. 404-8 (1968)
# ref: https://refractiveindex.info/?shelf=main&book=Y2O3&page=Nigara
# wavelength range: 0.25 - 9.6 μm

Y2O3 = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale / 9.6, max=um_scale / 0.25),
    [
        (1 / (0.1387 * um_scale), 0, 2.578),
        (1 / (22.936 * um_scale), 0, 3.935),
    ],
)

# ------------------------------------------------------------------
# undoped yttrium aluminum garnet (YAG) from D.E. Zelmon et al., Applied Optics, Vol. 37, 4933-5 (1998)
# ref: https://refractiveindex.info/?shelf=main&book=Y3Al5O12&page=Zelmon
# wavelength range: 0.4 - 5.0 μm

YAG = _lorentzian(
    1.0,
    mp.FreqRange(min=um_scale / 5.0, max=um_scale / 0.4),
    [
        (1 / (0.1088577052853862 * um_scale), 0, 2.28200),
        (1 / (16.814695953242804 * um_scale), 0, 3.27644),
    ],
)

# ------------------------------------------------------------------
# cadmium telluride (CdTe) from D.T.F. Marple, J. Applied Physics, Vol. 35, pp. 539-42 (1964)
# ref: https://refractiveindex.info/?shelf=main&book=CdTe&page=Marple
# wavelength range: 0.86 - 2.5 μm

CdTe = _lorentzian(
    5.68,
    mp.FreqRange(min=um_scale / 2.5, max=um_scale / 0.86),
    [(1 / (0.6049793384901669 * um_scale), 0, 1.53)],
)


# ------------------------------------------------------------------
# Registry of all materials defined above, keyed by name. Enables programmatic
# lookup, e.g. materials_library["Si"] or get_material("Si").

materials_library = {
    name: obj for name, obj in list(globals().items()) if isinstance(obj, mp.Medium)
}


def get_material(name):
    """Return the mp.Medium registered in this library under `name`."""
    return materials_library[name]


__all__ = sorted(materials_library) + [
    "um_scale",
    "eV_um_scale",
    "materials_library",
    "get_material",
]
