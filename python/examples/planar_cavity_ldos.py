# Computes the Purcell enhancement factor of a horizontal dipole in a 3D
# homogeneous dielectric cavity with lossless metallic walls on two sides.
# The simulated result is compared with the analytic theory from
# I. Abram et al., IEEE J. Quantum Electronics, Vol. 34, pp. 71-76 (1998).

import meep as mp
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


resolution = 50  # pixels/μm
dpml = 0.5       # thickness of PML
L = 6.0          # length of non-PML region
n = 2.4          # refractive index of surrounding medium
wvl = 1.0        # wavelength (in vacuum)

fcen = 1/wvl
sources = [mp.Source(src=mp.GaussianSource(fcen,fwidth=0.2*fcen),
                     component=mp.Ex,
                     center=mp.Vector3())]

symmetries = [mp.Mirror(direction=mp.X,phase=-1),
              mp.Mirror(direction=mp.Y),
              mp.Mirror(direction=mp.Z)]


def bulk_ldos():
    s = L+2*dpml
    cell_size = mp.Vector3(s,s,s)

    pml_layers = [mp.PML(dpml)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=pml_layers,
                        sources=sources,
                        symmetries=symmetries,
                        default_material=mp.Medium(index=n))

    sim.run(mp.dft_ldos(fcen,0,1),
            until_after_sources=mp.stop_when_fields_decayed(20,
                                                            mp.Ex,
                                                            mp.Vector3(),
                                                            1e-6))

    return sim.ldos_data[0]


def cavity_ldos(sz):
    sxy = L+2*dpml
    cell_size = mp.Vector3(sxy,sxy,sz)

    boundary_layers = [mp.PML(dpml,direction=mp.X),
                       mp.PML(dpml,direction=mp.Y)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=boundary_layers,
                        sources=sources,
                        symmetries=symmetries,
                        default_material=mp.Medium(index=n))

    sim.run(mp.dft_ldos(fcen,0,1),
            until_after_sources=mp.stop_when_fields_decayed(20,
                                                            mp.Ex,
                                                            mp.Vector3(),
                                                            1e-6))

    return sim.ldos_data[0]


if __name__ == '__main__':
    ldos_bulk = bulk_ldos()
    print("ldos_bulk:, {:.6f}".format(ldos_bulk))

    # units of wavelength in medium
    cavity_thickness = np.arange(0.50,2.55,0.05)

    gap = cavity_thickness*wvl/n

    ldos_cavity = np.zeros(len(cavity_thickness))
    for idx,g in enumerate(gap):
        ldos_cavity[idx] = cavity_ldos(g)
        print("ldos_cavity:, {:.3f}, {:.6f}".format(g,ldos_cavity[idx]))

    # Purcell enhancement factor (relative to bulk medium)
    pe_meep = ldos_cavity/ldos_bulk

    # equation 7 of reference
    pe_theory = (3*np.fix(cavity_thickness+0.5)/(4*cavity_thickness) +
                 (4*np.power(np.fix(cavity_thickness+0.5),3) -
                  np.fix(cavity_thickness+0.5)) /
                 (16*np.power(cavity_thickness,3)))

    if mp.am_master():
        plt.plot(cavity_thickness,pe_meep,'b-',label='Meep')
        plt.plot(cavity_thickness,pe_theory,'r-',label='theory')
        plt.plot(cavity_thickness,np.ones(len(cavity_thickness)),'k--')
        plt.xlabel('cavity thickness')
        plt.ylabel('Purcell enhancement factor (relative to bulk)')
        plt.title("horizontal point dipole at λ=1.0 μm in a cavity with"
                  "\n n=2.4 and lossless metallic walls on two sides")
        plt.axis([0.5,2.5,0.4,3.1])
        plt.legend()
        plt.savefig('cavity_purcell_factor_vs_thickness',
                    bbox_inches='tight')
