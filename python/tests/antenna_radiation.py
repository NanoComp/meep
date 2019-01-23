from __future__ import division

import meep as mp
import math
import numpy as np
import unittest

## compute the Poynting flux of an Ez-polarized dipole point source
## from the fields in 3 arrangements:
##      (1) bounding box of the near fields
##      (2) bounding circle of the far fields
##      (3) bounding box of the far fields

class TestAntennaRadiation(unittest.TestCase):

    def test_farfield(self):
        resolution = 50
        sxy = 4
        dpml = 1
        cell = mp.Vector3(sxy+2*dpml,sxy+2*dpml,0)

        pml_layers = mp.PML(dpml)

        fcen = 1.0
        df = 0.4

        sources = mp.Source(src=mp.GaussianSource(fcen,fwidth=df),
                            center=mp.Vector3(),
                            component=mp.Ez)

        symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y)]

        sim = mp.Simulation(cell_size=cell,
                            resolution=resolution,
                            sources=[sources],
                            symmetries=symmetries,
                            boundary_layers=[pml_layers])

        nearfield_box = sim.add_near2far(fcen, 0, 1,
                                         mp.Near2FarRegion(mp.Vector3(y=0.5*sxy), size=mp.Vector3(sxy)),
                                         mp.Near2FarRegion(mp.Vector3(y=-0.5*sxy), size=mp.Vector3(sxy), weight=-1),
                                         mp.Near2FarRegion(mp.Vector3(0.5*sxy), size=mp.Vector3(y=sxy)),
                                         mp.Near2FarRegion(mp.Vector3(-0.5*sxy), size=mp.Vector3(y=sxy), weight=-1))

        flux_box = sim.add_flux(fcen, 0, 1,
                                mp.FluxRegion(mp.Vector3(y=0.5*sxy), size=mp.Vector3(sxy)),
                                mp.FluxRegion(mp.Vector3(y=-0.5*sxy), size=mp.Vector3(sxy), weight=-1),
                                mp.FluxRegion(mp.Vector3(0.5*sxy), size=mp.Vector3(y=sxy)),
                                mp.FluxRegion(mp.Vector3(-0.5*sxy), size=mp.Vector3(y=sxy), weight=-1))

        sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(), 1e-8))

        near_flux = mp.get_fluxes(flux_box)[0]

        r = 1000/fcen      # circle of radius 1000 wavelengths
        npts = 100         # number of points in [0,2*pi) range of angles
        E = np.zeros((npts,3),dtype=np.complex128)
        H = np.zeros((npts,3),dtype=np.complex128)
        for n in range(npts):
            ff = sim.get_farfield(nearfield_box,
                                  mp.Vector3(r*math.cos(2*math.pi*n/npts),
                                             r*math.sin(2*math.pi*n/npts)))
            E[n,:] = [np.conj(ff[j]) for j in range(3)]
            H[n,:] = [ff[j+3] for j in range(3)]

        Px = np.real(np.multiply(E[:,1],H[:,2])-np.multiply(E[:,2],H[:,1]))
        Py = np.real(np.multiply(E[:,2],H[:,0])-np.multiply(E[:,0],H[:,2]))
        Pz = np.real(np.multiply(E[:,0],H[:,1])-np.multiply(E[:,1],H[:,0]))
        Pr = np.sqrt(np.square(Px)+np.square(Py))
        far_flux_circle = np.sum(Pr)*2*np.pi*r/len(Pr)

        rr = 20/fcen      # square of side length 20 wavelengths
        far_flux_square = (nearfield_box.flux(mp.Y, mp.Volume(center=mp.Vector3(y=0.5*rr), size=mp.Vector3(rr)).swigobj, resolution)[0]
                           - nearfield_box.flux(mp.Y, mp.Volume(center=mp.Vector3(y=-0.5*rr), size=mp.Vector3(rr)).swigobj, resolution)[0]
                           + nearfield_box.flux(mp.X, mp.Volume(center=mp.Vector3(0.5*rr), size=mp.Vector3(y=rr)).swigobj, resolution)[0]
                           - nearfield_box.flux(mp.X, mp.Volume(center=mp.Vector3(-0.5*rr), size=mp.Vector3(y=rr)).swigobj, resolution)[0])

        print("flux:, {:.6f}, {:.6f}, {:.6f}".format(near_flux,far_flux_circle,far_flux_square))

        self.assertAlmostEqual(near_flux, far_flux_circle, places=2)
        self.assertAlmostEqual(far_flux_circle, far_flux_square, places=2)
        self.assertAlmostEqual(far_flux_square, near_flux, places=2)

if __name__ == '__main__':
    unittest.main()
