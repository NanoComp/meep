from __future__ import division

import math
import unittest

import numpy as np
import meep as mp


class TestAntennaRadiation(unittest.TestCase):

    def setUp(self):
        resolution = 50
        sxy = 4
        dpml = 1
        cell = mp.Vector3(sxy + 2 * dpml, sxy + 2 * dpml, 0)

        pml_layers = mp.PML(dpml)

        fcen = 1.0
        df = 0.4
        self.src_cmpt = mp.Ez

        sources = mp.Source(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(),
            component=self.src_cmpt
        )

        if self.src_cmpt == mp.Ex:
            symmetries = [mp.Mirror(mp.Y)]
        elif self.src_cmpt == mp.Ey:
            symmetries = [mp.Mirror(mp.X)]
        elif self.src_cmpt == mp.Ez:
            symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y)]

        self.sim = mp.Simulation(
            cell_size=cell,
            resolution=resolution,
            sources=[sources],
            symmetries=symmetries,
            boundary_layers=[pml_layers]
        )

        self.nearfield = self.sim.add_near2far(
            fcen,
            0,
            1,
            mp.Near2FarRegion(mp.Vector3(0, 0.5 * sxy), size=mp.Vector3(sxy)),
            mp.Near2FarRegion(mp.Vector3(0, -0.5 * sxy), size=mp.Vector3(sxy), weight=-1.0),
            mp.Near2FarRegion(mp.Vector3(0.5 * sxy), size=mp.Vector3(0, sxy)),
            mp.Near2FarRegion(mp.Vector3(-0.5 * sxy), size=mp.Vector3(0, sxy), weight=-1.0)
        )

    def test_farfield(self):

        expected = [
            (0j, 0j,
             0.013561672901190156 - 0.014417985982674887j,
             -0.007972091889832573 + 0.008474039259808384j,
             -0.010972504173533907 + 0.011663526883728901j, 0j),
            (0j, 0j,
             0.013580605154131353 - 0.014437805790483897j,
             -0.00727747760452799 + 0.007735510886150803j,
             -0.011467448148229137 + 0.01218937654274057j, 0j),
            (0j, 0j,
             0.013623909088910181 - 0.0144248508451555j,
             -0.006563855061778871 + 0.0069485969930299695j,
             -0.011939764884917667 + 0.012639701096677738j, 0j),
            (0j, 0j,
             0.01366551237281673 - 0.014382588988355682j,
             -0.005818861980069385 + 0.006123274021053104j,
             -0.012366016128341942 + 0.013012805294365515j, 0j),
            (0j, 0j,
             0.013680914029322 - 0.014330104120260687j,
             -0.00503658263321274 + 0.005274864846985408j,
             -0.012721298907853191 + 0.013322780730790889j, 0j),
            (0j, 0j,
             0.013666069573084385 - 0.014288938849212545j,
             -0.004223325784559204 + 0.004415251197679649j,
             -0.01299830976630439 + 0.0135885339986728j, 0j),
            (0j, 0j,
             0.013636004904245425 - 0.014270749446931172j,
             -0.003391403389081595 + 0.003548797079942323j,
             -0.013208708739843122 + 0.013821337222637264j, 0j),
            (0j, 0j,
             0.013609442340491029 - 0.014273698854090303j,
             -0.0025503938754594265 + 0.0026744729386901805j,
             -0.013369492605877786 + 0.014019798855923278j, 0j),
            (0j, 0j,
             0.013595569402532206 - 0.014287439152148014j,
             -0.0017041572263699332 + 0.001790574472119783j,
             -0.013489487541122221 + 0.014173702684166192j, 0j),
            (0j, 0j,
             0.01359208053488311 - 0.014300667598824968j,
             -8.535506874711307e-4 + 8.978801355605307e-4j,
             -0.01356639388716244 + 0.01427136892007339j, 0j)
        ]

        self.sim.run(until_after_sources=mp.stop_when_fields_decayed(50, self.src_cmpt, mp.Vector3(), 1e-8))
        r = 1000
        npts = 100

        result = []
        for n in range(npts):
            ff = self.sim.get_farfield(
                self.nearfield,
                mp.Vector3(r * math.cos(2 * math.pi * (n / npts)),
                           r * math.sin(2 * math.pi * (n / npts))))
            result.append(ff)

        np.testing.assert_allclose(expected, result[-10:])

if __name__ == '__main__':
    unittest.main()
