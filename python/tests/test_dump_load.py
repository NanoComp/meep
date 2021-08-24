import itertools
import os
import re
import sys
import unittest
import warnings
import h5py
import numpy as np
import meep as mp

try:
    unicode
except NameError:
    unicode = str


class TestLoadDump(unittest.TestCase):

    fname_base = re.sub(r'\.py$', '', os.path.split(sys.argv[0])[1])
    fname = fname_base + '-ez-000200.00.h5'

    def setUp(self):
        print("Running {}".format(self._testMethodName))

    @classmethod
    def setUpClass(cls):
        cls.temp_dir = mp.make_output_directory()
        print("Saving temp files to dir: {}".format(cls.temp_dir))

    @classmethod
    def tearDownClass(cls):
        # mp.delete_directory(cls.temp_dir)
        pass

    # Tests various combinations of dumping/loading structure & chunk layout.
    def _load_dump_structure(self, chunk_file=False, chunk_sim=False, single_parallel_file=True):
        from meep.materials import Al
        resolution = 50
        cell = mp.Vector3(5, 5)
        sources = mp.Source(src=mp.GaussianSource(1, fwidth=0.2), center=mp.Vector3(), component=mp.Ez)
        one_by_one = mp.Vector3(1, 1, mp.inf)
        geometry = [mp.Block(material=Al, center=mp.Vector3(), size=one_by_one),
                    mp.Block(material=mp.Medium(epsilon=13), center=mp.Vector3(1), size=one_by_one)]
        pml_layers = [mp.PML(0.5)]

        symmetries = [mp.Mirror(mp.Y)]

        sim1 = mp.Simulation(resolution=resolution,
                             cell_size=cell,
                             boundary_layers=pml_layers,
                             geometry=geometry,
                             symmetries=symmetries,
                             sources=[sources])

        sample_point = mp.Vector3(0.12, -0.29)
        ref_field_points = []

        def get_ref_field_point(sim):
            p = sim.get_field_point(mp.Ez, sample_point)
            ref_field_points.append(p.real)

        sim1.run(mp.at_every(5, get_ref_field_point), until=50)

        dump_dirname = os.path.join(self.temp_dir, 'test_load_dump_structure')
        sim1.dump(dump_dirname, dump_structure=True, dump_fields=False, single_parallel_file=single_parallel_file)

        dump_chunk_fname = None
        chunk_layout = None
        if chunk_file:
            dump_chunk_fname = os.path.join(dump_dirname, 'chunk_layout.h5')
            sim1.dump_chunk_layout(dump_chunk_fname)
            chunk_layout = dump_chunk_fname
        if chunk_sim:
            chunk_layout = sim1

        sim = mp.Simulation(resolution=resolution,
                            cell_size=cell,
                            boundary_layers=pml_layers,
                            sources=[sources],
                            symmetries=symmetries,
                            chunk_layout=chunk_layout)
        sim.load(dump_dirname, load_structure=True, load_fields=False, single_parallel_file=single_parallel_file)
        field_points = []

        def get_field_point(sim):
            p = sim.get_field_point(mp.Ez, sample_point)
            field_points.append(p.real)

        sim.run(mp.at_every(5, get_field_point), until=50)

        for ref_pt, pt in zip(ref_field_points, field_points):
            self.assertAlmostEqual(ref_pt, pt)

    @unittest.skip("disable_temporarily_for_debugging")
    def test_load_dump_structure(self):
        self._load_dump_structure()

    @unittest.skipIf(not mp.with_mpi(), "MPI specific test")
    @unittest.skip("disable_temporarily_for_debugging")
    def test_load_dump_structure_sharded(self):
        self._load_dump_structure(single_parallel_file=False)

    @unittest.skip("disable_temporarily_for_debugging")
    def test_load_dump_chunk_layout_file(self):
        self._load_dump_structure(chunk_file=True)

    @unittest.skip("disable_temporarily_for_debugging")
    def test_load_dump_chunk_layout_sim(self):
        self._load_dump_structure(chunk_sim=True)

    # Tests dumping/loading of fields & structure.
    def _load_dump_fields(self, single_parallel_file=True):
        from meep.materials import Al
        resolution = 50
        cell = mp.Vector3(5, 5)
        sources = mp.Source(src=mp.GaussianSource(1, fwidth=0.2), center=mp.Vector3(), component=mp.Ez)
        one_by_one = mp.Vector3(1, 1, mp.inf)
        geometry = [mp.Block(material=Al, center=mp.Vector3(), size=one_by_one),
                    mp.Block(material=mp.Medium(epsilon=13), center=mp.Vector3(1), size=one_by_one)]
        pml_layers = [mp.PML(0.5)]

        symmetries = [mp.Mirror(mp.Y)]

        sim1 = mp.Simulation(resolution=resolution,
                             cell_size=cell,
                             boundary_layers=pml_layers,
                             geometry=geometry,
                             symmetries=symmetries,
                             sources=[sources])

        sample_point = mp.Vector3(0.12, -0.29)
        ref_field_points = {}

        def get_ref_field_point(sim):
            p = sim.get_field_point(mp.Ez, sample_point)
            ref_field_points[sim.meep_time()] = p.real

        # First run until t=25 and save structure/fields
        sim1.run(mp.at_every(5, get_ref_field_point), until=25)

        dump_dirname = os.path.join(self.temp_dir, 'test_load_dump_fields')
        sim1.dump(dump_dirname, dump_structure=True, dump_fields=True, single_parallel_file=single_parallel_file)

        # Then continue running another 25 until t=50
        sim1.run(mp.at_every(5, get_ref_field_point), until=25)
        print('ref_field_points = ' + str(ref_field_points))

        # Now create a new simulation and try restoring state.
        sim = mp.Simulation(resolution=resolution,
                            cell_size=cell,
                            boundary_layers=pml_layers,
                            sources=[sources],
                            symmetries=symmetries,
                            chunk_layout=sim1)
        # Just restore structure first.
        sim.load(dump_dirname, load_structure=True, load_fields=False, single_parallel_file=single_parallel_file)
        field_points = {}

        def get_field_point(sim):
            p = sim.get_field_point(mp.Ez, sample_point)
            field_points[sim.meep_time()] = p.real

        # Run until t=1 so that fields are initialized.
        sim.run(mp.at_every(5, get_field_point), until=1)

        # Now load the fields.
        sim.load(dump_dirname, load_structure=False, load_fields=True, single_parallel_file=single_parallel_file)
        sim.run(mp.at_every(5, get_field_point), until=25)
        print('field_points = ' + str(field_points))

        for t, v in field_points.items():
            self.assertAlmostEqual(ref_field_points[t], v)

    def test_load_dump_fields(self):
        self._load_dump_fields()

    @unittest.skipIf(not mp.with_mpi(), "MPI specific test")
    @unittest.skip("disable_temporarily_for_debugging")
    def test_load_dump_fields_sharded(self):
        self._load_dump_fields(single_parallel_file=False)


if __name__ == '__main__':
    unittest.main()
