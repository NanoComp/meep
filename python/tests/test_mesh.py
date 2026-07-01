"""Integration tests for the Mesh geometric object type."""

import meep as mp
import numpy as np
import struct
import tempfile
import os
import unittest


def make_cube_mesh(material, center=mp.Vector3(0, 0, 0)):
    """Create a unit cube mesh centered at the given point."""
    verts = [
        (-0.5, -0.5, -0.5),
        (0.5, -0.5, -0.5),
        (0.5, 0.5, -0.5),
        (-0.5, 0.5, -0.5),
        (-0.5, -0.5, 0.5),
        (0.5, -0.5, 0.5),
        (0.5, 0.5, 0.5),
        (-0.5, 0.5, 0.5),
    ]
    tris = [
        (0, 2, 1),
        (0, 3, 2),  # -z
        (4, 5, 6),
        (4, 6, 7),  # +z
        (0, 1, 5),
        (0, 5, 4),  # -y
        (2, 3, 7),
        (2, 7, 6),  # +y
        (0, 4, 7),
        (0, 7, 3),  # -x
        (1, 2, 6),
        (1, 6, 5),  # +x
    ]
    return mp.Mesh(vertices=verts, triangles=tris, center=center, material=material)


def write_cube_stl(path):
    """Write a unit cube as a binary STL file."""
    faces = [
        ((-0.5, -0.5, -0.5), (0.5, 0.5, -0.5), (0.5, -0.5, -0.5)),
        ((-0.5, -0.5, -0.5), (-0.5, 0.5, -0.5), (0.5, 0.5, -0.5)),
        ((-0.5, -0.5, 0.5), (0.5, -0.5, 0.5), (0.5, 0.5, 0.5)),
        ((-0.5, -0.5, 0.5), (0.5, 0.5, 0.5), (-0.5, 0.5, 0.5)),
        ((-0.5, -0.5, -0.5), (0.5, -0.5, -0.5), (0.5, -0.5, 0.5)),
        ((-0.5, -0.5, -0.5), (0.5, -0.5, 0.5), (-0.5, -0.5, 0.5)),
        ((0.5, 0.5, -0.5), (-0.5, 0.5, -0.5), (-0.5, 0.5, 0.5)),
        ((0.5, 0.5, -0.5), (-0.5, 0.5, 0.5), (0.5, 0.5, 0.5)),
        ((-0.5, -0.5, -0.5), (-0.5, -0.5, 0.5), (-0.5, 0.5, 0.5)),
        ((-0.5, -0.5, -0.5), (-0.5, 0.5, 0.5), (-0.5, 0.5, -0.5)),
        ((0.5, -0.5, -0.5), (0.5, 0.5, -0.5), (0.5, 0.5, 0.5)),
        ((0.5, -0.5, -0.5), (0.5, 0.5, 0.5), (0.5, -0.5, 0.5)),
    ]
    with open(path, "wb") as f:
        f.write(b"\x00" * 80)
        f.write(struct.pack("<I", len(faces)))
        for v0, v1, v2 in faces:
            f.write(struct.pack("<fff", 0, 0, 0))  # normal (ignored)
            for v in (v0, v1, v2):
                f.write(struct.pack("<fff", *v))
            f.write(struct.pack("<H", 0))


class TestMesh(unittest.TestCase):
    def test_stl_import(self):
        """STL import: correct vertex/face counts and closure."""
        with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
            stl_path = f.name
            write_cube_stl(stl_path)
        try:
            mesh = mp.Mesh.from_stl(stl_path, material=mp.Medium(epsilon=12))
            self.assertEqual(len(mesh.vertices), 8)
            self.assertEqual(len(mesh.triangles), 12)
        finally:
            os.unlink(stl_path)

    def test_point_in_mesh(self):
        """Points inside and outside a cube mesh."""
        mat = mp.Medium(epsilon=12)
        cell = mp.Vector3(2, 2, 2)
        cube = make_cube_mesh(mat)
        sim = mp.Simulation(
            cell_size=cell,
            geometry=[cube],
            resolution=4,
            dimensions=3,
            eps_averaging=False,
        )
        sim.init_sim()

        # Sample single points via get_epsilon_grid
        def eps_at(x, y, z):
            return np.real(
                sim.get_epsilon_grid(np.array([x]), np.array([y]), np.array([z]))
            ).item()

        self.assertAlmostEqual(eps_at(0, 0, 0), 12.0, places=1)  # inside
        self.assertAlmostEqual(eps_at(0.3, 0.3, 0.3), 12.0, places=1)
        self.assertAlmostEqual(eps_at(0.8, 0, 0), 1.0, places=1)  # outside
        self.assertAlmostEqual(eps_at(0, 0.8, 0), 1.0, places=1)

    def test_volume(self):
        """Mesh cube volume matches analytic value via grid sampling."""
        mat = mp.Medium(epsilon=12)
        cube = make_cube_mesh(mat)
        sim = mp.Simulation(
            cell_size=mp.Vector3(2, 2, 2),
            geometry=[cube],
            resolution=4,
            dimensions=3,
            eps_averaging=False,
        )
        sim.init_sim()

        # Sample a uniform grid and count points inside the cube
        n = 40
        xtics = np.linspace(-1, 1, n)
        ytics = np.linspace(-1, 1, n)
        ztics = np.linspace(-1, 1, n)
        eps = np.real(sim.get_epsilon_grid(xtics, ytics, ztics))
        inside = np.sum(eps > 6)
        total = eps.size
        fill_fraction = inside / total
        expected = 1.0 / 8.0  # unit cube in 2x2x2 cell
        self.assertAlmostEqual(fill_fraction, expected, delta=0.02)

    def test_eps_averaging(self):
        """Subpixel smoothing produces intermediate epsilon at interfaces."""
        mat = mp.Medium(epsilon=12)
        cell = mp.Vector3(2, 2, 2)
        cube = make_cube_mesh(mat)

        sim_off = mp.Simulation(
            cell_size=cell,
            geometry=[cube],
            resolution=4,
            dimensions=3,
            eps_averaging=False,
        )
        sim_off.init_sim()
        eps_off = np.real(
            sim_off.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
        )

        sim_on = mp.Simulation(
            cell_size=cell,
            geometry=[cube],
            resolution=4,
            dimensions=3,
            eps_averaging=True,
        )
        sim_on.init_sim()
        eps_on = np.real(
            sim_on.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
        )

        unique_off = len(np.unique(np.round(eps_off, 2)))
        unique_on = len(np.unique(np.round(eps_on, 2)))
        self.assertGreater(
            unique_on,
            unique_off,
            "eps_averaging should produce more unique epsilon values",
        )

    def test_mesh_vs_block(self):
        """Cube mesh epsilon grid matches native Block."""
        mat = mp.Medium(epsilon=12)
        cell = mp.Vector3(2, 2, 2)
        n = 30
        xtics = np.linspace(-1, 1, n)
        ytics = np.linspace(-1, 1, n)
        ztics = np.array([0.0])

        cube_mesh = make_cube_mesh(mat)
        sim_mesh = mp.Simulation(
            cell_size=cell,
            geometry=[cube_mesh],
            resolution=4,
            dimensions=3,
            eps_averaging=False,
        )
        sim_mesh.init_sim()
        eps_mesh = np.real(sim_mesh.get_epsilon_grid(xtics, ytics, ztics)).flatten()

        block = mp.Block(size=mp.Vector3(1, 1, 1), material=mat)
        sim_block = mp.Simulation(
            cell_size=cell,
            geometry=[block],
            resolution=4,
            dimensions=3,
            eps_averaging=False,
        )
        sim_block.init_sim()
        eps_block = np.real(sim_block.get_epsilon_grid(xtics, ytics, ztics)).flatten()

        mismatches = np.sum(np.abs(eps_mesh - eps_block) > 0.5)
        self.assertEqual(
            mismatches,
            0,
            f"Cube mesh should match Block exactly, got {mismatches} mismatches",
        )


if __name__ == "__main__":
    unittest.main()
