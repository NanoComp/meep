# visualization.py - Tests the visualization module. Checks 2D
# plotting of a waveguide with several sources, monitors, and
# boundary conditions. Checks for subdomain plots.
#
# Also tests the animation run function, mp4 output, jshtml output, and git output.
import os
import unittest
from subprocess import call

import matplotlib
import numpy as np

import meep as mp

# Make sure we have matplotlib installed

matplotlib.use("agg")  # Set backend for consistency and to pull pixels quickly
import io

from matplotlib import pyplot as plt


def hash_figure(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="raw")
    buf.seek(0)
    data = np.frombuffer(buf.getvalue(), dtype=np.uint8)
    return np.sum((data > np.mean(data)) + data)


def setup_sim(zDim=0):
    cell = mp.Vector3(16, 8, zDim)

    # A simple waveguide
    geometry = [
        mp.Block(
            mp.Vector3(mp.inf, 1, 1),
            center=mp.Vector3(),
            material=mp.Medium(epsilon=12),
        )
    ]

    # Add point sources
    sources = [
        mp.Source(
            mp.ContinuousSource(frequency=0.15),
            component=mp.Ez,
            center=mp.Vector3(-5, 0),
            size=mp.Vector3(0, 0, 2),
        ),
        mp.Source(
            mp.ContinuousSource(frequency=0.15),
            component=mp.Ez,
            center=mp.Vector3(0, 2),
            size=mp.Vector3(0, 0, 2),
        ),
        mp.Source(
            mp.ContinuousSource(frequency=0.15),
            component=mp.Ez,
            center=mp.Vector3(-1, 1),
            size=mp.Vector3(0, 0, 2),
        ),
        mp.Source(
            mp.ContinuousSource(frequency=0.15),
            component=mp.Ez,
            center=mp.Vector3(-2, -2, 1),
            size=mp.Vector3(0, 0, 0),
        ),
    ]

    # Add line sources
    sources += [
        mp.Source(
            mp.ContinuousSource(frequency=0.15),
            component=mp.Ez,
            size=mp.Vector3(0, 2, 2),
            center=mp.Vector3(-6, 0),
        ),
        mp.Source(
            mp.ContinuousSource(frequency=0.15),
            component=mp.Ez,
            size=mp.Vector3(0, 2, 2),
            center=mp.Vector3(0, 1),
        ),
    ]

    # Add plane sources
    sources += [
        mp.Source(
            mp.ContinuousSource(frequency=0.15),
            component=mp.Ez,
            size=mp.Vector3(2, 2, 2),
            center=mp.Vector3(-3, 0),
        ),
        mp.Source(
            mp.ContinuousSource(frequency=0.15),
            component=mp.Ez,
            size=mp.Vector3(2, 2, 2),
            center=mp.Vector3(0, -2),
        ),
    ]

    # Different pml layers
    pml_layers = [
        mp.PML(2.0, mp.X),
        mp.PML(1.0, mp.Y, mp.Low),
        mp.PML(1.5, mp.Y, mp.High),
    ]
    if zDim > 0:
        pml_layers += [mp.PML(1.5, mp.Z)]

    resolution = 10

    sim = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        geometry=geometry,
        sources=sources,
        resolution=resolution,
    )
    # Line monitor
    sim.add_flux(
        1,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(5, 0, 0), size=mp.Vector3(0, 4, 4), direction=mp.X
        ),
    )

    # Plane monitor
    sim.add_flux(
        1,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(2, 0, 0), size=mp.Vector3(4, 4, 4), direction=mp.X
        ),
    )

    return sim


def view_sim():
    sim = setup_sim(8)
    xy0 = mp.Volume(
        center=mp.Vector3(0, 0, 0), size=mp.Vector3(sim.cell_size.x, sim.cell_size.y, 0)
    )
    xy1 = mp.Volume(
        center=mp.Vector3(0, 0, 1), size=mp.Vector3(sim.cell_size.x, sim.cell_size.y, 0)
    )
    yz0 = mp.Volume(
        center=mp.Vector3(0, 0, 0), size=mp.Vector3(0, sim.cell_size.y, sim.cell_size.z)
    )
    yz1 = mp.Volume(
        center=mp.Vector3(1, 0, 0), size=mp.Vector3(0, sim.cell_size.y, sim.cell_size.z)
    )
    xz0 = mp.Volume(
        center=mp.Vector3(0, 0, 0), size=mp.Vector3(sim.cell_size.x, 0, sim.cell_size.z)
    )
    xz1 = mp.Volume(
        center=mp.Vector3(0, 1, 0), size=mp.Vector3(sim.cell_size.x, 0, sim.cell_size.z)
    )
    vols = [xy0, xy1, yz0, yz1, xz0, xz1]
    titles = ["xy0", "xy1", "yz0", "yz1", "xz0", "xz1"]
    xlabel = ["x", "x", "y", "y", "x", "x"]
    ylabel = ["y", "y", "z", "z", "z", "z"]
    for k in range(len(vols)):
        ax = plt.subplot(2, 3, k + 1)
        sim.plot2D(ax=ax, output_plane=vols[k])
        ax.set_xlabel(xlabel[k])
        ax.set_ylabel(ylabel[k])
        ax.set_title(titles[k])
    plt.tight_layout()
    plt.show()


class TestVisualization(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_dir = mp.make_output_directory()

    @classmethod
    def tearDownClass(cls):
        mp.delete_directory(cls.temp_dir)

    def test_plot2D(self):
        # Check plotting of geometry with several sources, monitors, and PMLs
        f = plt.figure()
        ax = f.gca()
        sim = setup_sim()
        ax = sim.plot2D(ax=ax)
        if mp.am_master():
            hash_figure(f)
            # self.assertAlmostEqual(hash_figure(f),10231488)

        # Check plotting of fields after timestepping
        f = plt.figure()
        ax = f.gca()
        sim.run(until=200)
        ax = sim.plot2D(ax=ax, fields=mp.Ez)
        if mp.am_master():
            hash_figure(f)
            # self.assertAlmostEqual(hash_figure(f),79786722)

        # Check output_plane feature
        f = plt.figure()
        ax = f.gca()
        vol = mp.Volume(center=mp.Vector3(), size=mp.Vector3(2, 2))
        ax = sim.plot2D(ax=ax, fields=mp.Ez, output_plane=vol)
        if mp.am_master():
            hash_figure(f)
            # self.assertAlmostEqual(hash_figure(f),68926258)

    @unittest.skipIf(call(["which", "ffmpeg"]) != 0, "ffmpeg is not installed")
    def test_animation_output(self):
        # ------------------------- #
        # Check over 2D domain
        # ------------------------- #

        sim = setup_sim()  # generate 2D simulation

        Animate = mp.Animate2D(
            sim=sim, fields=mp.Ez, realtime=False, normalize=False
        )  # Check without normalization
        Animate_norm = mp.Animate2D(
            sim=sim, fields=mp.Ez, realtime=False, normalize=True
        )  # Check with normalization

        # test both animation objects during same run
        sim.run(mp.at_every(1, Animate), mp.at_every(1, Animate_norm), until=5)

        # Test outputs
        Animate.to_mp4(
            5, os.path.join(self.temp_dir, "test_2D.mp4")
        )  # Check mp4 output
        Animate.to_gif(
            150, os.path.join(self.temp_dir, "test_2D.gif")
        )  # Check gif output
        Animate.to_jshtml(10)  # Check jshtml output
        Animate_norm.to_mp4(
            5, os.path.join(self.temp_dir, "test_2D_norm.mp4")
        )  # Check mp4 output
        Animate_norm.to_gif(
            150, os.path.join(self.temp_dir, "test_2D_norm.gif")
        )  # Check gif output
        Animate_norm.to_jshtml(10)  # Check jshtml output

        # ------------------------- #
        # Check over 3D domain
        # ------------------------- #
        sim = setup_sim(5)  # generate 2D simulation

        Animate_xy = mp.Animate2D(
            sim=sim, fields=mp.Ey, realtime=False, normalize=True
        )  # Check without normalization
        Animate_xz = mp.Animate2D(
            sim=sim, fields=mp.Ey, realtime=False, normalize=True
        )  # Check with normalization

        # test both animation objects during same run
        sim.run(
            mp.at_every(
                1,
                mp.in_volume(
                    mp.Volume(
                        center=mp.Vector3(),
                        size=mp.Vector3(sim.cell_size.x, sim.cell_size.y),
                    ),
                    Animate_xy,
                ),
            ),
            mp.at_every(
                1,
                mp.in_volume(
                    mp.Volume(
                        center=mp.Vector3(),
                        size=mp.Vector3(sim.cell_size.x, 0, sim.cell_size.z),
                    ),
                    Animate_xz,
                ),
            ),
            until=5,
        )

        # Test outputs
        Animate_xy.to_mp4(
            5, os.path.join(self.temp_dir, "test_3D_xy.mp4")
        )  # Check mp4 output
        Animate_xy.to_gif(
            150, os.path.join(self.temp_dir, "test_3D_xy.gif")
        )  # Check gif output
        Animate_xy.to_jshtml(10)  # Check jshtml output
        Animate_xz.to_mp4(
            5, os.path.join(self.temp_dir, "test_3D_xz.mp4")
        )  # Check mp4 output
        Animate_xz.to_gif(
            150, os.path.join(self.temp_dir, "test_3D_xz.gif")
        )  # Check gif output
        Animate_xz.to_jshtml(10)  # Check jshtml output

    """
    Travis does not play well with Mayavi
    def test_3D_mayavi(self):
        sim = setup_sim(4)
        sim.plot3D()
    """


if __name__ == "__main__":
    unittest.main()
