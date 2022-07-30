from collections import namedtuple
import warnings

import numpy as np

import meep as mp
from meep.geom import Vector3, init_do_averaging
from meep.source import EigenModeSource, check_positive
from meep.simulation import Simulation, Volume

## Typing imports
from matplotlib.axes import Axes
from typing import Callable, Union, Any, Iterable

# ------------------------------------------------------- #
# Visualization
# ------------------------------------------------------- #
# Contains all necessary visualization routines for use with
# pymeep and pympb.

# ------------------------------------------------------- #
# Functions used to define the default plotting parameters
# for the different plotting routines.

default_source_parameters = {
    "color": "r",
    "edgecolor": "r",
    "facecolor": "none",
    "hatch": "/",
    "linewidth": 2,
}

default_monitor_parameters = {
    "color": "b",
    "edgecolor": "b",
    "facecolor": "none",
    "hatch": "/",
    "linewidth": 2,
}

default_field_parameters = {
    "interpolation": "spline36",
    "cmap": "RdBu",
    "alpha": 0.8,
    "post_process": np.real,
}

default_eps_parameters = {
    "interpolation": "spline36",
    "cmap": "binary",
    "alpha": 1.0,
    "contour": False,
    "contour_linewidth": 1,
    "frequency": None,
    "resolution": None,
}

default_boundary_parameters = {
    "color": "g",
    "edgecolor": "g",
    "facecolor": "none",
    "hatch": "/",
}

default_volume_parameters = {
    "alpha": 1.0,
    "color": "k",
    "linestyle": "-",
    "linewidth": 1,
    "marker": ".",
    "edgecolor": "k",
    "facecolor": "none",
    "hatch": "/",
}

default_label_parameters = {"label_color": "r", "offset": 20, "label_alpha": 0.3}

# Used to remove the elements of a dictionary (dict_to_filter) that
# don't correspond to the keyword arguments of a particular
# function (func_with_kwargs.)
# Adapted from https://stackoverflow.com/questions/26515595/how-does-one-ignore-unexpected-keyword-arguments-passed-to-a-function/44052550
def filter_dict(dict_to_filter: dict, func_with_kwargs: Callable) -> dict:
    import inspect

    filter_keys = []
    try:
        # Python3 ...
        sig = inspect.signature(func_with_kwargs)
        filter_keys = [param.name for param in sig.parameters.values()]
    except:
        # Python2 ...
        filter_keys = inspect.getargspec(func_with_kwargs)[0]

    filtered_dict = {
        filter_key: dict_to_filter[filter_key]
        for filter_key in filter_keys
        if filter_key in dict_to_filter
    }
    return filtered_dict


# ------------------------------------------------------- #
# Routines to add legends to plot


def place_label(
    ax: Axes, label_text: str, x, y, centerx, centery, label_parameters: dict = None
) -> Axes:

    if label_parameters is None:
        label_parameters = default_label_parameters
    else:
        label_parameters = dict(default_label_parameters, **label_parameters)

    offset = label_parameters["offset"]
    alpha = label_parameters["label_alpha"]
    color = label_parameters["label_color"]

    if x > centerx:
        xtext = -offset
    else:
        xtext = offset
    if y > centery:
        ytext = -offset
    else:
        ytext = offset

    ax.annotate(
        label_text,
        xy=(x, y),
        xytext=(xtext, ytext),
        textcoords="offset points",
        ha="center",
        va="bottom",
        bbox=dict(boxstyle="round,pad=0.2", fc=color, alpha=alpha),
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0.5", color=color),
    )
    return ax


# ------------------------------------------------------- #
# Helper functions used to plot volumes on a 2D plane

# Returns the intersection points of two Volumes.
# Volumes must be a line, plane, or rectangular prism
# (since they are volume objects)
def intersect_volume_volume(volume1: Volume, volume2: Volume) -> list:
    # volume1 ............... [volume]
    # volume2 ............... [volume]

    # Represent the volumes by an "upper" and "lower" coordinate
    U1 = [
        volume1.center.x + volume1.size.x / 2,
        volume1.center.y + volume1.size.y / 2,
        volume1.center.z + volume1.size.z / 2,
    ]
    L1 = [
        volume1.center.x - volume1.size.x / 2,
        volume1.center.y - volume1.size.y / 2,
        volume1.center.z - volume1.size.z / 2,
    ]

    U2 = [
        volume2.center.x + volume2.size.x / 2,
        volume2.center.y + volume2.size.y / 2,
        volume2.center.z + volume2.size.z / 2,
    ]
    L2 = [
        volume2.center.x - volume2.size.x / 2,
        volume2.center.y - volume2.size.y / 2,
        volume2.center.z - volume2.size.z / 2,
    ]

    # Evaluate intersection
    U = np.min([U1, U2], axis=0)
    L = np.max([L1, L2], axis=0)

    # For single points we have to check manually
    if np.all(U - L == 0):
        if (not volume1.pt_in_volume(Vector3(*U))) or (
            not volume2.pt_in_volume(Vector3(*U))
        ):
            return []

    # Check for two volumes that don't intersect
    if np.any(U - L < 0):
        return []

    # Pull all possible vertices
    vertices = []
    for x_vals in [L[0], U[0]]:
        for y_vals in [L[1], U[1]]:
            for z_vals in [L[2], U[2]]:
                vertices.append(Vector3(x_vals, y_vals, z_vals))

    # Remove any duplicate points caused by coplanar lines
    vertices = [
        vertices[i] for i, x in enumerate(vertices) if x not in vertices[i + 1 :]
    ]

    return vertices


# All of the 2D plotting routines need an output plane over which to plot.
# The user has many options to specify this output plane. They can pass
# the output_plane parameter, which is a 2D volume object. They can specify
# a volume using in_volume, which stores the volume as a C volume, not a Python
# volume. They can also do nothing and plot the XY plane through Z=0.
#
# Not only do we need to check for all of these possibilities, but we also need
# to check if the user accidentally specifies a plane that stretches beyond the
# simulation domain.
def get_2D_dimensions(sim: Simulation, output_plane: Volume) -> (Vector3, Vector3):
    # Pull correct plane from user
    if output_plane:
        plane_center, plane_size = (output_plane.center, output_plane.size)
    elif sim.output_volume:
        plane_center, plane_size = mp.get_center_and_size(sim.output_volume)
    else:
        if (sim.dimensions == mp.CYLINDRICAL) or sim.is_cylindrical:
            plane_center, plane_size = (
                sim.geometry_center + Vector3(sim.cell_size.x / 2),
                sim.cell_size,
            )
        else:
            plane_center, plane_size = (sim.geometry_center, sim.cell_size)
    plane_volume = Volume(center=plane_center, size=plane_size)

    if plane_size.x != 0 and plane_size.y != 0 and plane_size.z != 0:
        raise ValueError("Plane volume must be 2D (a plane).")
    if (sim.dimensions == mp.CYLINDRICAL) or sim.is_cylindrical:
        center = sim.geometry_center + Vector3(sim.cell_size.x / 2)
        check_volume = Volume(center=center, size=sim.cell_size)
    else:
        check_volume = Volume(center=sim.geometry_center, size=sim.cell_size)
    vertices = intersect_volume_volume(check_volume, plane_volume)

    if len(vertices) == 0:
        raise ValueError(
            "The specified user volume is completely outside of the simulation domain."
        )

    intersection_vol = Volume(vertices=vertices)

    if (intersection_vol.size != plane_volume.size) or (
        intersection_vol.center != plane_volume.center
    ):
        warnings.warn(
            "The specified user volume is larger than the simulation domain and has been truncated."
        )

    sim_center, sim_size = (intersection_vol.center, intersection_vol.size)
    return sim_center, sim_size


def box_vertices(
    box_center: Vector3, box_size: Vector3, is_cylindrical: bool = False
) -> (float, float, float, float, float, float):
    # in cylindrical coordinates, radial (R) axis
    # is in the range (0,R) rather than (-R/2,+R/2)
    # as in Cartesian coordinates.
    if is_cylindrical:
        xmin = 0
        xmax = box_size.x
    else:
        xmin = box_center.x - 0.5 * box_size.x
        xmax = box_center.x + 0.5 * box_size.x
    ymin = box_center.y - 0.5 * box_size.y
    ymax = box_center.y + 0.5 * box_size.y
    zmin = box_center.z - 0.5 * box_size.z
    zmax = box_center.z + 0.5 * box_size.z

    return xmin, xmax, ymin, ymax, zmin, zmax


# ------------------------------------------------------- #
# actual plotting routines
def plot_volume(
    sim: Simulation,
    ax: Axes,
    volume: Volume,
    output_plane: Volume = None,
    plotting_parameters: dict = None,
    label: str = None,
) -> Axes:
    import matplotlib.patches as patches
    from matplotlib import pyplot as plt

    # Set up the plotting parameters
    if plotting_parameters is None:
        plotting_parameters = default_volume_parameters
    else:
        plotting_parameters = dict(default_volume_parameters, **plotting_parameters)

    # Get domain measurements
    sim_center, sim_size = get_2D_dimensions(sim, output_plane)

    plane = Volume(center=sim_center, size=sim_size)

    size = volume.size
    center = volume.center

    xmin, xmax, ymin, ymax, zmin, zmax = box_vertices(center, size, sim.is_cylindrical)

    # Add labels if requested
    if label is not None and mp.am_master():
        if sim_size.x == 0:
            ax = place_label(
                ax,
                label,
                center.y,
                center.z,
                sim_center.y,
                sim_center.z,
                label_parameters=plotting_parameters,
            )
        elif sim_size.y == 0:
            ax = place_label(
                ax,
                label,
                center.x,
                center.z,
                sim_center.x,
                sim_center.z,
                label_parameters=plotting_parameters,
            )
        elif sim_size.z == 0:
            ax = place_label(
                ax,
                label,
                center.x,
                center.y,
                sim_center.x,
                sim_center.y,
                label_parameters=plotting_parameters,
            )

    # Intersect plane with volume
    intersection = intersect_volume_volume(volume, plane)

    # Sort the points in a counter clockwise manner to ensure convex polygon is formed
    def sort_points(xy):
        xy = np.squeeze(xy)
        xy_mean = np.mean(xy, axis=0)
        theta = np.arctan2(xy[:, 1] - xy_mean[1], xy[:, 0] - xy_mean[0])
        return xy[np.argsort(theta, axis=0), :]

    if mp.am_master():
        # Point volume
        if len(intersection) == 1:
            point_args = {
                key: value
                for key, value in plotting_parameters.items()
                if key in ["color", "marker", "alpha", "linewidth"]
            }
            if sim_size.y == 0:
                ax.scatter(center.x, center.z, **point_args)
                return ax
            elif sim_size.x == 0:
                ax.scatter(center.y, center.z, **point_args)
                return ax
            elif sim_size.z == 0:
                ax.scatter(center.x, center.y, **point_args)
                return ax
            else:
                return ax

        # Line volume
        elif len(intersection) == 2:
            line_args = {
                key: value
                for key, value in plotting_parameters.items()
                if key in ["color", "linestyle", "linewidth", "alpha"]
            }
            # Plot YZ
            if sim_size.x == 0:
                ax.plot(
                    [a.y for a in intersection],
                    [a.z for a in intersection],
                    **line_args
                )
                return ax
            # Plot XZ
            elif sim_size.y == 0:
                ax.plot(
                    [a.x for a in intersection],
                    [a.z for a in intersection],
                    **line_args
                )
                return ax
            # Plot XY
            elif sim_size.z == 0:
                ax.plot(
                    [a.x for a in intersection],
                    [a.y for a in intersection],
                    **line_args
                )
                return ax
            else:
                return ax

        # Planar volume
        elif len(intersection) > 2:
            planar_args = {
                key: value
                for key, value in plotting_parameters.items()
                if key in ["edgecolor", "linewidth", "facecolor", "hatch", "alpha"]
            }
            # Plot YZ
            if sim_size.x == 0:
                ax.add_patch(
                    patches.Polygon(
                        sort_points([[a.y, a.z] for a in intersection]), **planar_args
                    )
                )
                return ax
            # Plot XZ
            elif sim_size.y == 0:
                ax.add_patch(
                    patches.Polygon(
                        sort_points([[a.x, a.z] for a in intersection]), **planar_args
                    )
                )
                return ax
            # Plot XY
            elif sim_size.z == 0:
                ax.add_patch(
                    patches.Polygon(
                        sort_points([[a.x, a.y] for a in intersection]), **planar_args
                    )
                )
                return ax
            else:
                return ax
        else:
            return ax
    return ax


def plot_eps(
    sim: Simulation,
    ax: Axes,
    output_plane: Volume = None,
    eps_parameters: dict = None,
    frequency: float = None,
) -> Axes:
    # consolidate plotting parameters
    if eps_parameters is None:
        eps_parameters = default_eps_parameters
    else:
        eps_parameters = dict(default_eps_parameters, **eps_parameters)

    # Determine a frequency to plot all epsilon
    if frequency is not None:
        warnings.warn(
            "The frequency parameter of plot2D has been deprecated. "
            "Use the frequency key of the eps_parameters dictionary instead."
        )
        eps_parameters["frequency"] = frequency
    if eps_parameters["frequency"] is None:
        try:
            # Continuous sources
            eps_parameters["frequency"] = sim.sources[0].frequency
        except:
            try:
                # Gaussian sources
                eps_parameters["frequency"] = sim.sources[0].src.frequency
            except:
                try:
                    # Custom sources
                    eps_parameters["frequency"] = sim.sources[0].src.center_frequency
                except:
                    # No sources
                    eps_parameters["frequency"] = 0

    # Get domain measurements
    sim_center, sim_size = get_2D_dimensions(sim, output_plane)

    xmin, xmax, ymin, ymax, zmin, zmax = box_vertices(
        sim_center, sim_size, sim.is_cylindrical
    )

    if eps_parameters["resolution"]:
        grid_resolution = eps_parameters["resolution"]
    else:
        grid_resolution = sim.resolution

    Nx = int((xmax - xmin) * grid_resolution + 1)
    Ny = int((ymax - ymin) * grid_resolution + 1)
    Nz = int((zmax - zmin) * grid_resolution + 1)

    if sim_size.x == 0:
        # Plot y on x axis, z on y axis (YZ plane)
        extent = [ymin, ymax, zmin, zmax]
        xlabel = "Y"
        ylabel = "Z"
        xtics = np.array([sim_center.x])
        ytics = np.linspace(ymin, ymax, Ny)
        ztics = np.linspace(zmin, zmax, Nz)
    elif sim_size.y == 0:
        # Plot x on x axis, z on y axis (XZ plane)
        extent = [xmin, xmax, zmin, zmax]
        if (sim.dimensions == mp.CYLINDRICAL) or sim.is_cylindrical:
            xlabel = "R"
        else:
            xlabel = "X"
        ylabel = "Z"
        xtics = np.linspace(xmin, xmax, Nx)
        ytics = np.array([sim_center.y])
        ztics = np.linspace(zmin, zmax, Nz)
    elif sim_size.z == 0:
        # Plot x on x axis, y on y axis (XY plane)
        extent = [xmin, xmax, ymin, ymax]
        xlabel = "X"
        ylabel = "Y"
        xtics = np.linspace(xmin, xmax, Nx)
        ytics = np.linspace(ymin, ymax, Ny)
        ztics = np.array([sim_center.z])
    else:
        raise ValueError("A 2D plane has not been specified...")

    eps_data = np.rot90(
        np.real(sim.get_epsilon_grid(xtics, ytics, ztics, eps_parameters["frequency"]))
    )

    if mp.am_master():
        if eps_parameters["contour"]:
            ax.contour(
                eps_data,
                0,
                levels=np.unique(eps_data),
                colors="black",
                origin="upper",
                extent=extent,
                linewidths=eps_parameters["contour_linewidth"],
            )
        else:
            ax.imshow(eps_data, extent=extent, **filter_dict(eps_parameters, ax.imshow))
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    return ax


def plot_boundaries(
    sim: Simulation,
    ax: Axes,
    output_plane: Volume = None,
    boundary_parameters: dict = None,
) -> Axes:
    # consolidate plotting parameters
    if boundary_parameters is None:
        boundary_parameters = default_boundary_parameters
    else:
        boundary_parameters = dict(default_boundary_parameters, **boundary_parameters)

    def get_boundary_volumes(thickness: float, direction: float, side) -> Volume:
        thickness = boundary.thickness

        xmin, xmax, ymin, ymax, zmin, zmax = box_vertices(
            sim.geometry_center, sim.cell_size, sim.is_cylindrical
        )

        if direction == mp.X and side == mp.Low:
            return Volume(
                center=Vector3(
                    xmin + 0.5 * thickness, sim.geometry_center.y, sim.geometry_center.z
                ),
                size=Vector3(thickness, sim.cell_size.y, sim.cell_size.z),
            )
        elif (direction == mp.X and side == mp.High) or direction == mp.R:
            return Volume(
                center=Vector3(
                    xmax - 0.5 * thickness, sim.geometry_center.y, sim.geometry_center.z
                ),
                size=Vector3(thickness, sim.cell_size.y, sim.cell_size.z),
            )
        elif direction == mp.Y and side == mp.Low:
            return Volume(
                center=Vector3(
                    sim.geometry_center.x, ymin + 0.5 * thickness, sim.geometry_center.z
                ),
                size=Vector3(sim.cell_size.x, thickness, sim.cell_size.z),
            )
        elif direction == mp.Y and side == mp.High:
            return Volume(
                center=Vector3(
                    sim.geometry_center.x, ymax - 0.5 * thickness, sim.geometry_center.z
                ),
                size=Vector3(sim.cell_size.x, thickness, sim.cell_size.z),
            )
        elif direction == mp.Z and side == mp.Low:
            xcen = sim.geometry_center.x
            if sim.is_cylindrical:
                xcen += 0.5 * sim.cell_size.x
            return Volume(
                center=Vector3(xcen, sim.geometry_center.y, zmin + 0.5 * thickness),
                size=Vector3(sim.cell_size.x, sim.cell_size.y, thickness),
            )
        elif direction == mp.Z and side == mp.High:
            xcen = sim.geometry_center.x
            if sim.is_cylindrical:
                xcen += 0.5 * sim.cell_size.x
            return Volume(
                center=Vector3(xcen, sim.geometry_center.y, zmax - 0.5 * thickness),
                size=Vector3(sim.cell_size.x, sim.cell_size.y, thickness),
            )
        else:
            raise ValueError("Invalid boundary type")

    import itertools

    for boundary in sim.boundary_layers:
        # boundary on all four sides
        if boundary.direction == mp.ALL and boundary.side == mp.ALL:
            if sim.dimensions == 1:
                dims = [mp.X]
            elif sim.dimensions == mp.CYLINDRICAL or sim.is_cylindrical:
                dims = [mp.X, mp.Z]
            elif sim.dimensions == 2:
                dims = [mp.X, mp.Y]
            elif sim.dimensions == 3:
                dims = [mp.X, mp.Y, mp.Z]
            else:
                raise ValueError("Invalid simulation dimensions")
            for permutation in itertools.product(dims, [mp.Low, mp.High]):
                if ((permutation[0] == mp.X) and (permutation[1] == mp.Low)) and (
                    sim.dimensions == mp.CYLINDRICAL or sim.is_cylindrical
                ):
                    continue
                vol = get_boundary_volumes(boundary.thickness, *permutation)
                ax = plot_volume(
                    sim, ax, vol, output_plane, plotting_parameters=boundary_parameters
                )
        # boundary on only two of four sides
        elif boundary.side == mp.ALL:
            for side in [mp.Low, mp.High]:
                if ((boundary.direction == mp.X) and (side == mp.Low)) and (
                    sim.dimensions == mp.CYLINDRICAL or sim.is_cylindrical
                ):
                    continue
                vol = get_boundary_volumes(boundary.thickness, boundary.direction, side)
                ax = plot_volume(
                    sim, ax, vol, output_plane, plotting_parameters=boundary_parameters
                )
        # boundary on just one side
        else:
            if ((boundary.direction == mp.X) and (boundary.side == mp.Low)) and (
                sim.dimensions == mp.CYLINDRICAL or sim.is_cylindrical
            ):
                continue
            vol = get_boundary_volumes(
                boundary.thickness, boundary.direction, boundary.side
            )
            ax = plot_volume(
                sim, ax, vol, output_plane, plotting_parameters=boundary_parameters
            )
    return ax


def plot_sources(
    sim: Simulation,
    ax: Axes,
    output_plane: Volume = None,
    labels: bool = False,
    source_parameters: dict = None,
):

    # consolidate plotting parameters
    if source_parameters is None:
        source_parameters = default_source_parameters
    else:
        source_parameters = dict(default_source_parameters, **source_parameters)

    label = "source" if labels else None

    for src in sim.sources:
        vol = Volume(center=src.center, size=src.size)
        ax = plot_volume(
            sim,
            ax,
            vol,
            output_plane,
            plotting_parameters=source_parameters,
            label=label,
        )
    return ax


def plot_monitors(
    sim: Simulation,
    ax: Axes,
    output_plane: Volume = None,
    labels: bool = False,
    monitor_parameters: dict = None,
) -> Axes:
    # consolidate plotting parameters
    if monitor_parameters is None:
        monitor_parameters = default_monitor_parameters
    else:
        monitor_parameters = dict(default_monitor_parameters, **monitor_parameters)

    label = "monitor" if labels else None

    for mon in sim.dft_objects:
        for reg in mon.regions:
            vol = Volume(center=reg.center, size=reg.size)
            ax = plot_volume(
                sim,
                ax,
                vol,
                output_plane,
                plotting_parameters=monitor_parameters,
                label=label,
            )
    return ax


def plot_fields(
    sim: Simulation,
    ax: Axes = None,
    fields=None,
    output_plane: Volume = None,
    field_parameters: dict = None,
) -> Union[Axes, Any]:
    if not sim._is_initialized:
        sim.init_sim()

    if fields is None:
        return ax

    if field_parameters is None:
        field_parameters = default_field_parameters
    else:
        field_parameters = dict(default_field_parameters, **field_parameters)

    # user specifies a field component
    if fields in [
        mp.Ex,
        mp.Ey,
        mp.Ez,
        mp.Er,
        mp.Ep,
        mp.Dx,
        mp.Dy,
        mp.Dz,
        mp.Hx,
        mp.Hy,
        mp.Hz,
    ]:
        # Get domain measurements
        sim_center, sim_size = get_2D_dimensions(sim, output_plane)

        xmin, xmax, ymin, ymax, zmin, zmax = box_vertices(
            sim_center, sim_size, sim.is_cylindrical
        )

        if sim_size.x == 0:
            # Plot y on x axis, z on y axis (YZ plane)
            extent = [ymin, ymax, zmin, zmax]
            xlabel = "Y"
            ylabel = "Z"
        elif sim_size.y == 0:
            # Plot x on x axis, z on y axis (XZ plane)
            extent = [xmin, xmax, zmin, zmax]
            if (sim.dimensions == mp.CYLINDRICAL) or sim.is_cylindrical:
                xlabel = "R"
            else:
                xlabel = "X"
            ylabel = "Z"
        elif sim_size.z == 0:
            # Plot x on x axis, y on y axis (XY plane)
            extent = [xmin, xmax, ymin, ymax]
            xlabel = "X"
            ylabel = "Y"
        fields = sim.get_array(center=sim_center, size=sim_size, component=fields)
    else:
        raise ValueError("Please specify a valid field component (mp.Ex, mp.Ey, ...")

    fields = field_parameters["post_process"](fields)
    if (sim.dimensions == mp.CYLINDRICAL) or sim.is_cylindrical:
        fields = np.flipud(fields)
    else:
        fields = np.rot90(fields)

    # Either plot the field, or return the array
    if ax:
        if mp.am_master():
            ax.imshow(fields, extent=extent, **filter_dict(field_parameters, ax.imshow))
        return ax
    else:
        return fields
    return ax


def plot2D(
    sim: Simulation,
    ax: Axes = None,
    output_plane: Volume = None,
    fields=None,
    labels: bool = False,
    eps_parameters: dict = None,
    boundary_parameters: dict = None,
    source_parameters: dict = None,
    monitor_parameters: dict = None,
    field_parameters: dict = None,
    frequency: float = None,
    plot_eps_flag: bool = True,
    plot_sources_flag: bool = True,
    plot_monitors_flag: bool = True,
    plot_boundaries_flag: bool = True,
) -> Axes:

    # Ensure a figure axis exists
    if ax is None and mp.am_master():
        from matplotlib import pyplot as plt

        ax = plt.gca()

    # validate the output plane to ensure proper 2D coordinates
    sim_center, sim_size = get_2D_dimensions(sim, output_plane)
    output_plane = Volume(center=sim_center, size=sim_size)

    # Plot geometry
    if plot_eps_flag:
        ax = plot_eps(
            sim,
            ax,
            output_plane=output_plane,
            eps_parameters=eps_parameters,
            frequency=frequency,
        )

    # Plot boundaries
    if plot_boundaries_flag:
        ax = plot_boundaries(
            sim, ax, output_plane=output_plane, boundary_parameters=boundary_parameters
        )

    # Plot sources
    if plot_sources_flag:
        ax = plot_sources(
            sim,
            ax,
            output_plane=output_plane,
            labels=labels,
            source_parameters=source_parameters,
        )

    # Plot monitors
    if plot_monitors_flag:
        ax = plot_monitors(
            sim,
            ax,
            output_plane=output_plane,
            labels=labels,
            monitor_parameters=monitor_parameters,
        )

    # Plot fields
    if fields is not None:
        ax = plot_fields(
            sim,
            ax,
            fields,
            output_plane=output_plane,
            field_parameters=field_parameters,
        )

    return ax


def plot3D(sim: Simulation):
    from mayavi import mlab

    if sim.dimensions < 3:
        raise ValueError("Simulation must have 3 dimensions to visualize 3D")

    xmin, xmax, ymin, ymax, zmin, zmax = box_vertices(
        sim.geometry_center, sim.cell_size
    )

    Nx = int(sim.cell_size.x * sim.resolution) + 1
    Ny = int(sim.cell_size.y * sim.resolution) + 1
    Nz = int(sim.cell_size.z * sim.resolution) + 1

    xtics = np.linspace(xmin, xmax, Nx)
    ytics = np.linspace(ymin, ymax, Ny)
    ztics = np.linspace(zmin, zmax, Nz)

    eps_data = sim.get_epsilon_grid(xtics, ytics, ztics)
    s = mlab.contour3d(eps_data, colormap="YlGnBu")
    return s


def visualize_chunks(sim: Simulation):
    if sim.structure is None:
        sim.init_sim()

    import matplotlib.pyplot as plt
    import matplotlib.cm
    import matplotlib.colors

    if sim.structure.gv.dim == 2:
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    else:
        from matplotlib.collections import PolyCollection

    vols = sim.structure.get_chunk_volumes()
    owners = sim.structure.get_chunk_owners()

    def plot_box(box, proc, fig, ax):
        if sim.structure.gv.dim == 2:
            low = Vector3(box.low.x, box.low.y, box.low.z)
            high = Vector3(box.high.x, box.high.y, box.high.z)
            points = [low, high]

            x_len = Vector3(high.x) - Vector3(low.x)
            y_len = Vector3(y=high.y) - Vector3(y=low.y)
            xy_len = Vector3(high.x, high.y) - Vector3(low.x, low.y)

            points += [low + x_len]
            points += [low + y_len]
            points += [low + xy_len]
            points += [high - x_len]
            points += [high - y_len]
            points += [high - xy_len]
            points = np.array([np.array(v) for v in points])

            edges = [
                [points[0], points[2], points[4], points[3]],
                [points[1], points[5], points[7], points[6]],
                [points[0], points[3], points[5], points[7]],
                [points[1], points[4], points[2], points[6]],
                [points[3], points[4], points[1], points[5]],
                [points[0], points[7], points[6], points[2]],
            ]

            faces = Poly3DCollection(edges, linewidths=1, edgecolors="k")
            color_with_alpha = matplotlib.colors.to_rgba(chunk_colors[proc], alpha=0.2)
            faces.set_facecolor(color_with_alpha)
            ax.add_collection3d(faces)

            # Plot the points themselves to force the scaling of the axes
            ax.scatter(points[:, 0], points[:, 1], points[:, 2], s=0)
        else:
            low = Vector3(box.low.x, box.low.y)
            high = Vector3(box.high.x, box.high.y)
            points = [low, high]

            x_len = Vector3(high.x) - Vector3(low.x)
            y_len = Vector3(y=high.y) - Vector3(y=low.y)

            points += [low + x_len]
            points += [low + y_len]
            points = np.array([np.array(v)[:-1] for v in points])

            edges = [[points[0], points[2], points[1], points[3]]]

            faces = PolyCollection(edges, linewidths=1, edgecolors="k")
            color_with_alpha = matplotlib.colors.to_rgba(chunk_colors[proc])
            faces.set_facecolor(color_with_alpha)
            ax.add_collection(faces)

            # Plot the points themselves to force the scaling of the axes
            ax.scatter(points[:, 0], points[:, 1], s=0)

    if mp.am_master():
        fig = plt.figure()
        ax = fig.add_subplot(
            111, projection="3d" if sim.structure.gv.dim == 2 else None
        )
        chunk_colors = matplotlib.cm.rainbow(np.linspace(0, 1, mp.count_processors()))

        for i, v in enumerate(vols):
            plot_box(mp.gv2box(v.surroundings()), owners[i], fig, ax)

        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect("equal")

        cell_box = mp.gv2box(sim.structure.gv.surroundings())
        if sim.structure.gv.dim == 2:
            ax.set_xlim3d(left=cell_box.low.x, right=cell_box.high.x)
            ax.set_ylim3d(bottom=cell_box.low.y, top=cell_box.high.y)
            ax.set_zlim3d(bottom=cell_box.low.z, top=cell_box.high.z)
            ax.set_zlabel("z")
        else:
            ax.set_xlim(left=cell_box.low.x, right=cell_box.high.x)
            ax.set_ylim(bottom=cell_box.low.y, top=cell_box.high.y)

        plt.tight_layout()
        plt.show()


# ------------------------------------------------------- #
# Animate2D
# ------------------------------------------------------- #
# An extensive run function used to visualize the fields
# of a 2D simulation after every specified time step.
# ------------------------------------------------------- #
# Required arguments
# sim ................. [Simulation object]
# fields .............. [mp.Ex, mp.Ey, ..., mp. Hz]
# ------------------------------------------------------- #
# Optional arguments
# f ................... [matplotlib figure object]
# realtime ............ [bool] Update plot in each step
# normalize ........... [bool] saves fields to normalize
#                       after simulation ends.
# plot_modifiers ...... [list] additional functions to
#                       modify plot
# customization_args .. [dict] other customization args
#                       to pass to plot2D()


class Animate2D:
    """
    A class used to record the fields during timestepping (i.e., a [`run`](#run-functions)
    function). The object is initialized prior to timestepping by specifying the
    simulation object and the field component. The object can then be passed to any
    [step-function modifier](#step-function-modifiers). For example, one can record the
    $E_z$ fields at every one time unit using:

    ```py
    animate = mp.Animate2D(sim,
                           fields=mp.Ez,
                           realtime=True,
                           field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'none'},
                           boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3})

    sim.run(mp.at_every(1,animate),until=25)
    ```

    By default, the object saves each frame as a PNG image into memory (not disk). This is
    typically more memory efficient than storing the actual fields. If the user sets the
    `normalize` argument, then the object will save the actual field information as a
    NumPy array to be normalized for post processing. The fields of a figure can also be
    updated in realtime by setting the `realtime` flag. This does not work for
    IPython/Jupyter notebooks, however.

    Once the simulation is run, the animation can be output as an interactive JSHTML
    object, an mp4, or a GIF.

    Multiple `Animate2D` objects can be initialized and passed to the run function to
    track different volume locations (using `mp.in_volume`) or field components.
    """

    def __init__(
        self,
        sim,
        fields,
        f=None,
        realtime=False,
        normalize=False,
        plot_modifiers=None,
        **customization_args
    ):
        """
        Construct an `Animate2D` object.

        + **`sim`** — Simulation object.

        + **`fields`** — Field component to record at each time instant.

        + **`f=None`** — Optional `matplotlib` figure object that the routine will update
          on each call. If not supplied, then a new one will be created upon
          initialization.

        + **`realtime=False`** — Whether or not to update a figure window in realtime as
          the simulation progresses. Disabled by default. Not compatible with
          IPython/Jupyter notebooks.

        + **`normalize=False`** — Records fields at each time step in memory in a NumPy
          array and then normalizes the result by dividing by the maximum field value at a
          single point in the cell over all the time snapshots.

        + **`plot_modifiers=None`** — A list of functions that can modify the figure's
          `axis` object. Each function modifier accepts a single argument, an `axis`
          object, and must return that same axis object. The following modifier changes
          the `xlabel`:

        ```py
          def mod1(ax):
              ax.set_xlabel('Testing')
              return ax

          plot_modifiers = [mod1]
        ```

        + **`**customization_args`** — Customization keyword arguments passed to
          `plot2D()` (i.e. `labels`, `eps_parameters`, `boundary_parameters`, etc.)
        """

        self.fields = fields

        if f:
            self.f = f
            self.ax = self.f.gca()
        elif mp.am_master():
            from matplotlib import pyplot as plt

            self.f = plt.figure()
            self.ax = self.f.gca()
        else:
            self.f = None
            self.ax = None

        self.realtime = realtime
        self.normalize = normalize
        self.plot_modifiers = plot_modifiers
        self.customization_args = customization_args

        self.cumulative_fields = []
        self._saved_frames = []

        self.frame_format = "png"  # format in which each frame is saved in memory
        self.codec = "h264"  # encoding of mp4 video
        self.default_mode = "loop"  # html5 video control mode

        self.init = False

        # Needed for step functions
        self.__code__ = namedtuple("gna_hack", ["co_argcount"])
        self.__code__.co_argcount = 2

    def __call__(self, sim, todo):
        from matplotlib import pyplot as plt

        if todo == "step":
            # Initialize the plot
            if not self.init:
                filtered_plot2D = filter_dict(self.customization_args, plot2D)
                ax = sim.plot2D(ax=self.ax, fields=self.fields, **filtered_plot2D)
                # Run the plot modifier functions
                if self.plot_modifiers:
                    for k in range(len(self.plot_modifiers)):
                        ax = self.plot_modifiers[k](self.ax)
                # Store the fields
                if mp.am_master():
                    fields = ax.images[-1].get_array()
                    self.ax = ax
                    self.w, self.h = self.f.get_size_inches()
                self.init = True
            else:
                # Update the plot
                filtered_plot_fields = filter_dict(self.customization_args, plot_fields)
                fields = sim.plot_fields(fields=self.fields, **filtered_plot_fields)
                if mp.am_master():
                    self.ax.images[-1].set_data(fields)
                    self.ax.images[-1].set_clim(
                        vmin=0.8 * np.min(fields), vmax=0.8 * np.max(fields)
                    )

            if self.realtime and mp.am_master():
                # Redraw the current figure if requested
                plt.pause(0.05)

            if self.normalize and mp.am_master():
                # Save fields as a numpy array to be normalized
                # and saved later.
                self.cumulative_fields.append(fields)
            elif mp.am_master():
                # Capture figure as a png, but store the png in memory
                # to avoid writing to disk.
                self.grab_frame()
            return
        elif todo == "finish":

            # Normalize the frames, if requested, and export
            if self.normalize and mp.am_master():
                if mp.verbosity.meep > 0:
                    print("Normalizing field data...")
                fields = np.array(self.cumulative_fields) / np.max(
                    np.abs(self.cumulative_fields), axis=(0, 1, 2)
                )
                for k in range(len(self.cumulative_fields)):
                    self.ax.images[-1].set_data(fields[k, :, :])
                    self.ax.images[-1].set_clim(vmin=-0.8, vmax=0.8)
                    self.grab_frame()

            return

    @property
    def frame_size(self):
        # A tuple ``(width, height)`` in pixels of a movie frame.
        # modified from matplotlib library
        w, h = self.f.get_size_inches()
        return int(w * self.f.dpi), int(h * self.f.dpi)

    def grab_frame(self):
        # Saves the figures frame to memory.
        # modified from matplotlib library
        from io import BytesIO

        bin_data = BytesIO()
        self.f.savefig(bin_data, format=self.frame_format)
        # imgdata64 = base64.encodebytes(bin_data.getvalue()).decode('ascii')
        self._saved_frames.append(bin_data.getvalue())

    def _embedded_frames(self, frame_list, frame_format):
        # converts frame data stored in memory to html5 friendly format
        # frame_list should be a list of base64-encoded png files
        # modified from matplotlib
        import base64

        template = '  frames[{0}] = "data:image/{1};base64,{2}"\n'
        return "\n" + "".join(
            template.format(
                i,
                frame_format,
                base64.encodebytes(frame_data).decode("ascii").replace("\n", "\\\n"),
            )
            for i, frame_data in enumerate(frame_list)
        )

    def to_jshtml(self, fps):
        """
        Outputs an interactable JSHTML animation object that is embeddable in Jupyter
        notebooks. The object is packaged with controls to manipulate the video's
        playback. User must specify a frame rate `fps` in frames per second.
        """
        # Exports a javascript enabled html object that is
        # ready for jupyter notebook embedding.
        # modified from matplotlib/animation.py code.

        # Only works with Python3 and matplotlib > 3.1.0
        from distutils.version import LooseVersion
        import matplotlib

        if LooseVersion(matplotlib.__version__) < LooseVersion("3.1.0"):
            print("-------------------------------")
            print(
                "Warning: JSHTML output is not supported with your current matplotlib build. Consider upgrading to 3.1.0+"
            )
            print("-------------------------------")
            return
        if mp.am_master():
            from uuid import uuid4
            from matplotlib._animation_data import (
                DISPLAY_TEMPLATE,
                INCLUDED_FRAMES,
                JS_INCLUDE,
                STYLE_INCLUDE,
            )

            # save the frames to an html file
            fill_frames = self._embedded_frames(self._saved_frames, self.frame_format)
            Nframes = len(self._saved_frames)
            mode_dict = dict(once_checked="", loop_checked="", reflect_checked="")
            mode_dict[self.default_mode + "_checked"] = "checked"

            interval = 1000 // fps

            html_string = ""
            html_string += JS_INCLUDE
            html_string += STYLE_INCLUDE
            html_string += DISPLAY_TEMPLATE.format(
                id=uuid4().hex,
                Nframes=Nframes,
                fill_frames=fill_frames,
                interval=interval,
                **mode_dict
            )
            return JS_Animation(html_string)

    def to_gif(self, fps, filename):
        """
        Generates and outputs a GIF file of the animation with the filename, `filename`,
        and the frame rate, `fps`. Note that GIFs are significantly larger than mp4 videos
        since they don't use any compression. Artifacts are also common because the GIF
        format only supports 256 colors from a _predefined_ color palette. Requires
        `ffmpeg`.
        """
        # Exports a gif of the recorded animation
        # requires ffmpeg to be installed
        # modified from the matplotlib library
        if mp.am_master():
            from subprocess import Popen, PIPE
            from io import TextIOWrapper, BytesIO

            FFMPEG_BIN = "ffmpeg"
            command = [
                FFMPEG_BIN,
                "-f",
                "image2pipe",  # force piping of rawvideo
                "-vcodec",
                self.frame_format,  # raw input codec
                "-s",
                "%dx%d" % (self.frame_size),
                "-r",
                str(fps),  # frame rate in frames per second
                "-i",
                "pipe:",  # The input comes from a pipe
                "-vcodec",
                "gif",  # output gif format
                "-r",
                str(fps),  # frame rate in frames per second
                "-y",
                "-vf",
                "pad=width=ceil(iw/2)*2:height=ceil(ih/2)*2",
                "-an",
                filename,  # output filename
            ]
            if mp.verbosity.meep > 0:
                print("Generating GIF...")
            proc = Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            for i in range(len(self._saved_frames)):
                proc.stdin.write(self._saved_frames[i])
            out, err = proc.communicate()  # pipe in images
            proc.stdin.close()
            proc.wait()
        return

    def to_mp4(self, fps, filename):
        """
        Generates and outputs an mp4 video file of the animation with the filename,
        `filename`, and the frame rate, `fps`. Default encoding is h264 with yuv420p
        format. Requires `ffmpeg`.
        """
        # Exports an mp4 of the recorded animation
        # requires ffmpeg to be installed
        # modified from the matplotlib library
        if mp.am_master():
            from subprocess import Popen, PIPE
            from io import TextIOWrapper, BytesIO

            FFMPEG_BIN = "ffmpeg"
            command = [
                FFMPEG_BIN,
                "-f",
                "image2pipe",  # force piping of rawvideo
                "-vcodec",
                self.frame_format,  # raw input codec
                "-s",
                "%dx%d" % (self.frame_size),
                #'-pix_fmt', self.frame_format,
                "-r",
                str(fps),  # frame rate in frames per second
                "-i",
                "pipe:",  # The input comes from a pipe
                "-vcodec",
                self.codec,  # output mp4 format
                "-pix_fmt",
                "yuv420p",
                "-r",
                str(fps),  # frame rate in frames per second
                "-y",
                "-vf",
                "pad=width=ceil(iw/2)*2:height=ceil(ih/2)*2",
                "-an",
                filename,  # output filename
            ]
            if mp.verbosity.meep > 0:
                print("Generating MP4...")
            proc = Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            for i in range(len(self._saved_frames)):
                proc.stdin.write(self._saved_frames[i])
            out, err = proc.communicate()  # pipe in images
            proc.stdin.close()
            proc.wait()
        return

    def reset(self):
        self.cumulative_fields = []
        self.ax = None
        self.f = None

    def set_figure(self, f):
        self.f = f


# ------------------------------------------------------- #
# JS_Animation
# ------------------------------------------------------- #
# A helper class used to make jshtml animations embed
# seamlessly within Jupyter notebooks.


class JS_Animation:
    def __init__(self, jshtml):
        self.jshtml = jshtml

    def _repr_html_(self):
        return self.jshtml

    def get_jshtml(self):
        return self.jshtml
