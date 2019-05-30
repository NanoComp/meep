from collections import namedtuple
from collections import OrderedDict
from collections import Sequence
import warnings

import numpy as np

import meep as mp
from meep.geom import Vector3, init_do_averaging
from meep.source import EigenModeSource, check_positive
    

# ------------------------------------------------------- #
# Visualization
# ------------------------------------------------------- #
# Contains all necesarry visualation routines for use with
# pymeep and pympb.

# ------------------------------------------------------- #
# Functions used to define the default plotting parameters
# for the different plotting routines.

default_source_parameters = {
        'color':'r',
        'edgecolor':'r',
        'facecolor':'none',
        'hatch':'/',
        'linewidth':2
    }

default_monitor_parameters = {
        'color':'b',
        'edgecolor':'b',
        'facecolor':'none',
        'hatch':'/',
        'linewidth':2
    }

default_field_parameters = {
        'interpolation':'spline36',
        'cmap':'RdBu',
        'alpha':0.6
        }

default_eps_parameters = {
        'interpolation':'spline36', 
        'cmap':'binary',
        'alpha':1.0
    }

default_boundary_parameters = {
        'color':'g',
        'edgecolor':'g',
        'facecolor':'none',
        'hatch':'/'
    }

default_volume_parameters = {
        'alpha':1.0,
        'color':'k',
        'linestyle':'-',
        'linewidth':1,
        'marker':'.',
        'edgecolor':'k',
        'facecolor':'none',
        'hatch':'/'
    }

default_label_parameters = {
    'label_color':'r',
    'offset':20,
    'label_alpha':0.3
}

# ------------------------------------------------------- #
# Routines to add legends to plot

def place_label(ax,label_text,x,y,centerx,centery,label_parameters=None):

    label_parameters = default_label_parameters if label_parameters is None else dict(default_label_parameters, **label_parameters)

    offset = label_parameters['offset']
    alpha = label_parameters['label_alpha']
    color = label_parameters['label_color']

    if x > centerx:
        xtext = -offset
    else:
        xtext = offset
    if y > centery:
        ytext = -offset
    else:
        ytext = offset
    
    ax.annotate(label_text, xy=(x, y), xytext=(xtext,ytext), 
            textcoords='offset points', ha='center', va='bottom',
            bbox=dict(boxstyle='round,pad=0.2', fc=color, alpha=alpha),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5', 
                            color=color))
    return ax

# ------------------------------------------------------- #
# Helper functions used to plot volumes on a 2D plane

def intersect_plane_line(plane_0,plane_n,line_0,line_1):
    # Find the intersection point of a plane and line:
    # http://www.ambrsoft.com/TrigoCalc/Plan3D/PlaneLineIntersection_.htm
    # plane_0 ........ [Vector3] origin of plane
    # plane_n ........ [Vector3] normal vector of plane
    # line_0 ......... [Vector3] first point of line
    # line_1 ......... [Vector3] second point of line
    # If the line segment is parallel with the plane, then the first
    # and last points (line_0 and line_1) are returned in a list.

    # Plane coefficients
    D = -plane_0.dot(plane_n)
    A = plane_n.x
    B = plane_n.y
    C = plane_n.z

    # Line coefficients
    v = line_1 - line_0
    a = v.x
    b = v.y
    c = v.z
    x1 = line_0.x
    y1 = line_0.y
    z1 = line_0.z

    den = A*a + B*b + C*c

    # parallel case
    if den == 0:
        # coplanar
        if (line_0-plane_0).dot(plane_n) == 0:
            return [line_0,line_1]
        # just parallel
        else:
            return None 

    pt = Vector3()
    pt.x = x1 - a*(A*x1 + B*y1 + C*z1 + D) / den
    pt.y = y1 - b*(A*x1 + B*y1 + C*z1 + D) / den
    pt.z = z1 - c*(A*x1 + B*y1 + C*z1 + D) / den
    return pt

def intersect_volume_plane(volume,plane):
    # returns the vertices that correspond to the polygon, 
    # line, or single point of the intersection of a volume
    # object and a plane
    # volume ......... [Volume] volume object
    # plane .......... [Volume] volume object of the plane

    # Get normal vector of plane
    if plane.size.x == 0:
        plane_n = Vector3(x=1)
    elif plane.size.y == 0:
        plane_n = Vector3(y=1)
    elif plane.size.z == 0: 
        plane_n = Vector3(z=1)
    else:
        raise ValueError("plane volume must have a nonzero dimension")

    # Get origin of plane
    plane_0 = plane.center

    intersection_vertices = []
    edges = volume.get_edges()
    for ce in edges:          
        pt = intersect_plane_line(plane_0,plane_n,ce[0],ce[1])
        if isinstance(pt,(list,)):
            for pt_iter in pt:
                if (pt_iter is not None) and volume.pt_in_volume(pt_iter):
                    intersection_vertices.append(pt_iter)  
        else:
            if (pt is not None) and volume.pt_in_volume(pt):
                intersection_vertices.append(pt)
    # For point sources, check if point lies on plane
    if (volume.size.x == volume.size.y == volume.size.z == 0) and plane.pt_in_volume(volume.size):
        intersection_vertices.append(volume.size)
    return intersection_vertices

# ------------------------------------------------------- #
# actual plotting routines

def plot_volume(sim,ax,volume,output_plane=None,plotting_parameters=None,label=None):
    import matplotlib.patches as patches
    from matplotlib import pyplot as plt
    from meep.simulation import Volume

    # Set up the plotting parameters
    plotting_parameters = default_volume_parameters if plotting_parameters is None else dict(default_volume_parameters, **plotting_parameters)

    # Get domain measurements
    if output_plane:
        sim_center, sim_size = (output_plane.center, output_plane.size)
    elif sim.output_volume:
        sim_center, sim_size = mp.get_center_and_size(sim.output_volume)
    else:
        sim_center, sim_size = (sim.geometry_center, sim.cell_size)
    
    plane = Volume(center=sim_center,size=sim_size)

    # Pull volume parameters
    size = volume.size
    center = volume.center

    xmax = center.x+size.x/2
    xmin = center.x-size.x/2
    ymax = center.y+size.y/2
    ymin = center.y-size.y/2
    zmax = center.z+size.z/2
    zmin = center.z-size.z/2

    # Add labels if requested
    if label is not None:
        if sim_size.x == 0:
            ax = place_label(ax,label,center.y,center.z,sim_center.y,sim_center.z,label_parameters=plotting_parameters)
        elif sim_size.y == 0:
            ax = place_label(ax,label,center.x,center.z,sim_center.x,sim_center.z,label_parameters=plotting_parameters)
        elif sim_size.z == 0:
            ax = place_label(ax,label,center.x,center.y,sim_center.x,sim_center.y,label_parameters=plotting_parameters)

    # Intersect plane with volume
    intersection = intersect_volume_plane(volume,plane)

    # Sort the points in a counter clockwise manner to ensure convex polygon is formed
    def sort_points(xy):
        xy = np.squeeze(xy)
        theta = np.arctan2(xy[:,1],xy[:,0])
        return xy[np.argsort(theta,axis=0)]

    # Point volume
    if len(intersection) == 1:
        point_args = {key:value for key, value in plotting_parameters.items() if key in ['color','marker','alpha','linewidth']}
        if sim_center.y == center.y:
            ax.scatter(center.x,center.z, **point_args)
            return ax
        elif sim_center.x == center.x:
            ax.scatter(center.y,center.z, **point_args)
            return ax
        elif sim_center.z == center.z:
            ax.scatter(center.x,center.y, **point_args)
            return ax
        else:
            return ax
    
    # Line volume
    elif len(intersection) == 2:
        line_args = {key:value for key, value in plotting_parameters.items() if key in ['color','linestyle','linewidth','alpha']}
        # Plot YZ
        if sim_size.x == 0:
            ax.plot([a.y for a in intersection],[a.z for a in intersection], **line_args)
            return ax
        #Plot XZ
        elif sim_size.y == 0:
            ax.plot([a.x for a in intersection],[a.z for a in intersection], **line_args)
            return ax
        # Plot XY
        elif sim_size.z == 0:
            ax.plot([a.x for a in intersection],[a.y for a in intersection], **line_args)
            return ax
        else:
            return ax
    
    # Planar volume
    elif len(intersection) > 2:
        planar_args = {key:value for key, value in plotting_parameters.items() if key in ['edgecolor','linewidth','facecolor','hatch','alpha']}
        # Plot YZ
        if sim_size.x == 0:
            ax.add_patch(patches.Polygon(sort_points([[a.y,a.z] for a in intersection]), **planar_args))
            return ax
        #Plot XZ
        elif sim_size.y == 0:
            ax.add_patch(patches.Polygon(sort_points([[a.x,a.z] for a in intersection]), **planar_args))
            return ax
        # Plot XY
        elif sim_size.z == 0:
            ax.add_patch(patches.Polygon(sort_points([[a.x,a.y] for a in intersection]), **planar_args))
            return ax
        else:
            return ax
    else:
        return ax

def plot_eps(sim,ax,output_plane=None,eps_parameters=None):
    
    # consolidate plotting parameters
    eps_parameters = default_eps_parameters if eps_parameters is None else dict(default_eps_parameters, **eps_parameters)
    
    # Get domain measurements
    if output_plane:
        sim_center, sim_size = (output_plane.center, output_plane.size)
    elif sim.output_volume:
        sim_center, sim_size = mp.get_center_and_size(sim.output_volume)
    else:
        sim_center, sim_size = (sim.geometry_center, sim.cell_size)

    xmin = sim_center.x - sim_size.x/2
    xmax = sim_center.x + sim_size.x/2
    ymin = sim_center.y - sim_size.y/2
    ymax = sim_center.y + sim_size.y/2
    zmin = sim_center.z - sim_size.z/2
    zmax = sim_center.z + sim_size.z/2

    center = Vector3(sim_center.x,sim_center.y,sim_center.z)
    cell_size = Vector3(sim_size.x,sim_size.y,sim_size.z)

    if sim_size.x == 0:
        # Plot y on x axis, z on y axis (YZ plane)
        extent = [ymin,ymax,zmin,zmax]
        xlabel = 'Y'
        ylabel = 'Z'
    elif sim_size.y == 0:
        # Plot x on x axis, z on y axis (XZ plane)
        extent = [xmin,xmax,zmin,zmax]
        xlabel = 'X'
        ylabel = 'Z'
    elif sim_size.z == 0:
        # Plot x on x axis, y on y axis (XY plane)
        extent = [xmin,xmax,ymin,ymax]
        xlabel = 'X'
        ylabel = 'Y'
    else:
        raise ValueError("A 2D plane has not been specified...")

    eps_data = np.rot90(np.real(sim.get_array(center=center, size=cell_size, component=mp.Dielectric)))
    ax.imshow(eps_data, extent=extent, **eps_parameters)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return ax

def plot_boundaries(sim,ax,output_plane=None,boundary_parameters=None):

    # consolidate plotting parameters
    boundary_parameters = default_boundary_parameters if boundary_parameters is None else dict(default_boundary_parameters, **boundary_parameters)

    def get_boundary_volumes(thickness,direction,side):
        from meep.simulation import Volume

        thickness = boundary.thickness

        # Get domain measurements
        sim_center, sim_size = (sim.geometry_center, sim.cell_size)

        xmin = sim_center.x - sim_size.x/2
        xmax = sim_center.x + sim_size.x/2
        ymin = sim_center.y - sim_size.y/2
        ymax = sim_center.y + sim_size.y/2
        zmin = sim_center.z - sim_size.z/2
        zmax = sim_center.z + sim_size.z/2

        cell_x = sim.cell_size.x
        cell_y = sim.cell_size.y
        cell_z = sim.cell_size.z

        if direction == mp.X and side == mp.Low:
            return Volume(center=Vector3(xmin+thickness/2,sim.geometry_center.y,sim.geometry_center.z),
            size=Vector3(thickness,cell_y,cell_z))
        elif direction == mp.X and side == mp.High:
            return Volume(center=Vector3(xmax-thickness/2,sim.geometry_center.y,sim.geometry_center.z),
            size=Vector3(thickness,cell_y,cell_z))
        elif direction == mp.Y and side == mp.Low:
            return Volume(center=Vector3(sim.geometry_center.x,ymin+thickness/2,sim.geometry_center.z),
            size=Vector3(cell_x,thickness,cell_z))
        elif direction == mp.Y and side == mp.High:
            return Volume(center=Vector3(sim.geometry_center.x,ymax-thickness/2,sim.geometry_center.z),
            size=Vector3(cell_x,thickness,cell_z))
        elif direction == mp.Z and side == mp.Low:
            return Volume(center=Vector3(sim.geometry_center.x,sim.geometry_center.y,zmin+thickness/2),
            size=Vector3(cell_x,cell_y,thickness))
        elif direction == mp.Z and side == mp.High:
            return Volume(center=Vector3(sim.geometry_center.x,sim.geometry_center.y,zmax-thickness/2),
            size=Vector3(cell_x,cell_y,thickness))
        else:
            raise ValueError("Invalid boundary type")
    
    import itertools
    for boundary in sim.boundary_layers:
        # All 4 side are the same
        if boundary.direction == mp.ALL and boundary.side == mp.ALL:
            if sim.dimensions == 1:
                dims = [mp.X]
            elif sim.dimensions == 2:
                dims = [mp.X,mp.Y]
            elif sim.dimensions == 3:
                dims = [mp.X,mp.Y,mp.Z]
            else:
                raise ValueError("Invalid simulation dimensions")
            for permutation in itertools.product(dims, [mp.Low, mp.High]):
                vol = get_boundary_volumes(boundary.thickness,*permutation)
                ax = plot_volume(sim,ax,vol,output_plane,plotting_parameters=boundary_parameters)
        # 2 sides are the same
        elif boundary.side == mp.ALL:
            for side in [mp.Low, mp.High]:
                vol = get_boundary_volumes(boundary.thickness,boundary.direction,side)
                ax = plot_volume(sim,ax,vol,output_plane,plotting_parameters=boundary_parameters)
        # only one side
        else:
            vol = get_boundary_volumes(boundary.thickness,boundary.direction,boundary.side)
            ax = plot_volume(sim,ax,vol,output_plane,plotting_parameters=boundary_parameters)
    return ax

def plot_sources(sim,ax,output_plane=None,labels=False,source_parameters=None):
    from meep.simulation import Volume

    # consolidate plotting parameters
    source_parameters = default_source_parameters if source_parameters is None else dict(default_source_parameters, **source_parameters)

    label = 'source' if labels else None

    for src in sim.sources:
        vol = Volume(center=src.center,size=src.size)
        ax = plot_volume(sim,ax,vol,output_plane,plotting_parameters=source_parameters,label=label)
    return ax

def plot_monitors(sim,ax,output_plane=None,labels=False,monitor_parameters=None):
    from meep.simulation import Volume

     # consolidate plotting parameters
    monitor_parameters = default_monitor_parameters if monitor_parameters is None else dict(default_monitor_parameters, **monitor_parameters)

    label = 'monitor' if labels else None

    for mon in sim.dft_objects:
        for reg in mon.regions:
            vol = Volume(center=reg.center,size=reg.size)
            ax = plot_volume(sim,ax,vol,output_plane,plotting_parameters=monitor_parameters,label=label)
    return ax

def plot_fields(sim,ax=None,fields=None,output_plane=None,field_parameters=None):

    if fields is None:
        return ax
    
    # user specifies a field component
    if fields in [mp.Ex, mp.Ey, mp.Ez, mp.Hx, mp.Hy, mp.Hz]:
        # Get domain measurements
        if output_plane:
            sim_center, sim_size = (output_plane.center, output_plane.size)
        elif sim.output_volume:
            sim_center, sim_size = mp.get_center_and_size(sim.output_volume)
        else:
            sim_center, sim_size = (sim.geometry_center, sim.cell_size)
        
        xmin = sim_center.x - sim_size.x/2
        xmax = sim_center.x + sim_size.x/2
        ymin = sim_center.y - sim_size.y/2
        ymax = sim_center.y + sim_size.y/2
        zmin = sim_center.z - sim_size.z/2
        zmax = sim_center.z + sim_size.z/2

        center = Vector3(sim_center.x,sim_center.y,sim_center.z)
        cell_size = Vector3(sim_size.x,sim_size.y,sim_size.z)

        if sim_size.x == 0:
            # Plot y on x axis, z on y axis (YZ plane)
            extent = [ymin,ymax,zmin,zmax]
            xlabel = 'Y'
            ylabel = 'Z'
        elif sim_size.y == 0:
            # Plot x on x axis, z on y axis (XZ plane)
            extent = [xmin,xmax,zmin,zmax]
            xlabel = 'X'
            ylabel = 'Z'
        elif sim_size.z == 0:
            # Plot x on x axis, y on y axis (XY plane)
            extent = [xmin,xmax,ymin,ymax]
            xlabel = 'X'
            ylabel = 'Y'
        fields = sim.get_array(center=center, size=cell_size, component=fields)
    else:
        raise ValueError('Please specify a valid field component (mp.Ex, mp.Ey, ...')
    
    # Either plot the field, or return the array
    if ax:
        field_parameters = default_field_parameters if field_parameters is None else dict(default_field_parameters, **field_parameters)
        ax.imshow(np.rot90(fields), extent=extent, **field_parameters)
        return ax
    else:
        return np.rot90(fields)
    return ax

def plot2D(sim,ax=None, output_plane=None, fields=None, labels=False,
            eps_parameters=None,boundary_parameters=None,
            source_parameters=None,monitor_parameters=None,
            field_parameters=None):
    
    if not mp.am_master():
        return

    if not sim._is_initialized:
        sim.init_sim()

    if ax is None:
        from matplotlib import pyplot as plt
        ax = plt.gca()

    # User incorrectly specified a 3D output plane
    if output_plane and (output_plane.size.x != 0) and (output_plane.size.y != 0) and (output_plane.size.z != 0):
        raise ValueError("output_plane must be a 2 dimensional volume (a plane).")
    # User forgot to specify a 2D output plane for a 3D simulation
    elif output_plane is None and (sim.cell_size.x != 0) and (sim.cell_size.y != 0) and (sim.cell_size.z != 0):
        raise ValueError("For 3D simulations, you must specify an output_plane.")
    
    # Plot geometry
    ax = plot_eps(sim,ax,output_plane=output_plane,eps_parameters=eps_parameters)
    
    # Plot boundaries
    ax = plot_boundaries(sim,ax,output_plane=output_plane,boundary_parameters=boundary_parameters)
            
    # Plot sources
    ax = plot_sources(sim,ax,output_plane=output_plane,labels=labels,source_parameters=source_parameters)
        
    # Plot monitors
    ax = plot_monitors(sim,ax,output_plane=output_plane,labels=labels,monitor_parameters=monitor_parameters)

    # Plot fields
    ax = plot_fields(sim,ax,fields,output_plane=output_plane,field_parameters=field_parameters)

    return ax

def plot3D(sim):
    from mayavi import mlab

    if not sim._is_initialized:
        sim.init_sim()
        
    if sim.dimensions < 3:
        raise ValueError("Simulation must have 3 dimensions to visualize 3D")
    
    eps_data = sim.get_epsilon()
    s = mlab.contour3d(eps_data, colormap="YlGnBu")
    return s

def visualize_chunks(sim):
    if sim.structure is None:
        sim.init_sim()
    
    import matplotlib.pyplot as plt
    import matplotlib.cm
    import matplotlib.colors
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

    vols = sim.structure.get_chunk_volumes()
    owners = sim.structure.get_chunk_owners()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    chunk_colors = matplotlib.cm.rainbow(np.linspace(0, 1, mp.count_processors()))

    def plot_box(box, proc):
        low = mp.Vector3(box.low.x, box.low.y, box.low.z)
        high = mp.Vector3(box.high.x, box.high.y, box.high.z)
        points = [low, high]

        x_len = mp.Vector3(high.x) - mp.Vector3(low.x)
        y_len = mp.Vector3(y=high.y) - mp.Vector3(y=low.y)
        xy_len = mp.Vector3(high.x, high.y) - mp.Vector3(low.x, low.y)

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
            [points[0], points[7], points[6], points[2]]
        ]

        faces = Poly3DCollection(edges, linewidths=1, edgecolors='k')
        color_with_alpha = matplotlib.colors.to_rgba(chunk_colors[proc], alpha=0.2)
        faces.set_facecolor(color_with_alpha)
        ax.add_collection3d(faces)

        # Plot the points themselves to force the scaling of the axes
        ax.scatter(points[:, 0], points[:, 1], points[:, 2], s=0)

    if mp.am_master():
        for i, v in enumerate(vols):
            plot_box(mp.gv2box(v.surroundings()), owners[i])
        ax.set_aspect('auto')
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

class Animate2D(object):
    def __init__(self,sim,fields,f=None,realtime=True,normalize=False,
    plot_modifiers=None,**customization_args):
        self.fields = fields

        from matplotlib import pyplot as plt
        from matplotlib import animation

        if f:
            self.f = f
        else:
            self.f = plt.figure()

        self.realtime = realtime
        self.normalize = normalize
        self.plot_modifiers = plot_modifiers
        self.customization_args = customization_args

        self.ax = None

        self.cumulative_fields = []
        self._saved_frames = []

        self.frame_format = 'png' # format in which each frame is saved in memory
        self.codec = 'h264' # encoding of mp4 video
        self.default_mode = 'loop' # html5 video control mode
        
        # Needed for step functions
        self.__code__ = namedtuple('gna_hack',['co_argcount'])
        self.__code__.co_argcount=2
    
    def __call__(self,sim,todo):
        from matplotlib import pyplot as plt
        
        if todo == 'step':

            # Initialize the plot
            if self.ax is None:
                self.ax = sim.plot2D(ax=self.f.gca(),fields=self.fields,**self.customization_args)
                # Run the plot modifier functions
                if self.plot_modifiers:
                    for k in range(len(self.plot_modifiers)):
                        self.ax = self.plot_modifiers[k](self.ax)
                # Store the fields
                self.w, self.h = self.f.get_size_inches()
                fields = self.ax.images[-1].get_array()
            else:                
                # Update the plot
                fields = sim.plot_fields(fields=self.fields)
                self.ax.images[-1].set_data(fields)
                self.ax.images[-1].set_clim(vmin=0.8*np.min(fields), vmax=0.8*np.max(fields))
            
            if self.realtime:
                # Redraw the current figure if requested
                plt.pause(0.05)

            if self.normalize:
                # Save fields as a numpy array to be normalized 
                # and saved later.
                self.cumulative_fields.append(fields)
            else:
                # Capture figure as a png, but store the png in memory
                # to avoid writing to disk.
                self.grab_frame()
            return
        elif todo == 'finish':
            
            # Normalize the frames, if requested, and export
            if self.normalize:
                print("Normalizing field data...")
                fields = np.array(self.cumulative_fields) / np.max(np.abs(self.cumulative_fields),axis=(0,1,2))
                for k in range(len(self.cumulative_fields)):
                    self.ax.images[-1].set_data(fields[k,:,:])
                    self.ax.images[-1].set_clim(vmin=-0.8, vmax=0.8)
                    self.grab_frame()
            
            plt.close(self.f)
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
        #imgdata64 = base64.encodebytes(bin_data.getvalue()).decode('ascii')
        self._saved_frames.append(bin_data.getvalue())

    def _embedded_frames(self, frame_list, frame_format):
        # converts frame data stored in memory to html5 friendly format
        # frame_list should be a list of base64-encoded png files
        # modified from matplotlib
        import base64
        template = '  frames[{0}] = "data:image/{1};base64,{2}"\n'
        return "\n" + "".join(
            template.format(i, frame_format, base64.encodebytes(frame_data).decode('ascii').replace('\n', '\\\n'))
            for i, frame_data in enumerate(frame_list))

    def to_jshtml(self,fps):
        # Exports a javascript enabled html object that is
        # ready for jupyter notebook embedding.
        # modified from matplotlib/animation.py code.
        from uuid import uuid4

        # save the frames to an html file
        fill_frames = self._embedded_frames(self._saved_frames, self.frame_format)
        Nframes = len(self._saved_frames)
        mode_dict = dict(once_checked='',
                         loop_checked='',
                         reflect_checked='')
        mode_dict[self.default_mode + '_checked'] = 'checked'

        interval = 1000 // fps
        
        html_string = ""
        html_string += JS_INCLUDE
        html_string += STYLE_INCLUDE
        html_string += DISPLAY_TEMPLATE.format(id=uuid4().hex,
                                             Nframes=Nframes,
                                             fill_frames=fill_frames,
                                             interval=interval,
                                             **mode_dict)
        return JS_Animation(html_string)

    def to_gif(self,fps,filename):
        # Exports a gif of the recorded animation
        # requires ffmpeg to be installed
        # modified from the matplotlib library
        from subprocess import Popen, PIPE
        from io import TextIOWrapper, BytesIO
        FFMPEG_BIN = 'ffmpeg'
        command = [FFMPEG_BIN, 
                    '-f', 'image2pipe', # force piping of rawvideo
                    '-vcodec', self.frame_format, # raw input codec
                    '-s', '%dx%d' % (self.frame_size),
                    '-r', str(fps), # frame rate in frames per second
                    '-i', 'pipe:', # The input comes from a pipe
                    '-vcodec', 'gif', # output gif format
                    '-r', str(fps), # frame rate in frames per second
                    '-y', 
                    '-vf', 'pad=width=ceil(iw/2)*2:height=ceil(ih/2)*2',
                    '-an', filename  # output filename
        ]

        print("Generating GIF...")
        proc = Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        for i in range(len(self._saved_frames)): 
            proc.stdin.write(self._saved_frames[i])
        out, err = proc.communicate() # pipe in images
        proc.stdin.close()
        proc.wait()
        return

    def to_mp4(self,fps,filename):
        # Exports an mp4 of the recorded animation
        # requires ffmpeg to be installed
        # modified from the matplotlib library
        from subprocess import Popen, PIPE
        from io import TextIOWrapper, BytesIO
        FFMPEG_BIN = 'ffmpeg'
        command = [FFMPEG_BIN, 
                    '-f', 'image2pipe', # force piping of rawvideo
                    '-vcodec', self.frame_format, # raw input codec
                    '-s', '%dx%d' % (self.frame_size),
                    #'-pix_fmt', self.frame_format,
                    '-r', str(fps), # frame rate in frames per second
                    '-i', 'pipe:', # The input comes from a pipe
                    '-vcodec', self.codec, # output mp4 format
                    '-pix_fmt','yuv420p',
                    '-r', str(fps), # frame rate in frames per second
                    '-y', 
                    '-vf', 'pad=width=ceil(iw/2)*2:height=ceil(ih/2)*2',
                    '-an', filename  # output filename
        ]
        print("Generating MP4...")
        proc = Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        for i in range(len(self._saved_frames)): 
            proc.stdin.write(self._saved_frames[i])
        out, err = proc.communicate() # pipe in images
        proc.stdin.close()
        proc.wait()
        return
    
    def reset(self):
        self.cumulative_fields = []
        self.ax = None
        self.f = None

    def set_figure(self,f):
        self.f = f

# ------------------------------------------------------- #
# JS_Animation
# ------------------------------------------------------- #
# A helper class used to make jshtml animations embed 
# seamlessly within Jupyter notebooks.

class JS_Animation():
    def __init__(self,jshtml):
        self.jshtml = jshtml
    def _repr_html_(self):
        return self.jshtml
    def get_jshtml(self):
        return self.jshtml


# ------------------------------------------------------- #
# JSHTML Templates
# ------------------------------------------------------- #
# Template functions grabbed from matplotlib/_animation_data.py
# https://github.com/matplotlib/matplotlib/blob/aac57266803ff09c8c9dc10e1bb2bd389d79e42d/lib/matplotlib/_animation_data.py
# Included here since python2 only supports matplotlib <3.0


# Javascript template for HTMLWriter
JS_INCLUDE = """
<link rel="stylesheet"
href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/
css/font-awesome.min.css">
<script language="javascript">
  function isInternetExplorer() {
    ua = navigator.userAgent;
    /* MSIE used to detect old browsers and Trident used to newer ones*/
    return ua.indexOf("MSIE ") > -1 || ua.indexOf("Trident/") > -1;
  }
  /* Define the Animation class */
  function Animation(frames, img_id, slider_id, interval, loop_select_id){
    this.img_id = img_id;
    this.slider_id = slider_id;
    this.loop_select_id = loop_select_id;
    this.interval = interval;
    this.current_frame = 0;
    this.direction = 0;
    this.timer = null;
    this.frames = new Array(frames.length);
    for (var i=0; i<frames.length; i++)
    {
     this.frames[i] = new Image();
     this.frames[i].src = frames[i];
    }
    var slider = document.getElementById(this.slider_id);
    slider.max = this.frames.length - 1;
    if (isInternetExplorer()) {
        // switch from oninput to onchange because IE <= 11 does not conform
        // with W3C specification. It ignores oninput and onchange behaves
        // like oninput. In contrast, Mircosoft Edge behaves correctly.
        slider.setAttribute('onchange', slider.getAttribute('oninput'));
        slider.setAttribute('oninput', null);
    }
    this.set_frame(this.current_frame);
  }
  Animation.prototype.get_loop_state = function(){
    var button_group = document[this.loop_select_id].state;
    for (var i = 0; i < button_group.length; i++) {
        var button = button_group[i];
        if (button.checked) {
            return button.value;
        }
    }
    return undefined;
  }
  Animation.prototype.set_frame = function(frame){
    this.current_frame = frame;
    document.getElementById(this.img_id).src =
            this.frames[this.current_frame].src;
    document.getElementById(this.slider_id).value = this.current_frame;
  }
  Animation.prototype.next_frame = function()
  {
    this.set_frame(Math.min(this.frames.length - 1, this.current_frame + 1));
  }
  Animation.prototype.previous_frame = function()
  {
    this.set_frame(Math.max(0, this.current_frame - 1));
  }
  Animation.prototype.first_frame = function()
  {
    this.set_frame(0);
  }
  Animation.prototype.last_frame = function()
  {
    this.set_frame(this.frames.length - 1);
  }
  Animation.prototype.slower = function()
  {
    this.interval /= 0.7;
    if(this.direction > 0){this.play_animation();}
    else if(this.direction < 0){this.reverse_animation();}
  }
  Animation.prototype.faster = function()
  {
    this.interval *= 0.7;
    if(this.direction > 0){this.play_animation();}
    else if(this.direction < 0){this.reverse_animation();}
  }
  Animation.prototype.anim_step_forward = function()
  {
    this.current_frame += 1;
    if(this.current_frame < this.frames.length){
      this.set_frame(this.current_frame);
    }else{
      var loop_state = this.get_loop_state();
      if(loop_state == "loop"){
        this.first_frame();
      }else if(loop_state == "reflect"){
        this.last_frame();
        this.reverse_animation();
      }else{
        this.pause_animation();
        this.last_frame();
      }
    }
  }
  Animation.prototype.anim_step_reverse = function()
  {
    this.current_frame -= 1;
    if(this.current_frame >= 0){
      this.set_frame(this.current_frame);
    }else{
      var loop_state = this.get_loop_state();
      if(loop_state == "loop"){
        this.last_frame();
      }else if(loop_state == "reflect"){
        this.first_frame();
        this.play_animation();
      }else{
        this.pause_animation();
        this.first_frame();
      }
    }
  }
  Animation.prototype.pause_animation = function()
  {
    this.direction = 0;
    if (this.timer){
      clearInterval(this.timer);
      this.timer = null;
    }
  }
  Animation.prototype.play_animation = function()
  {
    this.pause_animation();
    this.direction = 1;
    var t = this;
    if (!this.timer) this.timer = setInterval(function() {
        t.anim_step_forward();
    }, this.interval);
  }
  Animation.prototype.reverse_animation = function()
  {
    this.pause_animation();
    this.direction = -1;
    var t = this;
    if (!this.timer) this.timer = setInterval(function() {
        t.anim_step_reverse();
    }, this.interval);
  }
</script>
"""


# Style definitions for the HTML template
STYLE_INCLUDE = """
<style>
.animation {
    display: inline-block;
    text-align: center;
}
input[type=range].anim-slider {
    width: 374px;
    margin-left: auto;
    margin-right: auto;
}
.anim-buttons {
    margin: 8px 0px;
}
.anim-buttons button {
    padding: 0;
    width: 36px;
}
.anim-state label {
    margin-right: 8px;
}
.anim-state input {
    margin: 0;
    vertical-align: middle;
}
</style>
"""


# HTML template for HTMLWriter
DISPLAY_TEMPLATE = """
<div class="animation">
  <img id="_anim_img{id}">
  <div class="anim-controls">
    <input id="_anim_slider{id}" type="range" class="anim-slider"
           name="points" min="0" max="1" step="1" value="0"
           oninput="anim{id}.set_frame(parseInt(this.value));"></input>
    <div class="anim-buttons">
      <button onclick="anim{id}.slower()"><i class="fa fa-minus"></i></button>
      <button onclick="anim{id}.first_frame()"><i class="fa fa-fast-backward">
          </i></button>
      <button onclick="anim{id}.previous_frame()">
          <i class="fa fa-step-backward"></i></button>
      <button onclick="anim{id}.reverse_animation()">
          <i class="fa fa-play fa-flip-horizontal"></i></button>
      <button onclick="anim{id}.pause_animation()"><i class="fa fa-pause">
          </i></button>
      <button onclick="anim{id}.play_animation()"><i class="fa fa-play"></i>
          </button>
      <button onclick="anim{id}.next_frame()"><i class="fa fa-step-forward">
          </i></button>
      <button onclick="anim{id}.last_frame()"><i class="fa fa-fast-forward">
          </i></button>
      <button onclick="anim{id}.faster()"><i class="fa fa-plus"></i></button>
    </div>
    <form action="#n" name="_anim_loop_select{id}" class="anim-state">
      <input type="radio" name="state" value="once" id="_anim_radio1_{id}"
             {once_checked}>
      <label for="_anim_radio1_{id}">Once</label>
      <input type="radio" name="state" value="loop" id="_anim_radio2_{id}"
             {loop_checked}>
      <label for="_anim_radio2_{id}">Loop</label>
      <input type="radio" name="state" value="reflect" id="_anim_radio3_{id}"
             {reflect_checked}>
      <label for="_anim_radio3_{id}">Reflect</label>
    </form>
  </div>
</div>
<script language="javascript">
  /* Instantiate the Animation class. */
  /* The IDs given should match those used in the template above. */
  (function() {{
    var img_id = "_anim_img{id}";
    var slider_id = "_anim_slider{id}";
    var loop_select_id = "_anim_loop_select{id}";
    var frames = new Array({Nframes});
    {fill_frames}
    /* set a timeout to make sure all the above elements are created before
       the object is initialized. */
    setTimeout(function() {{
        anim{id} = new Animation(frames, img_id, slider_id, {interval},
                                 loop_select_id);
    }}, 0);
  }})()
</script>
"""


INCLUDED_FRAMES = """
  for (var i=0; i<{Nframes}; i++){{
    frames[i] = "{frame_dir}/frame" + ("0000000" + i).slice(-7) +
                ".{frame_format}";
  }}
"""