import numpy as np
import matplotlib

# ------------------------------------------------------- #
# Animation
# ------------------------------------------------------- #
# A "mixin" class that extends the Simulation class. It 
# contains of the necesarry visualization routines.

class Animation(object):
    def intersect_plane_line(self,plane_0,plane_n,line_0,line_1):
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
            return [line_0,line_1]

        pt = Vector3()
        pt.x = x1 - a*(A*x1 + B*y1 + C*z1 + D) / den
        pt.y = y1 - b*(A*x1 + B*y1 + C*z1 + D) / den
        pt.z = z1 - c*(A*x1 + B*y1 + C*z1 + D) / den
        return pt

    def intersect_volume_plane(self,volume,plane):
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
            pt = self.intersect_plane_line(plane_0,plane_n,ce[0],ce[1])
            if isinstance(pt,(list,)):
                for pt_iter in pt:
                    if (pt_iter is not None) and volume.pt_in_volume(pt_iter):
                        intersection_vertices.append(pt_iter)  
            else:
                if (pt is not None) and volume.pt_in_volume(pt):
                    intersection_vertices.append(pt)
        # For point sources, check if point lies on plane
        if (volume.size.x == volume.center.y == volume.center.z == 0) and plane.pt_in_volume(volume.size):
            intersection_vertices.append(volume.size)
        return intersection_vertices
    
    def plot_volume(self,ax,volume,x,y,z,color='r'):
        # Import libraries
        try:
            import matplotlib.patches as patches
        except ImportError:
            warnings.warn("matplotlib is required for visualize_domain", ImportWarning)
            return

        # Pull parameters
        size = volume.size
        center = volume.center

        xmax = center.x+size.x/2
        xmin = center.x-size.x/2
        ymax = center.y+size.y/2
        ymin = center.y-size.y/2
        zmax = center.z+size.z/2
        zmin = center.z-size.z/2

        # Get intersecting plane
        if x is not None:
            plane = Volume(center=Vector3(x,self.geometry_center.y,self.geometry_center.z),size=Vector3(0,self.cell_size.y,self.cell_size.z))
        elif y is not None:
            plane = Volume(center=Vector3(self.geometry_center.x,y,self.geometry_center.z),size=Vector3(self.cell_size.x,0,self.cell_size.z))
        elif z is not None:
            plane = Volume(center=Vector3(self.geometry_center.x,self.geometry_center.y,z),size=Vector3(self.cell_size.x,self.cell_size.y,0))

        # Intersect plane with volume
        intersection = self.intersect_volume_plane(volume,plane)

        # Sort the points in a counter clockwise manner to ensure convex polygon is formed
        def sort_points(xy):
            xy = np.squeeze(xy)
            theta = np.arctan2(xy[:,1],xy[:,0])
            return xy[np.argsort(theta,axis=0)]

        # Point volume
        if len(intersection) == 1:
            if x == center.x:
                ax.scatter(center.x,center.z, color=color)
                return ax
            elif y == center.y:
                ax.scatter(center.y,center.z, color=color)
                return ax
            elif z == center.z:
                ax.scatter(center.x,center.y, color=color)
                return ax
            else:
                return ax
        
        # Line volume
        elif len(intersection) == 2:
            # Plot YZ
            if x is not None:
                ax.plot([a.y for a in intersection],[a.z for a in intersection], color=color)
                return ax
            #Plot XZ
            elif y is not None:
                ax.plot([a.x for a in intersection],[a.z for a in intersection], color=color)
                return ax
            # Plot XY
            elif z is not None:
                ax.plot([a.x for a in intersection],[a.y for a in intersection], color=color)
                return ax
            else:
                return ax
        
        # Planar volume
        elif len(intersection) > 2:
            # Plot YZ
            if x is not None:
                ax.add_patch(patches.Polygon(sort_points([[a.y,a.z] for a in intersection]), edgecolor=color, facecolor='none', hatch='/'))
                return ax
            #Plot XZ
            elif y is not None:
                ax.add_patch(patches.Polygon(sort_points([[a.x,a.z] for a in intersection]), edgecolor=color, facecolor='none', hatch='/'))
                return ax
            # Plot XY
            elif z is not None:
                ax.add_patch(patches.Polygon(sort_points([[a.x,a.y] for a in intersection]), edgecolor=color, facecolor='none', hatch='/'))
                return ax
            else:
                return ax
        else:
            return ax
    
    def plot_eps(self,ax,x,y,z,labels):
        # Get domain measurements
        xmin = self.geometry_center.x - self.cell_size.x/2
        xmax = self.geometry_center.x + self.cell_size.x/2
        ymin = self.geometry_center.y - self.cell_size.y/2
        ymax = self.geometry_center.y + self.cell_size.y/2
        zmin = self.geometry_center.z - self.cell_size.z/2
        zmax = self.geometry_center.z + self.cell_size.z/2

        if x is not None:
            # Plot y on x axis, z on y axis (YZ plane)
            center = Vector3(x,self.geometry_center.y,self.geometry_center.z)
            cell_size = Vector3(0,self.cell_size.y,self.cell_size.z)
            extent = [ymin,ymax,zmin,zmax]
            xlabel = 'Y'
            ylabel = 'Z'
        elif y is not None:
            # Plot x on x axis, z on y axis (XZ plane)
            center = Vector3(self.geometry_center.x,y,self.geometry_center.z)
            cell_size = Vector3(self.cell_size.x,0,self.cell_size.z)
            extent = [xmin,xmax,zmin,zmax]
            xlabel = 'X'
            ylabel = 'Z'
        elif z is not None:
            # Plot x on x axis, y on y axis (XY plane)
            center = Vector3(self.geometry_center.x,self.geometry_center.y,z)
            cell_size = Vector3(self.cell_size.x,self.cell_size.y)
            extent = [xmin,xmax,ymin,ymax]
            xlabel = 'X'
            ylabel = 'Y'

        eps_data = np.rot90(np.real(self.get_array(center=center, size=cell_size, component=mp.Dielectric)))
        ax.imshow(eps_data, interpolation='spline36', cmap='binary', extent=extent)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        return ax

    def plot_fields(self,ax,x,y,z,fields):
        if fields is not None:
            # user specifies a field component
            if fields in [mp.Ex, mp.Ey, mp.Ez, mp.Hx, mp.Hy, mp.Hz]:
                # Get domain measurements
                xmin = self.geometry_center.x - self.cell_size.x/2
                xmax = self.geometry_center.x + self.cell_size.x/2
                ymin = self.geometry_center.y - self.cell_size.y/2
                ymax = self.geometry_center.y + self.cell_size.y/2
                zmin = self.geometry_center.z - self.cell_size.z/2
                zmax = self.geometry_center.z + self.cell_size.z/2

                if x is not None:
                    # Plot y on x axis, z on y axis (YZ plane)
                    center = Vector3(x,self.geometry_center.y,self.geometry_center.z)
                    cell_size = Vector3(0,self.cell_size.y,self.cell_size.z)
                    extent = [ymin,ymax,zmin,zmax]
                    xlabel = 'Y'
                    ylabel = 'Z'
                elif y is not None:
                    # Plot x on x axis, z on y axis (XZ plane)
                    center = Vector3(self.geometry_center.x,y,self.geometry_center.z)
                    cell_size = Vector3(self.cell_size.x,0,self.cell_size.z)
                    extent = [xmin,xmax,zmin,zmax]
                    xlabel = 'X'
                    ylabel = 'Z'
                elif z is not None:
                    # Plot x on x axis, y on y axis (XY plane)
                    center = Vector3(self.geometry_center.x,self.geometry_center.y,z)
                    cell_size = Vector3(self.cell_size.x,self.cell_size.y)
                    extent = [xmin,xmax,ymin,ymax]
                    xlabel = 'X'
                    ylabel = 'Y'
                fields = self.get_array(center=center, size=cell_size, component=fields)
            # Plot whether user specifies a field component, or supplies an array
            ax.imshow(np.rot90(fields), interpolation='spline36', cmap='RdBu', alpha=0.6, extent=extent)
            return ax
        else:
            return ax

    def plot_boundaries(self,ax,x,y,z):
        # Import libraries
        try:
            import matplotlib.patches as patches
        except ImportError:
            warnings.warn("matplotlib is required for visualize_domain", ImportWarning)
            return
        
        def get_boundary_volumes(thickness,direction,side):
            thickness = boundary.thickness

            xmin = self.geometry_center.x - self.cell_size.x/2
            xmax = self.geometry_center.x + self.cell_size.x/2
            ymin = self.geometry_center.y - self.cell_size.y/2
            ymax = self.geometry_center.y + self.cell_size.y/2
            zmin = self.geometry_center.z - self.cell_size.z/2
            zmax = self.geometry_center.z + self.cell_size.z/2

            cell_x = self.cell_size.x
            cell_y = self.cell_size.y
            cell_z = self.cell_size.z

            if direction == mp.X and side == mp.Low:
                return Volume(center=Vector3(xmin+thickness/2,self.geometry_center.y,self.geometry_center.z),
                size=Vector3(thickness,cell_y,cell_z))
            elif direction == mp.X and side == mp.High:
                return Volume(center=Vector3(xmax-thickness/2,self.geometry_center.y,self.geometry_center.z),
                size=Vector3(thickness,cell_y,cell_z))
            elif direction == mp.Y and side == mp.Low:
                return Volume(center=Vector3(self.geometry_center.x,ymin+thickness/2,self.geometry_center.z),
                size=Vector3(cell_x,thickness,cell_z))
            elif direction == mp.Y and side == mp.High:
                return Volume(center=Vector3(self.geometry_center.x,ymax-thickness/2,self.geometry_center.z),
                size=Vector3(cell_x,thickness,cell_z))
            elif direction == mp.Z and side == mp.Low:
                return Volume(center=Vector3(self.geometry_center.x,self.geometry_center.y,zmin+thickness/2),
                size=Vector3(cell_x,cell_y,thickness))
            elif direction == mp.Z and side == mp.High:
                return Volume(center=Vector3(self.geometry_center.x,self.geometry_center.y,zmax-thickness/2),
                size=Vector3(cell_x,cell_y,thickness))
            else:
                raise ValueError("Invalid boundary type")
        
        import itertools
        for boundary in self.boundary_layers:
            # All 4 side are the same
            if boundary.direction == mp.ALL and boundary.side == mp.ALL:
                if self.dimensions == 1:
                    dims = [mp.X]
                elif self.dimensions == 2:
                    dims = [mp.X,mp.Y]
                elif self.dimensions == 3:
                    dims = [mp.X,mp.Y,mp.Z]
                else:
                    raise ValueError("Invalid simulation dimensions")
                for permutation in itertools.product(dims, [mp.Low, mp.High]):
                    vol = get_boundary_volumes(boundary.thickness,*permutation)
                    ax = self.plot_volume(ax,vol,x,y,z,color='g')
            # 2 sides are the same
            elif boundary.side == mp.ALL:
                for side in [mp.Low, mp.High]:
                    vol = get_boundary_volumes(boundary.thickness,*permutation)
                    ax = self.plot_volume(ax,vol,x,y,z,color='g')
            # only one side
            else:
                vol = get_boundary_volumes(boundary.thickness,boundary.direction,boundary.side)
                ax = self.plot_volume(ax,vol,x,y,z,color='g')
        return ax

    def plot_sources(self,ax,x,y,z):
        for src in self.sources:
            vol = Volume(center=src.center,size=src.size)
            ax = self.plot_volume(ax,vol,x,y,z,'r')
        return ax
    
    def plot_monitors(self,ax,x,y,z):        
        for mon in self.dft_objects:
            for reg in mon.regions:
                vol = Volume(center=reg.center,size=reg.size)
                ax = self.plot_volume(ax,vol,x,y,z,'b')
        return ax

    def plot2D(self,ax=None,x=None,y=None,z=0,fields=None,labels=True):
        if not self._is_initialized:
            self.init_sim()

        if ax is None:
            ax = plt.gca()
        
        # Plot geometry
        ax = self.plot_eps(ax,x,y,z,labels)
        
        # Plot boundaries
        ax = self.plot_boundaries(ax,x,y,z)
                
        # Plot sources
        ax = self.plot_sources(ax,x,y,z)
            
        # Plot monitors
        ax = self.plot_monitors(ax,x,y,z)

        # Plot fields
        ax = self.plot_fields(ax,x,y,z,fields)

        return ax

# ------------------------------------------------------- #
# Animate2D
# ------------------------------------------------------- #
# An extensive run function used to visualize the fields
# of a 2D simulation after every specified time step. 
# Can also be used to create running animations.

class Animate2D(object):
    def __init__(self,sim,f,fields,x=None,y=None,z=0,interval=None,update_plot=True):
        self.f = f
        self.fields = fields
        self.x = x
        self.y = y
        self.z = z
        self.interval = interval
        self.update_plot = update_plot
        self.ax = sim.plot2D(ax=self.f.gca(),x=self.x,y=self.y,z=self.z,fields=mp.Ez)
        self.next_plot_time=sim.round_time() + interval
        self.cumulative_fields = []

        self.__code__ = namedtuple('gna_hack',['co_argcount'])
        self.__code__.co_argcount=2

        import matplotlib.pyplot as plt
        self.plt = plt
    
    def __call__(self,sim,todo):
        if todo == 'step':
            if sim.round_time()<self.next_plot_time:
                return
            self.next_plot_time = sim.round_time() + self.interval
            self.ax.images[-1].remove()
            self.ax = sim.plot_fields(ax=self.ax,x=self.x,y=self.y,z=self.z,fields=self.fields)
            self.cumulative_fields.append(self.ax.images[-1].get_array())
            if self.update_plot:
                self.plt.pause(0.05)
            return
        elif todo == 'finish':
            return
    def animate_results(self,fps=20):
        from matplotlib import animation
        num_frames = len(self.cumulative_fields)

        # Normalize fields
        fields = np.array(self.cumulative_fields) / np.max(self.cumulative_fields,axis=(0,1,2))
        def animate(i):
            return self.ax.images[-1].set_data(fields[i])

        anim = animation.FuncAnimation(self.f, animate,frames=num_frames, interval=int(1/fps * 1000), blit=False)
        return anim