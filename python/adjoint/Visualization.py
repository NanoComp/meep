###################################################
# visualize.py -- various routines for visualizing
# the inputs and outputs of pymeep calculations in
# simple standardized ways
###################################################
import warnings
import numpy as np
import sympy
from collections import namedtuple
import time
import datetime
from multiprocessing import Process, Pipe
from matplotlib import ticker

import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

from matplotlib.collections import PolyCollection, LineCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import meep as mp

from .ObjectiveFunction import global_dft_cell_names

########################################################
# the next few routines are some of my older general-purpose
# utilities for working with DFT cells; for that purpose
# they have been rendered obsolete by the DFTCell class
# in the adjoint module, but I reproduce them in this
# visualization module just so it can be independent of
# the adjoint innards.
########################################################
EHTransverse=[ [mp.Ey, mp.Ez, mp.Hy, mp.Hz],
               [mp.Ez, mp.Ex, mp.Hz, mp.Hx],
               [mp.Ex, mp.Ey, mp.Hx, mp.Hy] ]
Exyz=[mp.Ex, mp.Ey, mp.Ez]
Hxyz=[mp.Hx, mp.Hy, mp.Hz]
EHxyz=Exyz + Hxyz
origin=mp.Vector3()

def abs2(z):
    return np.real(np.conj(z)*z)

def dft_cell_type(cell):
        return 'flux' if isinstance(cell, mp.simulation.DftFlux) else 'fields'

def dft_cell_name(nc):
    try:
        name=global_dft_cell_names[nc]
        if name.endswith('_flux'):
            name=name[0:-5]
        if name.endswith('_fields'):
            name=name[0:-7]
        return name.replace('_','\_')
    except NameError:
        return 'flux_' + str(nc)

def is_flux_cell(cell):
    return dft_cell_type(cell)=='flux'


def add_dft_cell(sim, region, fcen, df=0, nfreq=1):
    if hasattr(region,"direction"):
        return sim.add_flux(fcen,df,nfreq,region)
    return sim.add_dft_fields(EHxyz, fcen, df,nfreq, where=region)


def flux_line(x0, y0, length, dir):
    size=length*(xhat if dir is mp.Y else yhat)
    return mp.FluxRegion(center=mp.Vector3(x0,y0), size=size, direction=dir)


def tangential_components(normal_dir):
    if normal_dir in range(mp.X, mp.Z+1):
        return EHTransverse[normal_dir-mp.X]
    raise ValueError("invalid normal_dir={} in tangential_components".format(normal_dir))


# return a list of the field components stored in a DFT cell
def dft_cell_components(cell):
    if dft_cell_type(cell)=='flux':
        return tangential_components(cell.normal_direction)
    elif cell.num_components==6:
        return EHxyz
    elif cell.num_components==3:
        return Exyz
    else:
        raise ValueError("internal error in dft_cell_components")


# like get_dft_array, but 'zero-padded:' when the DFT cell
# does not have data for the requested component (perhaps
# because it vanishes identically by symmetry), this routine
# returns an array of the expected dimensions with all zero
# entries, instead of a rank-0 array that prints out as a
# single nonsensical floating-point number which is the
# not-very-user-friendly behavior of core pymeep.
def get_dft_array_zp(sim, cell, c, nf=0, w=None):
    array = sim.get_dft_array(cell, c, nf)
    if len(np.shape(array))==0:
        if w is None:
            _,_,_,w = sim.get_dft_array_metadata(dft_cell=cell)
        array=np.zeros(np.shape(w))
    return array


def unpack_dft_cell(sim,cell,nf=0):
    (x,y,z,w)=sim.get_dft_array_metadata(dft_cell=cell)
    cEH = dft_cell_components(cell)
    EH = [ get_dft_array_zp(sim, cell, c, nf, w) for c in cEH ]
    return x, y, z, w, cEH, EH


def unpack_dft_fields(sim,cell,nf=0):
    _,_,_,_,_,EH=unpack_dft_cell(sim,cell,nf=nf)
    return EH

###################################################
###################################################
###################################################
def eigenfield_tag(ncell, mode, freq):
    return "C{}.M{}.F{:.6f}".format(ncell,mode,freq)

def lookup_eigenfields(ncell, mode, freq, eigencache):
    tag=eigenfield_tag(ncell, mode, freq)
    if eigencache and tag in eigencache:
        return eigencache[tag]
    raise ValueError("eigenfields for tag {} not found in cache".format(tag))

def get_eigenfields(sim, dft_cells, ncell, mode, nf=0, eigencache=None):

    # look for data in cache
    cell=dft_cells[ncell]

########################################################
# this routine configures some global matplotlib settings
# (as distinct from the plot-specific settings handled
# by the plot_options dicts). if you intend to set your
# own values of these parameters, pass set_rmecParams=False
# to visualize_sim.
########################################################
def set_meep_rcParams():
    plt.rc('xtick')
    plt.rc('ytick')
    plt.rc('font', size=40)
    plt.rc('text', usetex=True)
    matplotlib.rcParams['axes.labelsize']='large'
    matplotlib.rcParams['axes.titlesize']='large'
    matplotlib.rcParams['axes.titlepad']=20

########################################################
# general-purpose default option values, customized for
# specific types of plots below, and everything settable
# via set_plot_option
########################################################
def_plot_options={ 'line_color'     : [1.0,0.0,1.0],
                   'line_width'     : 4.0,
                   'line_style'     : '-',
                   'boundary_width' : 2.0,
                   'boundary_color' : [0.0,1.0,0.0],
                   'boundary_style' : '--',
                   'z_offset'       : 0.0,
                   'zmin'           : 0.60,
                   'zmax'           : 1.00,
                   'fill_color'     : 'none',
                   'alpha'          : 1.0,
                   'cmap'           : matplotlib.cm.plasma,
                   'fontsize'       : 40,
                   'colorbar_shrink': 0.60,
                   'colorbar_pad'   : 0.04,
                   'num_contours'   : 100
                  }

#--------------------------------------------------
# options for permittivity plots
#--------------------------------------------------
def_eps_options={ **def_plot_options,
                  'cmap':matplotlib.cm.Blues,
                  'line_width': 0.00, 'colorbar_shrink':0.75,
                  'plot_method':'contourf', # or 'pcolormesh' or 'imshow'
                  'num_contours':100
                }

def_src_options={ **def_plot_options,
                  'line_width':4.0, 'line_color':[0.0,1.0,1.0],
                  'fontsize':0, 'zmin':0.0, 'zmax':0.0
                }

def_pml_options={ **def_plot_options,
                  'boundary_color':'none', 'boundary_width':0.0,
                  'fill_color': 0.75*np.ones(3), 'alpha':0.25
                }

def_flux_options={ **def_plot_options,
                   'boundary_color':[0.0,0.0,0.0], 'boundary_width':2.0,
                   'boundary_style':'--',
                   'line_color': [1.0,0.0,1.0], 'line_width':4.0
                 }

def_field_options={ **def_plot_options,
                    'line_width': 0.0, 'alpha': 0.5, 'z_offset':0.5
                  }

def_dft_options={ **def_plot_options }


def set_plot_default(option, value, type=None):
    which=       def_eps_options    if type=='eps'      \
            else def_src_options    if type=='src'      \
            else def_pml_optionss   if type=='pml'      \
            else def_flux_options   if type=='flux'     \
            else def_fields_options if type=='fields'   \
            else def_plot_options
    which[option]=value


###################################################
###################################################
###################################################
def get_text_size(fig, label, fontsize):
  r = fig.canvas.get_renderer()
  t = plt.text(0.5, 0.5, label, fontsize=fontsize)
  bb = t.get_window_extent(renderer=r)
  return bb.width, bb.height

###################################################
# visualize epsilon distribution. eps_min, eps_max
# are optional upper and lower clipping bounds.
# If we are in 3D, the permittivity is plotted using
# plot_surface.
# Otherwise, the permittivity is plotted using
# imshow (if use_imshow==True) or using pcolormesh
# (by default).
###################################################
def plot_eps(sim, eps_min=None, eps_max=None, options=None, plot3D=False, fig=None):

    options       = options if options else def_eps_options
    cmap          = options['cmap']
    edgecolors    = options['line_color']
    linewidth     = options['line_width']
    alpha         = options['alpha']
    fontsize      = options['fontsize']
    plot_method   = options['plot_method']
    num_contours  = options['num_contours']
    shading       = 'gouraud' if linewidth==0.0 else 'none'
    interpolation = 'gaussian' if linewidth==0.0 else 'none'

    eps=np.transpose(sim.get_epsilon())
    if eps_min is not None or eps_max is not None:
        eps=np.clip(eps,eps_min,eps_max)
    eps_min, eps_max=np.min(eps), np.max(eps)
    (x,y,z,w)=sim.get_array_metadata()
    extent=(min(x), max(x), min(y), max(y))

    if fig is None:
        fig = plt.gcf()
        ax  = fig.gca(projection='3d') if plot3D else fig.gca()
    else:
        ax  = fig.gca()
    cb = None
    if plot3D:  # check if 3D plot
        X, Y = np.meshgrid(x, y)
        zmin = 0.0
        zmax = max(sim.cell_size.x, sim.cell_size.y)
        Z0   = zmin + options['z_offset']*(zmax-zmin)
        img  = ax.contourf(X, Y, eps, num_contours, zdir='z', offset=Z0,
                           vmin=eps_min, vmax=eps_max, cmap=cmap, alpha=alpha)
        ax.set_zlim3d(zmin, zmax)
        ax.set_zticks([])
        pad=options['colorbar_pad']
        shrink=options['colorbar_shrink']
        cb=fig.colorbar(img, shrink=shrink, pad=pad)
    elif plot_method=='imshow':
        img = plt.imshow(np.transpose(eps), extent=extent, cmap=cmap,
                         interpolation=interpolation, alpha=alpha)
    elif plot_method=='pcolormesh':
        img = plt.pcolormesh(x,y,np.transpose(eps), cmap=cmap, shading=shading,
                             edgecolors=edgecolors, linewidth=linewidth, alpha=alpha)
    else:
        X, Y = np.meshgrid(x, y)
        img  = ax.contourf(X, Y, eps, num_contours, vmin=eps_min, vmax=eps_max,
                           cmap=cmap, alpha=alpha)

    ax.set_xlabel(r'$x$', fontsize=fontsize, labelpad=0.25*fontsize)
    ax.set_ylabel(r'$y$', fontsize=fontsize, labelpad=fontsize, rotation=0)
    ax.tick_params(axis='both', labelsize=0.75*fontsize)
    cb=cb if cb else fig.colorbar(img)
    cb.set_label(r'$\epsilon$',fontsize=1.5*fontsize,rotation=0,labelpad=0.5*fontsize)
    cb.ax.tick_params(labelsize=0.75*fontsize)
    cb.locator = ticker.MaxNLocator(nbins=5)
    cb.update_ticks()

##################################################
# plot_volume() adds a polygon representing a given
# mp.Volume to the current 2D or 3D plot.
##################################################
def plot_volume(sim, vol=None, center=None, size=None,
                options=None, plot3D=False, label=None):

    options = options if options else def_plot_options

    fig=plt.gcf()
    ax=fig.gca(projection='3d') if plot3D else fig.gca()

    if vol:
       center, size = vol.center, vol.size
    v0=np.array([center.x, center.y])
    dx,dy=np.array([0.5*size.x,0.0]), np.array([0.0,0.5*size.y])
    if plot3D:
        zmin,zmax = ax.get_zlim3d()
        z0 = zmin + options['z_offset']*(zmax-zmin)

    ##################################################
    # add polygon(s) to the plot to represent the volume
    ##################################################
    def add_to_plot(c):
        ax.add_collection3d(c,zs=z0,zdir='z') if plot3D else ax.add_collection(c)

    if size.x==0.0 or size.y==0.0:    # zero thickness, plot as line
        polygon = [ v0+dx+dy, v0-dx-dy ]
        add_to_plot( LineCollection( [polygon], colors=options['line_color'],
                                     linewidths=options['line_width'],
                                     linestyles=options['line_style']
                                   )
                   )
    else:
        if options['fill_color'] is not 'none': # first copy: faces, no edges
            polygon = np.array([v0+dx+dy, v0-dx+dy, v0-dx-dy, v0+dx-dy])
            pc=PolyCollection( [polygon], linewidths=0.0)
            pc.set_color(options['fill_color'])
            pc.set_alpha(options['alpha'])
            add_to_plot(pc)
        if options['boundary_width']>0.0: # second copy: edges, no faces
            closed_polygon = np.array([v0+dx+dy, v0-dx+dy, v0-dx-dy, v0+dx-dy, v0+dx+dy])
            lc=LineCollection([closed_polygon])
            lc.set_linestyle(options['boundary_style'])
            lc.set_linewidth(options['boundary_width'])
            lc.set_edgecolor(options['boundary_color'])
            add_to_plot(lc)

    ######################################################################
    # attempt to autodetermine text rotation and alignment
    ######################################################################
    if label:
        x0, y0, r, h, v = np.mean(ax.get_xlim()),np.mean(ax.get_ylim()), 0, 'center', 'center'
        if size.y==0.0:
            v = 'top' if center.y>y0 else 'bottom'
        elif size.x==0.0:
            r, h = (270,'left') if center.x>x0 else (90,'right')
        if plot3D:
            ax.text(center.x, center.y, z0, label, rotation=r,
                    fontsize=options['fontsize'], color=options['line_color'],
                    horizontalalignment=h, verticalalignment=v)
        else:
            ax.text(center.x, center.y, label, rotation=r,
                    fontsize=options['fontsize'], color=options['line_color'],
                    horizontalalignment=h, verticalalignment=v)


##################################################
# Plot one or more curves,
##################################################
def plot_data_curves(sim,center=None,size=None,superpose=True,
                     data=None, labels=None, options=None,
                     dmin=None, dmax=None):

    if size.x>0 and size.y>0:
        msg="plot_data_curves: expected zero-width region, got {}x{} (skipping)"
        warnings.warn(msg.format(size.x,size.y),RuntimeWarning)
        return
    if np.ndim(data[0])!=1:
        msg="plot_data_curves: expected 1D data arrays, got {} (skipping)"
        warnings.warn(msg.format(np.shape(data[0])),RuntimeWarning)
        return

    options=options if options else def_flux_options
    lw=options['line_width']
    lc=options['line_color']
    zeta_min=options['zmin']
    zeta_max=options['zmax']
    draw_baseline=(options['boundary_width']>0.0)

    kwargs=dict()
    if 'line_color' in options:
       kwargs['color']=options['line_color']
    if 'line_width' in options:
       kwargs['linewidth']=options['line_width']
    if 'line_style' in options:
       kwargs['linestyle']=options['line_style']

    # construct horizontal axis
    ii=1 if size.x==0 else 0
    hstart,hend = (center-0.5*size).__array__()[0:2], (center+0.5*size).__array__()[0:2]
    hmin,hmax=hstart[ii],hend[ii]
    haxis = np.linspace(hmin, hmax, len(data[0]))

    # if we are superposing the curves onto a simulation-geometry
    # visualization plot, construct the appropriate mapping that
    # squeezes the full vertical extent of the curve into the
    # z-axis interval [zmin, zmax]
    if superpose:
        ax=plt.gcf().gca(projection='3d')
        (zfloor,zceil)=ax.get_zlim()
        zmin=zfloor + zeta_min*(zceil-zfloor)
        zmax=zfloor + zeta_max*(zceil-zfloor)
        z0,dz=0.5*(zmax+zmin),(zmax-zmin)
        dmin=dmin if dmin else np.min(data)
        dmax=dmax if dmax else np.max(data)
        d0,dd=0.5*(dmax+dmin),(dmax-dmin)
        zs = center[1-ii]
        zdir='x' if size.x==0 else 'y'
        if draw_baseline:
            lc=LineCollection( [[hstart,hend]], colors=options['boundary_color'],
                                linewidths=options['boundary_width'],
                                linestyles=options['boundary_style']
                             )
            ax.add_collection3d(lc,zs=z0,zdir='z')

    for n in range(len(data)):
        kwargs['label']=None if not labels else labels[n]
        if superpose:
            ax.plot(haxis,z0+(data[n]-d0)*dz/dd, zs=zs, zdir=zdir, **kwargs)
        else:
            plt.plot(haxis,data[n],**kwargs)

##################################################
##################################################
##################################################
def visualize_source_distribution(sim, superpose=True, options=None):

    if not mp.am_master():
        return
    options=options if options else def_src_options

    for ns,s in enumerate(sim.sources):
        sc,ss=s.center,s.size
        J2=sum([abs2(sim.get_source_slice(c,center=sc,size=ss)) for c in Exyz])
#       M2=sum([abs2(sim.get_source_slice(c,center=sc,size=ss)) for c in Hxyz])
        if superpose==False:
            if ns==0:
                plt.ion()
                plt.figure()
                plt.title('Source regions')
            fig.subplot(len(sim.sources),1,ns+1)
            fig.title('Currents in source region {}'.format(ns))
#        plot_data_curves(sim,superpose,[J2,M2],labels=['||J||','||M||'],
#                         styles=['bo-','rs-'],center=sc,size=ssu
        plot_data_curves(sim,center=sc,size=ss, superpose=superpose,
                         data=[J2], labels=['J'], options=options)

##################################################
##################################################
##################################################
def visualize_dft_flux(sim, superpose=True, flux_cells=[],
                       options=None, nf=0):

    if not mp.am_master():
        return
    options=options if options else def_flux_options

    # first pass to get arrays of poynting flux strength for all cells
    if len(flux_cells)==0:
        flux_cells=[cell for cell in sim.dft_objects if is_flux_cell(cell)]
    flux_arrays=[]
    for cell in flux_cells:    # first pass to compute flux data
        (x,y,z,w,c,EH)=unpack_dft_cell(sim,cell,nf=nf)
        flux_arrays.append( 0.25*np.real(w*(np.conj(EH[0])*EH[3] - np.conj(EH[1])*EH[2])) )

    # second pass to plot
    for n, cell in enumerate(flux_cells):    # second pass to plot
        if superpose==False:
            if n==0:
                plt.figure()
                plt.title('Poynting flux')
            plt.subplot(len(flux_cells),1,n)
            plt.gca().set_title('Flux cell {}'.format(n))
        cn,sz=mp.get_center_and_size(cell.where)
        max_flux=np.amax(flux_arrays)
        plot_data_curves(sim, center=cn, size=sz, data=[flux_arrays[n]],
                         superpose=superpose, options=options,
                         labels=['flux through cell {}'.format(n)],
                         dmin=-max_flux,dmax=max_flux)

##################################################
##################################################
##################################################
def fc_name(c,which):
    name=mp.component_name(c)
    return name if which=='scattered' else str(name[0].upper())+str(name[1])

def field_func_array(fexpr,x,y,z,w,cEH,EH):
    if fexpr=='re(Ex)':
        return np.real(EH[0])
    if fexpr=='im(Ex)':
        return np.imag(EH[0])
    if fexpr=='re(Ey)':
        return np.real(EH[1])
    if fexpr=='im(Ey)':
        return np.imag(EH[1])
    if fexpr=='re(Ez)':
        return np.real(EH[2])
    if fexpr=='im(Ez)':
        return np.imag(EH[2])
    if fexpr=='re(Hx)':
        return np.real(EH[3])
    if fexpr=='im(Hx)':
        return np.imag(EH[3])
    if fexpr=='re(Hy)':
        return np.real(EH[4])
    if fexpr=='im(Hy)':
        return np.imag(EH[4])
    if fexpr=='re(Hz)':
        return np.real(EH[5])
    if fexpr=='im(Hz)':
        return np.imag(EH[5])
    if fexpr=='abs2(H)':
        return abs2(EH[3]) + abs2(EH[4]) + abs2(EH[5])
    if True: # fexpr=='abs2(E)':
        return abs2(EH[0]) + abs2(EH[1]) + abs2(EH[2])

##################################################
# this routine intended to be called by
# visualize_dft_fields, not directly by user
##################################################
def plot_dft_fields(sim, field_cells=[], options=None, nf=0):

    options       = options if options else def_field_options
    cmap          = options['cmap']
    alpha         = options['alpha']
    num_contours  = options['num_contours']
    fontsize      = options['fontsize']

    for ncell, cell in enumerate(field_cells):
        (x,y,z,w,cEH,EH)=unpack_dft_cell(sim,cell,nf=nf)
        X, Y = np.meshgrid(x, y)
        plt.figure()
        plt.suptitle('DFT cell {}'.format(ncell+1))
        ops=['Re', 'Im', 'Abs']
        rows,cols=len(cEH),len(ops)
        for row, c in enumerate(cEH):
            for col, op in enumerate(ops):
                plt.subplot(rows, cols, row*cols + col + 1)
                ax=plt.gca()
                ax.set_title('{}({})'.format(op,mp.component_name(c)))
                ax.set_xlabel(r'$x$', fontsize=fontsize, labelpad=fontsize)
                ax.set_ylabel(r'$y$', fontsize=fontsize, labelpad=fontsize, rotation=0)
                ax.tick_params(axis='both', labelsize=0.75*fontsize)
                data=np.real(EH[row]) if op=='Re' else np.imag(EH[row]) if op=='Im' else np.abs(EH[row])
                img = ax.contourf(X,Y,np.transpose(data),num_contours, cmap=cmap,alpha=alpha)
                plt.colorbar(img)
        #cb.set_label(r'$\epsilon$',fontsize=1.5*fontsize,rotation=0,labelpad=0.5*fontsize)
        #cb.ax.tick_params(labelsize=0.75*fontsize)

    plt.show(False)
    plt.draw()
    return 0

##################################################
##################################################
##################################################
def visualize_dft_fields(sim, superpose=True, field_cells=[],
                         field_funcs=None, z_offsets=None, options=None, nf=0):

    if not mp.am_master():
        return

    if len(field_cells)==0:
        field_cells=[cl for cl in sim.dft_objects if dft_cell_type(cl)=='fields']
    if len(field_cells)==0:
        return

    if superpose and not isinstance(plt.gcf().gca(),axes3d.Axes3D):
        warnings.warn("visualize_dft_fields: non-3D plot, can't superpose.")
        superpose=False

    if not superpose:
        return plot_dft_fields(sim, field_cells, options, nf=nf)

    # the remainder of this routine is for the superposition case

    options       = options if options else def_field_options
    cmap          = options['cmap']
    alpha         = options['alpha']
    num_contours  = options['num_contours']
    fontsize      = options['fontsize']

    if field_funcs is None:
        field_funcs, z_offsets = ['abs2(E)'], [options['z_offset']]

    for n, cell in enumerate(field_cells):

        # skip if cell is entirely contained in another field_cell
        Contained=False
        for m, parent in enumerate(field_cells):
            if m!=n and parent.where.contains(cell.where):
                Contained=True
        if Contained:
            continue

        (x,y,z,w,cEH,EH)=unpack_dft_cell(sim,cell,nf=nf)
        X, Y = np.meshgrid(x, y)
        fig = plt.gcf()
        ax  = fig.gca(projection='3d')
        (zmin,zmax)=ax.get_zlim()
        for f,z in zip(field_funcs,z_offsets):
            data=field_func_array(f,x,y,z,w,cEH,EH)
            Z0=zmin + z*(zmax-zmin)
            img = ax.contourf(X, Y, np.transpose(data), num_contours,
                              cmap=cmap, alpha=alpha, zdir='z', offset=Z0)
            pad=options['colorbar_pad']
            shrink=options['colorbar_shrink']
            cb=fig.colorbar(img, shrink=shrink, pad=pad, orientation='horizontal')
            #cb.set_label(r'$|\mathbf{E}|^2$',fontsize=1.0*fontsize,rotation=0,labelpad=0.5*fontsize)
            cb.set_label(f,fontsize=1.0*fontsize,rotation=0,labelpad=0.5*fontsize)
            cb.ax.tick_params(labelsize=0.75*fontsize)
            cb.locator = ticker.MaxNLocator(nbins=5)
            cb.update_ticks()

    plt.show(False)
    plt.draw()

##################################################
##################################################
##################################################
def visualize_sim(sim, fig=None, plot3D=None,
                  eps_min=0.0, eps_max=None,
                  eps_options=None, src_options=None,
                  pml_options=None, dft_options=None,
                  flux_options=None, field_options=None,
                  set_rcParams=True, plot_dft_data=None):

    if not mp.am_master():
        return

    # if plot3D not specified, set it automatically: false
    # if we are plotting only the geometry (at the beginning
    # of a simulation), true if we are also plotting results
    # (at the end of a simulation).
    sources_finished = sim.round_time() > sim.fields.last_source_time()
    if plot3D is None:
        plot3D=sources_finished

    ######################################################
    # create figure and set some global parameters, unless
    #  the caller asked us not to
    ######################################################
    if fig is None:
        plt.ion
        fig=plt.gcf()
        fig.clf()
        if set_rcParams:
            set_meep_rcParams()
    ax = axes3d.Axes3D(fig) if plot3D else fig.gca()

    if not plot3D:
        ax.set_aspect('equal')
        plt.tight_layout()

    ##################################################
    # plot permittivity
    ##################################################
    eps_options = eps_options if eps_options else def_eps_options
    plot_eps(sim, eps_min=eps_min, eps_max=eps_max, options=eps_options, plot3D=plot3D)

    ###################################################
    ## plot source regions and optionally source amplitudes
    ###################################################
    src_options = src_options if src_options else def_src_options
    for ns,s in enumerate(sim.sources):
        plot_volume(sim, center=s.center, size=s.size,
                    options=src_options, plot3D=plot3D,
                    label=( None if src_options['fontsize']==0
                            else 'src' + ( '\_'+str(ns) if len(sim.sources)>1 else ''))
                   )
    if src_options['zmin']!=src_options['zmax']:
        visualize_source_distribution(sim, standalone=not plot3D, options=src_options)

    ###################################################
    ## plot PML regions
    ###################################################
    def_options=def_pml_options
    if pml_options is None:
        pml_options = { **def_pml_options,
                        'fill_color' : (0.25 if plot3D else 0.75)*np.ones(3)
                      }
    if sim.boundary_layers and hasattr(sim.boundary_layers[0],'thickness'):
        dpml    = sim.boundary_layers[0].thickness
        sx, sy  = sim.cell_size.x, sim.cell_size.y
        y0, x0  = mp.Vector3(0.0, 0.5*(sy-dpml)), mp.Vector3(0.5*(sx-dpml), 0.0)
        ns, ew  = mp.Vector3(sx-2*dpml, dpml),    mp.Vector3(dpml,sy)
        centers = [ y0, -1*y0, x0, -1*x0 ]   # north, south, east, west
        sizes   = [ ns,    ns, ew,    ew ]
        for c,s in zip(centers,sizes):
            plot_volume(sim, center=c, size=s,
                        options=pml_options, plot3D=plot3D)

    ######################################################################
    # plot DFT cell regions, with labels for flux cells.
    ######################################################################
    dft_options = dft_options if dft_options else def_dft_options
    for nc, c in enumerate(sim.dft_objects):
        plot_volume(sim,center=c.regions[0].center,size=c.regions[0].size,
                    options=dft_options, plot3D=plot3D,
                    label=dft_cell_name(nc) if dft_cell_type(c)=='flux' else None)

    ###################################################
    ###################################################
    ###################################################
    if plot_dft_data is None:
        plot_dft_data=sources_finished

    if plot_dft_data:
        visualize_dft_flux(sim, superpose=True, options=flux_options)
        visualize_dft_fields(sim, superpose=True, options=field_options)

    plt.show(False)
    plt.draw()
    return fig

######################################################################
# useful options:
#
# matplotlib.rcParams['axes.titlesize']='medium'
# plt.rc('font',size=20)
######################################################################
def plot_basis(opt_prob):

    x0=opt_prob.dft_cells[-1].center
    xyzw=opt_prob.dft_cells[-1].xyzw
    x,y,z,w = xyzw[0],xyzw[1],xyzw[2],xyzw[3]
    xyz=[mp.Vector3(xx,yy,zz) for xx in x for yy in y for zz in z]
    bmatrix=np.zeros(np.shape(w))

    rows,cols=opt_prob.basis.shape
    fig, axes = plt.subplots(rows,cols)
    plt.tight_layout()
    extent=(min(x), max(x), min(y), max(y))

    usetex=matplotlib.rcParams['text.usetex']

    if usetex:
        titles=['$b_{' + str(n) + '}: ' + bn[1:] for n,bn in enumerate(opt_prob.basis.tex_names)]
    else:
        titles=['b{}: {}'.format(n,bn) for n,bn in enumerate(opt_prob.basis.names)]


    for d in range(len(titles)):
        fig.sca(axes.flat[d])
        axes.flat[d].set_title(titles[d])
        axes.flat[d].set_aspect('equal')
        print('titles[d]={}'.format(titles[d]))
        it=np.nditer(w,flags=['f_index','multi_index'])
        while not it.finished:
            n, nn = it.index, it.multi_index
            bmatrix[nn] = opt_prob.basis(xyz[n]-x0)[d]
            it.iternext()
        img=plt.imshow(np.transpose(bmatrix),extent=extent,cmap=matplotlib.cm.Blues)
        fig.colorbar(img)

    plt.show(False)
    plt.draw()


######################################################################
######################################################################
######################################################################
LogFile='/tmp/PFCServer.log'
def message(msg):
    dt=datetime.datetime.now().strftime("%T ")
    with open(LogFile,'a') as f:
        f.write("{} {}\n".format(dt,msg))
    #os.system('zenity --info --text "{} {}" &'.format(dt,msg))

PFCRequest = namedtuple('PFCRequest', 't flist')

######################################################################
######################################################################
######################################################################
class PFCServer(object):

    def __init__(self, x, y, clist, cumulative=True,
                       poll_interval=None, options=None):
        self.X, self.Y     = np.meshgrid(x, y)
        self.extent        = (min(x), max(x), min(y), max(y))
        self.clist         = clist
        self.cumulative    = cumulative
        self.poll_interval = poll_interval if poll_interval else 500
        self.options       = options if options else def_field_options
        self.num_contours  = self.options['num_contours']
        self.cmap          = self.options['cmap']

    ######################################################################
    ######################################################################
    ######################################################################
    def __call__(self, pipe):
        self.pipe  = pipe
        self.fig, self.axes = plt.subplots(1,len(self.clist))
        self.imgs, self.cbs=[],[]
        for n, ax in enumerate(self.axes):
            self.fig.sca(ax)
            ax.set_title(mp.component_name(self.clist[n]))
            img = ax.contourf(self.X,self.Y,np.zeros(np.shape(self.X)),
                              self.num_contours, cmap=self.cmap)
#            img=plt.imshow(np.zeros(np.shape(np.transpose(self.X))), cmap=self.cmap,
#                           extent=self.extent)
#                           extent=self.extent)
            cb=self.fig.colorbar(img)
            self.imgs.append(img)
            self.cbs.append(cb)

        timer=self.fig.canvas.new_timer(interval=self.poll_interval)
        timer.add_callback(self.make_callback())
        timer.start()
        plt.show()

    ######################################################################
    ######################################################################
    ######################################################################
    def make_callback(self):

        def callback():
            request = self.pipe.recv()
            if request is None:
                return False
            self.fig.suptitle('t={:.8f}'.format(request.t))
            for n,f in enumerate(request.flist):
                ax=self.axes[n]
                self.fig.sca(ax)
                plt.cla()
                ax.set_title(mp.component_name(self.clist[n]))
                vmin, vmax = np.amin(f), np.amax(f)
                if self.cumulative:
                   lvmin,lvmax = self.cbs[n].get_clim()
                   vmin, vmax  = min(vmin,lvmin), max(vmax,lvmax)
#                self.imgs[n].set_data(np.transpose(f))
#                self.imgs[n].set_clim(vmin,vmax)
                img = ax.contourf(self.X,self.Y,np.transpose(f),
                                  self.num_contours, cmap=self.cmap,
                                  vmin=vmin,vmax=vmax)
                #img.set_clim(vmin,vmax)
                self.cbs[n]=self.fig.colorbar(img, cax=self.cbs[n].ax)
                self.cbs[n].set_clim(vmin,vmax)
                self.cbs[n].draw_all()
                # self.cbs[n].update_bruteforce(img)
                #self.cbs[n].remove()
                #self.cbs[n]=self.fig.colorbar(img)
                #self.cbs[n].draw_all()
                #self.cbs[n].update_bruteforce(img)
                #self.cbs[n].
                #self.cbs[n].set_clim(vmin=vmin,vmax=vmax)
                #self.cbs[n].update
                #self.cbs[n].draw_all()
                #self.cbs[n].remove()
                #self.cbs[n]=self.fig.colorbar(img)
            self.fig.canvas.draw()
            return True

        return callback

##################################################
##################################################
##################################################
def str_to_component(s):
    for c in EHxyz:
        if mp.component_name(c)==s.lower():
            return c
    raise ValueError("unknown field component " + s)

######################################################################
######################################################################
######################################################################
def plot_field_components(sim, components, vol=None, size=None, center=origin,
                          interval=1.0, options=None):

    clist=[str_to_component(c) if isinstance(c,str) else c for c in components]

    ##################################################
    # initialization: launch server process
    ##################################################
    closure = { 'next_plot_time': interval }
    if vol is None:
        vol=mp.Volume(center=center, size=(size if size else sim.cell_size) )
    (x,y,z,w)=sim.get_array_metadata(vol=vol)

    pipe_to_server, pipe_from_client = Pipe()
    server=PFCServer(x,y,clist,options=options)
    server_process = Process(target=server, args=(pipe_from_client,), daemon=True)
    server_process.start()

    def _step_func(sim, todo):
        if todo=='step' and sim.round_time()<closure['next_plot_time']:
            return
        closure['next_plot_time']+=interval
        flist=[np.real(sim.get_array(vol=vol, component=c)) for c in clist]
        request = PFCRequest(sim.round_time(), flist)
        pipe_to_server.send(request)
        #if todo=='finish': # send 'None' request to shutdown server
        #    pipe_to_server.send(None)

    return _step_func

######################################################################
# AdjointVisualizer is an entity that knows how to create and update
# various types of visualizations of meep simulations.
######################################################################
class AdjointVisualizer(object):

    ALL_CASES=['Geometry', 'ForwardTD', 'ForwardFD', 'AdjointTD', 'AdjointFD']

    def __init__(self, label_source_regions=False, eps_options=None, src_options=None, cases=None, set_rcParams=True):
        self.cases=['Geometry', 'ForwardFD']
        self.figs, self.axes, self.cbs=[0,0,0], [0,0,0], [0,0,0]
        self.eps_options = eps_options if eps_options else def_eps_options
        self.src_options = src_options if src_options else def_src_options
        if label_source_regions:
            self.src_options['fontsize']=def_plot_options['fontsize']
        self.plot_delay  = 0.1
        if set_rcParams:
            set_meep_rcParams()

        if 'Geometry' in self.cases:
            self.figs[0]=plt.figure(1)
            self.axes[0]=self.figs[0].gca()
            self.axes[0].set_title('Geometry')

        if 'ForwardFD' in self.cases:
            self.figs[1]=plt.figure(2)
            self.axes[1]=self.figs[1].gca(projection='3d')
            self.axes[1].set_title('Forward FD fields')

#        if 'AdjointFD' in self.cases:
#            fig=plt.figure(3)
#            ax=fig.gca(projection='3d')
#            ax.set_title('Adjoint FD fields')
#            self.afd_fig=fig

        plt.show(False)
        plt.draw()
        plt.pause(self.plot_delay)


    def update(self, sim, which):

        if which=='Geometry':
            plot3D=False
            fig=plt.figure(1)
            fig.clf()
            ax=fig.gca()
            title='Geometry'

        if which=='ForwardFD':
            plot3D=True
            fig=plt.figure(2)
            fig.clf()
            #ax=fig.gca(projection='3d')
            #ax.set_axis_off()
            title='Forward FD fields, t={:.3f}'.format(sim.round_time())

#        elif which=='AdjointFD':
#            plt.figure(3)
#            plot3D=True
#            fig=self.afd_fig.clf()
#            fig.clf()
#            ax=fig.gca(projection='3d')
#            title='Adjoint FD fields, t={:.3f}'.format(sim.round_time())

        visualize_sim(sim,plot3D=plot3D,fig=fig,eps_options=self.eps_options,src_options=self.src_options)
        plt.gcf().gca().set_title(title)
        plt.draw()
        plt.pause(self.plot_delay)
