###################################################
# visualize.py -- various routines for visualizing
# the inputs and outputs of pymeep calculations in
# simple standardized ways
###################################################
import warnings
from datetime import datetime as dt2
from collections import namedtuple
import numpy as np
import re
from multiprocessing import Process, Pipe

import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import PolyCollection, LineCollection
import matplotlib.cm

import meep as mp

from .Objective import get_dft_cell_names
#from . get_dft_cell_names

########################################################
# the next few routines are some of my older general-purpose
# utilities for working with DFT cells; for that purpose
# they have been rendered obsolete by the DFTCell class
# in the adjoint module, but I reproduce them in this
# visualization module just so it can be independent of
# the adjoint innards.
########################################################
xHat=mp.Vector3(1.0,0.0,0.0)
yHat=mp.Vector3(0.0,1.0,0.0)
zHat=mp.Vector3(0.0,0.0,1.0)
origin=mp.Vector3()
EHTransverse=[ [mp.Ey, mp.Ez, mp.Hy, mp.Hz],
               [mp.Ez, mp.Ex, mp.Hz, mp.Hx],
               [mp.Ex, mp.Ey, mp.Hx, mp.Hy] ]
Exyz=[mp.Ex, mp.Ey, mp.Ez]
Hxyz=[mp.Hx, mp.Hy, mp.Hz]
EHxyz=Exyz + Hxyz

def abs2(z):
    return np.real(np.conj(z)*z)

def dft_cell_type(cell):
        return 'flux' if isinstance(cell, mp.simulation.DftFlux) else 'fields'

def dft_cell_name(nc):
#     name=DFTCell.get_cell_names()[nc]
    dft_cell_names=get_dft_cell_names()
    if len(dft_cell_names) <= nc:
        return 'flux{}'.format(nc)
    name=get_dft_cell_names()[nc]
    if name.endswith('_flux'):
        name=name[0:-5]
    if name.endswith('_fields'):
        name=name[0:-7]
    return name.replace('_','\_')


def is_flux_cell(cell):
    return dft_cell_type(cell)=='flux'


def add_dft_cell(sim, region, fcen, df=0, nfreq=1):
    if hasattr(region,"direction"):
        return sim.add_flux(fcen,df,nfreq,region)
    return sim.add_dft_fields(EHxyz, fcen, df,nfreq, where=region)


def flux_line(x0, y0, length, dir):
    size=length*(xHat if dir is mp.Y else yHat)
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

########################################################
# this routine configures some global matplotlib settings
# (as distinct from the plot-specific settings handled
# by the plot_options dicts). if you intend to set your
# own values of these parameters, pass set_rcParams=False
# to e.g. visualize_sim and AdjointVisualizer
########################################################
def set_meep_rcParams():
    plt.rc('xtick')
    plt.rc('ytick')
#    plt.rc('font', size=20)
    plt.rc('text', usetex=True)
    matplotlib.rcParams['axes.labelsize']='medium'
    matplotlib.rcParams['axes.titlesize']='medium'
    #matplotlib.rcParams['axes.titlepad']=20

########################################################
# general-purpose default option values, customized for
# specific types of plots below, and everything settable
# via set_plot_default
########################################################
def_plot_options={ 'line_color'           : [1.0,0.0,1.0],
                   'line_width'           : 4.0,
                   'line_style'           : '-',
                   'boundary_width'       : 2.0,
                   'boundary_color'       : [0.0,1.0,0.0],
                   'boundary_style'       : '--',
                   'zrel'                 : 0.0,
                   'zrel_min'             : 0.60,
                   'zrel_max'             : 1.00,
                   'fill_color'           : 'none',
                   'alpha'                : 1.0,
                   'cmap'                 : matplotlib.cm.plasma,
                   'fontsize'             : 25,
                   'colorbar_shrink'      : 0.60,
                   'colorbar_pad'         : 0.04,
                   'colorbar_cannibalize' : True,
                   'num_contours'         : 100,
                   'plot_delay'           : 0.1,
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

#--------------------------------------------------
# options for source region visualization (default: cyan line, no label)
#--------------------------------------------------
def_src_options={ **def_plot_options,
                  'line_width':4.0, 'line_color':[0.0,1.0,1.0],
                  'fontsize':0, 'zrel_min':0.0, 'zrel_max':0.0
                }

#--------------------------------------------------
# options for PML visualization (default: grey semitransparent blocks)
#--------------------------------------------------
def_pml_options={ **def_plot_options,
                  'boundary_color':'none', 'boundary_width':0.0,
                  'fill_color': 0.75*np.ones(3), 'alpha':0.25
                }

#--------------------------------------------------
# options for flux monitor visualization (default: magenta line with label)
#--------------------------------------------------
def_flux_options={ **def_plot_options,
                   'boundary_color':[0.0,0.0,0.0], 'boundary_width':2.0,
                   'boundary_style':'--',
                   'line_color': [1.0,0.0,1.0], 'line_width':4.0
                 }

#--------------------------------------------------
# options for dft_field cell visualization (default: green dashed border, not filled)
#--------------------------------------------------
def_field_options={ **def_plot_options,
                    'line_width': 0.0,  'alpha': 0.5,
                    'plot_method': 'contourf',
                    'zrel_min':0.4, 'zrel_max':0.6
                  }

def_dft_options={ **def_plot_options }


def set_plot_default(option, value, type=None):
    which=       def_eps_options    if type=='eps'      \
            else def_src_options    if type=='src'      \
            else def_pml_options    if type=='pml'      \
            else def_flux_options   if type=='flux'     \
            else def_field_options  if type=='fields'   \
            else def_plot_options
    which[option]=value

def plot_pause():
    plt.pause(def_eps_options['plot_delay'])

###################################################
###################################################
###################################################
def get_text_size(fig, label, fontsize):
  r = fig.canvas.get_renderer()
  t = plt.text(0.5, 0.5, label, fontsize=fontsize)
  bb = t.get_window_extent(renderer=r)
  return bb.width, bb.height

###################################################
# routine to produce proper colorbars even with subplots
# gleeped from https://joseph-long.com/writing/colorbars/
###################################################
def happy_cb(img, axes):
    divider = make_axes_locatable(axes)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return axes.figure.colorbar(img, cax=cax)

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
        Z0   = zmin + options['zrel']*(zmax-zmin)
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

    ax.set_xlabel(r'$x$', fontsize=fontsize, labelpad=0.50*fontsize)
    ax.set_ylabel(r'$y$', fontsize=fontsize, labelpad=fontsize, rotation=0)
    ax.tick_params(axis='both', labelsize=0.75*fontsize)
    cb=cb if cb else fig.colorbar(img)
    #  cb.set_label(r'$\epsilon$',fontsize=1.5*fontsize,rotation=0,labelpad=0.5*fontsize)
    cb.ax.set_xlabel(r'$\epsilon$',fontsize=1.5*fontsize,rotation=0,labelpad=0.5*fontsize)
    cb.ax.tick_params(labelsize=0.75*fontsize)
    cb.locator = ticker.MaxNLocator(nbins=5)
    cb.update_ticks()

##################################################
# plot_volume() adds a polygon representing a given
# mp.Volume to the current 2D or 3D plot, with an
# optional text label.
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
        z0 = zmin + options['zrel']*(zmax-zmin)

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
        if options['fill_color'] != 'none': # first copy: faces, no edges
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
            v = 'bottom' if center.y>y0 else 'top'
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
    #lw=options['line_width']
    lc=options['line_color']
    zrel_min=options['zrel_min']
    zrel_max=options['zrel_max']
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
        zmin=zfloor + zrel_min*(zceil-zfloor)
        zmax=zfloor + zrel_max*(zceil-zfloor)
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
            plt.fig().subplot(len(sim.sources),1,ns+1)
            plt.fig().title('Currents in source region {}'.format(ns))
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
        max_flux=np.amax([np.amax(fa) for fa in flux_arrays])
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
##################################################
##################################################
def texify(expr):
    expr=re.sub(r'([eEhH])([xyz])',r'\1_\2',expr)
    expr=re.sub(r'e_','E_',expr)
    expr=re.sub(r'H_','H_',expr)
    expr=re.sub(r'abs\((.*)\)',r'|\1|',expr)
    expr=re.sub(r'abs2\((.*)\)',r'|\1|^2',expr)

    loglike=['Re','Im']
    for s in loglike:
        expr=re.sub(s,'\textrm{'+s+'}',expr)
    return r'$'+expr+'$'

##################################################
# this routine intended to be called by
# visualize_dft_fields, not directly by user
##################################################
def plot_dft_fields(sim, field_cells=[], field_funcs=None,
                         ff_arrays=None, options=None, nf=0):

    options       = options if options else def_field_options
    cmap          = options['cmap']
    edgecolors    = options['line_color']
    linewidth     = options['line_width']
    alpha         = options['alpha']
    fontsize      = options['fontsize']
    plot_method   = options['plot_method']
    num_contours  = options['num_contours']
    shading       = 'gouraud' if linewidth==0.0 else 'none'
    interpolation = 'gaussian' if linewidth==0.0 else 'none'

    for ncell, cell in enumerate(field_cells):
        (x,y,z,w,cEH,EH)=unpack_dft_cell(sim,cell,nf=nf)
        X, Y = np.meshgrid(x, y)
        if ff_arrays is None:
            plt.figure()
            plt.suptitle('DFT cell {}'.format(ncell+1))

        if field_funcs==None:
            ops=['Re', 'Im', 'abs']
            field_funcs=[texify(op+'('+mp.component_name(c)+')') for c in cEH for op in ops]
            rows,cols=len(cEH),len(ops)
        else:
            rows,cols=1,len(field_funcs)

        def op(F,index):
            return np.real(F) if op=='Re' else np.imag(F) if op=='Im' else np.abs(F)

        for row in range(rows):
            for col in range(cols):
                nplot = row*cols + col
                data = ff_arrays[nplot] if ff_arrays else op(EH[row],ops[col])
                plt.subplot(rows, cols, nplot+1)
                ax=plt.gca()
                ax.set_title(field_funcs[nplot])
                ax.set_xlabel(r'$x$', fontsize=fontsize, labelpad=0.5*fontsize)
                ax.set_ylabel(r'$y$', fontsize=fontsize, labelpad=fontsize, rotation=0)
                ax.tick_params(axis='both', labelsize=0.75*fontsize)
                ax.set_aspect('equal')
                plt.tight_layout()
                if plot_method=='imshow':
                    img = plt.imshow(np.transpose(data), extent=(min(x), max(x), min(y), max(y)),
                                     cmap=cmap, interpolation=interpolation, alpha=alpha)
                elif plot_method=='pcolormesh':
                    img = plt.pcolormesh(x,y,np.transpose(data), cmap=cmap, shading=shading,
                                         edgecolors=edgecolors, linewidth=linewidth, alpha=alpha)
                else:
                    img = ax.contourf(X,Y,np.transpose(data),num_contours, cmap=cmap,alpha=alpha)
                #cb=plt.colorbar(img,shrink=options['colorbar_shrink'], pad=options['colorbar_pad'])
                cb=happy_cb(img,ax)
                #cb.ax.set_xlabel(ff,fontsize=1.5*fontsize,rotation=0,labelpad=0.5*fontsize)
                cb.locator = ticker.MaxNLocator(nbins=5)
                cb.update_ticks()

    plt.tight_layout()
    plt.show(False)
    plt.draw()
    return 0

##################################################
##################################################
##################################################
def visualize_dft_fields(sim, superpose=True, field_cells=[], field_funcs=None,
                         ff_arrays=None, zrels=None, options=None, nf=0):

    if not mp.am_master():
        return

    if len(field_cells)==0:
        field_cells=[cl for cl in sim.dft_objects if dft_cell_type(cl)=='fields']
        full_cells=[cell for cell in field_cells if cell.regions[0].size==sim.cell_size]
        field_cells=full_cells if full_cells else field_cells

    if len(field_cells)==0:
        return

    if superpose and not isinstance(plt.gcf().gca(),axes3d.Axes3D):
        warnings.warn("visualize_dft_fields: non-3D plot, can't superpose.")
        superpose=False

    if not superpose:
        return plot_dft_fields(sim, field_cells, field_funcs, ff_arrays, options, nf=nf)

    # the remainder of this routine is for the superposition case

    options       = options if options else def_field_options
    cmap          = options['cmap']
    alpha         = options['alpha']
    num_contours  = options['num_contours']
    fontsize      = options['fontsize']

    if field_funcs is None:
        field_funcs = ['abs2(E)']
    if zrels is None:
        zrel_min, zrel_max, nz = options['zrel_min'], options['zrel_max'], len(field_funcs)
        zrels=[0.5*(zrel_min+zrel_max)] if nz==1 else np.linspace(zrel_min,zrel_max,nz)

    for n, cell in enumerate(field_cells):
        (x,y,z,w,cEH,EH)=unpack_dft_cell(sim,cell,nf=nf)
        X, Y = np.meshgrid(x, y)
        fig = plt.gcf()
        ax  = fig.gca(projection='3d')
        (zmin,zmax)=ax.get_zlim()
        for n,(ff,zrel) in enumerate(zip(field_funcs,zrels)):
            data = ff_arrays[n] if ff_arrays else field_func_array(ff,x,y,z,w,cEH,EH)
            z0   = zmin + zrel*(zmax-zmin)
            img  = ax.contourf(X, Y, np.transpose(data), num_contours,
                               cmap=cmap, alpha=alpha, zdir='z', offset=z0)
            pad=options['colorbar_pad']
            shrink=options['colorbar_shrink']
            if options['colorbar_cannibalize']:
                cax=fig.axes[-1]
                cb=plt.colorbar(img, cax=cax)
            else:
                cb=plt.colorbar(img, shrink=shrink, pad=pad, panchor=(0.0,0.5))
            #cb.set_label(ff,fontsize=1.0*fontsize,rotation=0,labelpad=0.5*fontsize)
            cb.ax.set_xlabel(texify(ff),fontsize=1.5*fontsize,rotation=0,labelpad=0.5*fontsize)
            cb.ax.tick_params(labelsize=0.75*fontsize)
            cb.locator = ticker.MaxNLocator(nbins=5)
            cb.update_ticks()
            cb.draw_all()

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
    if src_options['zrel_min']!=src_options['zrel_max']:
        visualize_source_distribution(sim, superpose=plot3D, options=src_options)

    ###################################################
    ## plot PML regions
    ###################################################
    if sim.boundary_layers and hasattr(sim.boundary_layers[0],'thickness'):
        dpml    = sim.boundary_layers[0].thickness
        sx, sy  = sim.cell_size.x, sim.cell_size.y
        y0, x0  = mp.Vector3(0.0, 0.5*(sy-dpml)), mp.Vector3(0.5*(sx-dpml), 0.0)
        ns, ew  = mp.Vector3(sx-2*dpml, dpml),    mp.Vector3(dpml,sy)
        centers = [ y0, -1*y0, x0, -1*x0 ]   # north, south, east, west
        sizes   = [ ns,    ns, ew,    ew ]
        for c,s in zip(centers,sizes):
            plot_volume(sim, center=c, size=s, plot3D=plot3D,
                        options=pml_options if pml_options else def_pml_options)

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

    if plot_dft_data==True or plot_dft_data=='flux':
        visualize_dft_flux(sim, superpose=True, options=flux_options)

    if plot_dft_data==True:
        visualize_dft_fields(sim, superpose=True, options=field_options)

    plt.show(False)
    plt.draw()
    return fig

# ######################################################################
# # useful options:
# #
# # matplotlib.rcParams['axes.titlesize']='medium'
# # plt.rc('font',size=20)
# ######################################################################
def plot_basis(opt_prob):
    print("plot_basis temporarily not implemented")


#     x0=opt_prob.dft_cells[-1].center
#     xyzw=opt_prob.dft_cells[-1].xyzw
#     x,y,z,w = xyzw[0],xyzw[1],xyzw[2],xyzw[3]
#     xyz=[mp.Vector3(xx,yy,zz) for xx in x for yy in y for zz in z]
#     bmatrix=np.zeros(np.shape(w))
#
#     rows,cols=opt_prob.basis.shape
#     fig, axes = plt.subplots(rows,cols)
#     plt.tight_layout()
#     extent=(min(x), max(x), min(y), max(y))
#
#     usetex=matplotlib.rcParams['text.usetex']
#
#     if usetex:
#         titles=['$b_{' + str(n) + '}: ' + bn[1:] for n,bn in enumerate(opt_prob.basis.tex_names)]
#     else:
#         titles=['b{}: {}'.format(n,bn) for n,bn in enumerate(opt_prob.basis.names)]
#
#
#     for d in range(len(titles)):
#         fig.sca(axes.flat[d])
#         axes.flat[d].set_title(titles[d])
#         axes.flat[d].set_aspect('equal')
#         print('titles[d]={}'.format(titles[d]))
#         it=np.nditer(w,flags=['f_index','multi_index'])
#         while not it.finished:
#             n, nn = it.index, it.multi_index
#             bmatrix[nn] = opt_prob.basis(xyz[n]-x0)[d]
#             it.iternext()
#         img=plt.imshow(np.transpose(bmatrix),extent=extent,cmap=matplotlib.cm.Blues)
#         fig.colorbar(img)
#
#     plt.show(False)
#     plt.draw()
#
#
######################################################################
######################################################################
######################################################################
LogFile='/tmp/PFCServer.log'
def message(msg):
    tm=dt2.now().strftime("%T ")
    with open(LogFile,'a') as f:
        f.write("{} {}\n".format(tm,msg))
    #os.system('zenity --info --text "{} {}" &'.format(dt,msg))

PFCRequest = namedtuple('PFCRequest', 't flist')

######################################################################
# 'AnimateFieldEvolution' server class ###############################
######################################################################
class AFEServer(object):

    def __init__(self, x, y, components, field_options=None, poll_interval=None):
        self.X, self.Y     = np.meshgrid(x, y)
        self.extent        = (min(x), max(x), min(y), max(y))
        self.components    = components
        self.cumulative    = True; # cumulative
        self.poll_interval = poll_interval if poll_interval else 500
        self.options       = field_options if field_options else def_field_options
        self.num_contours  = self.options['num_contours']
        self.alpha         = self.options['alpha']
        self.cmap          = self.options['cmap']

    ######################################################################
    ######################################################################
    ######################################################################
    def __call__(self, pipe):
        self.pipe = pipe
        self.fig, self.axes = plt.subplots(1,len(self.components))
        if len(self.components)==1:
            self.axes=[self.axes]
        self.imgs, self.cbs = [], []
        self.fix_clim=True
        self.clim=[0.0,0.025]
        for n, ax in enumerate(self.axes):
            self.fig.sca(ax)
#            ax.set_title(mp.component_name(self.components[n]))
            img = plt.contourf(self.X,self.Y,np.zeros(np.shape(self.X)),
                               self.num_contours, cmap=self.cmap)
            ax.set_xlabel(r'$x$')
            ax.set_ylabel(r'$y$')
            cb = plt.colorbar(img) # , #orientation='horizontal',
            cb.set_clim(self.clim[0],self.clim[1])
            cb.set_ticks(np.linspace(self.clim[0],self.clim[1],3))
            cb.draw_all()
            ax.set_aspect('equal')
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
#            self.fig.suptitle('t={:.8f}'.format(request.t))
            for n,f in enumerate(request.flist):
                ax=self.axes[n]
                self.fig.sca(ax)
                plt.cla()
#                ax.set_title(mp.component_name(self.components[n]) + ' @ t={:.2f} )
                ax.set_title(r'$|E|^2, t={:.2f}$'.format(request.t))
                ax.set_xlabel(r'$x$')
                ax.set_ylabel(r'$y$')
#mp.component_name(self.components[n]) + ' @ t={:.2f} )
#                cb=self.cbs[n]
#                clim = cb.get_clim()
                f=f*f
                if not self.fix_clim:
                    self.clim[0] = min(self.clim[0],np.amin(f))
                    self.clim[1] = max(self.clim[1],np.amax(f))
                img = ax.contourf(self.X,self.Y,np.transpose(f),
                                  self.num_contours, cmap=self.cmap,
                                  vmin=self.clim[0],vmax=self.clim[1])
                ax.set_aspect('equal')
                self.cbs[n]=self.fig.colorbar(img, cax=self.cbs[n].ax)
                self.cbs[n].set_clim(self.clim[0],self.clim[1])
                self.cbs[n].set_ticks(np.linspace(self.clim[0],self.clim[1],3))
                self.cbs[n].draw_all()
            self.fig.canvas.draw()
            return True

        return callback

##################################################
##################################################
##################################################
def to_component(c):
    if c in EHxyz:
        return c
    try:
        EHxyz_names=['ex','ey','ez','hx','hy','hz']
        return EHxyz[EHxyz_names.index(c.lower())]
    except:
        return ValueError("unknown field component " + c)

######################################################################
# AFEClient (where AFE stands for 'Animate Field Evolution') is a
# class whose constructor launches an AFEServer process,
# the __call__ method of AFEClient becomes a MEEP step function,
# each time it is hod of AFEClient becomes a MEEP step function;
######################################################################
class AFEClient(object):

    ##################################################
    # the class constructor launches an AFEServer,
    # which initiates a plot and then stands by
    # awaiting updated field data
    ##################################################
    def __init__(self, sim, clist, interval=1.0,
                 vol=None, size=None, center=origin,
                 field_options=None):

        self.components=[to_component(c) for c in clist]
        self.vol=vol if vol else mp.Volume(center=center,
                                           size=(size if size else sim.cell_size) )
        (x,y,z,w)=sim.get_array_metadata(vol=self.vol)
        self.interval=interval
        self.next_plot_time=sim.round_time() + interval

        server=AFEServer(x,y,self.components,field_options=field_options)
        self.pipe_to_server, pipe_from_client = Pipe()
        Process(target=server,daemon=True, args=(pipe_from_client,)).start()

        # hack, explain and delete me
        self.__code__ = namedtuple('gna_hack',['co_argcount'])
        self.__code__.co_argcount=2


    ##################################################
    # the __call__ method has the appropriate prototype
    # to allow an instance of AFECLient to serve as a
    # meep step function; on each invocation it
    # fetches arrays of the current instantaneous values
    # of the time-domain MEEP fields and sends these to
    # the server for visualization
    ##################################################
    def __call__(self, sim, todo):
        #if todo=='finish': # send 'None' request to shutdown server
        #    pipe_to_server.send(None)
        #    return
        if sim.round_time()<self.next_plot_time:
            return
        self.next_plot_time = sim.round_time() + self.interval
        fslices=[np.real(sim.get_array(vol=self.vol, component=c)) for c in self.components]
        self.pipe_to_server.send( PFCRequest(sim.round_time(),fslices) )

#####################################################################
# This routine attempts to position the current matplotlib window on
# the user's monitor in a way that makes sense for the indexth in a
# series of plot windows. It relies on the 'screeninfo' module
# available from pypl.
#####################################################################
def set_figure_geometry(index):

    window, n = plt.get_current_fig_manager().window, index+1
######################################################################
    xmax, ymax, dx, dy=5120, 1440, 840, 670
    window.setGeometry(xmax-n*dx,ymax-dy,dx,dy)
    return
######################################################################

    try:
        # attempt to fetch size of monitor
        from screeninfo import get_monitors
        xmax, ymax=get_monitors()[0].width, get_monitors()[0].height
    except:
        warnings.warn('could not find screeninfo module; run % pip install screeninfo')
        return      # failed to fetch current monitor size

    # attempt to fetch size of current window
    if callable(getattr(window,'winfo_width',None)):    # Tkagg
        dx,dy=window.winfo_width(),window.winfo_height()
        window.wm_geometry('+{}+{}'.format(xmax-n*dx,ymax-dy))
    elif callable(getattr(window,'width',None)):        # qt{4,5}agg
        dx,dy=window.width(), window.height()
        window.setGeometry(xmax-n*dx,ymax-dy,dx,dy)
    elif getattr(window,'Size',None):                   # works for wxagg
        dx,dy=window.Size
        window.Move(xmax-n*dx, ymax-dy)
    else:
        warnings.warn('failed to auto-position matplotlib windows: unknown backend')


######################################################################
# AdjointVisualizer is an entity that knows how to create and update
# various types of MEEP simulation visualizations that are useful
# to monitor during optimization sessions.
######################################################################
class AdjointVisualizer(object):

    # TODO: add additional cases to support real-time plotting of time-domain fields
    # using AFEClient/AFEServer
    ALL_CASES=['Geometry', 'Forward', 'Adjoint']

    ##################################################
    ##################################################
    ##################################################
    def __init__(self, cases=ALL_CASES, eps_options=None, src_options=None, set_rcParams=True):

        self.cases=cases
        self.eps_options = eps_options if eps_options else def_eps_options
        self.src_options = src_options if src_options else def_src_options
        if set_rcParams:
            set_meep_rcParams()

        for case in cases:
            index = cases.index(case)
            plt.figure(index + 1)
            plot_pause()
            set_figure_geometry(index)
            if case=='Geometry':
                plt.gcf().gca().set_title(case)
            elif case=='Forward':
                plt.gcf().gca(projection='3d').set_title(case)
            else:
                plt.gcf().gca().set_title('dfdEpsilon')

        plt.show(False)
        plt.draw()
        plot_pause()

    ##################################################
    ##################################################
    ##################################################
    def update(self, sim, case, dfdEps=None):

        if case not in self.cases:
            warnings.warn("unknown case {} in AdjointVisualizer.update (ignoring)".format(case))
            return

        index = self.cases.index(case)
        fig=plt.figure(index+1)
        fig.clf()
        plot3D = (case=='Forward')

        if case in ['Forward', 'Adjoint']:
            title='{} FD fields, t={:.3f}'.format(case,sim.round_time())
        else: # case=='Geometry'
            title=case

        if case=='Adjoint':

            ff_ria= [ r'Re $\left(\frac{\partial f}{\partial\epsilon}\right)$',
                      r'Im $\left(\frac{\partial f}{\partial\epsilon}\right)$',
                      r'  $\left|\frac{\partial f}{\partial\epsilon}\right|$']
            fa_ria= [  np.real(dfdEps), np.imag(dfdEps), np,abs(dfdEps) ]

            plot_ria = False;  # -->True to plot real, imag, abs components
            if plot_ria:
               plt.suptitle(title)
               ffs, fas = ff_ria, fa_ria
            else:
               ffs, fas = [ff_ria[0] + ',' + title], [fa_ria[0]]
            visualize_dft_fields(sim,superpose=False,field_cells=[sim.dft_objects[-1]],
                                 field_funcs = ffs, ff_arrays=fas)
        else:
            visualize_sim(sim,plot3D=plot3D,fig=fig,eps_options=self.eps_options,src_options=self.src_options,
                          plot_dft_data=False if case=='Geometry' else True if case =='Forward' else 'flux')
            plt.gcf().gca().set_title(title)

        if case=='Geometry':
            plt.tight_layout()

        plt.draw()
        plot_pause()
