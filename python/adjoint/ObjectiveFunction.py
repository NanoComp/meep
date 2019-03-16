###############################################
# ObjectiveFunction.py --- routines for evaluating
# objective quantities and objective functions for
# adjoint-based optimization in meep
#####################################################
import datetime
from collections import namedtuple
import numpy as np
import sympy

import meep as mp

######################################################################
# various global options affecting the adjoint solver, user-tweakable
# via command-line arguments or set_adjoint_option()
######################################################################
adjoint_options={
 'dft_reltol':          1.0e-6,
 'dft_timeout':         10.0,
 'dft_interval':        0.25,
 'verbosity':           'default',
 'visualize':           False,
 'logfile':             None,
 'plot_pause':          0.01,
 'animate_components':  None,
 'animate_interval':    1.0
}

##################################################
# some convenient constants and typedefs
##################################################
xHat=mp.Vector3(1.0,0.0,0.0)
yHat=mp.Vector3(0.0,1.0,0.0)
zHat=mp.Vector3(0.0,0.0,1.0)
origin=mp.Vector3()

EHTransverse=[ [mp.Ey, mp.Ez, mp.Hy, mp.Hz],
               [mp.Ez, mp.Ex, mp.Hz, mp.Hx],
               [mp.Ex, mp.Ey, mp.Hx, mp.Hy] ]
Exyz=[mp.Ex, mp.Ey, mp.Ez]
Hxyz=[mp.Hx, mp.Hy, mp.Hz]
EHxyz=Exyz+Hxyz

# GridInfo stores the extents and resolution of the full Yee grid in a MEEP
# simulation; this is the minimal information needed to compute array metadata.
GridInfo = namedtuple('GridInfo', ['size', 'res'])

##################################################
# miscellaneous utilities
##################################################
#CU=1.0+0.0j    # complex unity, to force numbers to be treated as complex


def abs2(z):
    return np.real(np.conj(z)*z)


def unit_vector(n,N):
    return np.array([1.0 if i==n else 0.0 for i in range(N)])


# error relative to magnitude, returns value in range [0,2]
def rel_diff(a,b):
    return (      2.0    if np.isinf(a) or np.isinf(b)      \
             else 0.0    if a==0.0 and b==0.0               \
             else abs(a-b)/max(abs(a),abs(b))               \
           )


def log(msg):
    if not mp.am_master() or adjoint_options['logfile'] is None:
        return
    dt=datetime.datetime.now().strftime("%T ")
    with open(adjoint_options['logfile'],'a') as f:
        f.write("{} {}\n".format(dt,msg))

#####################################################################
# A FluxLine is essentially a 2D mp.FluxRegion with a convenient constructor
# prototype and a user-specified label to facilitate identification.
######################################################################
FluxLineT=namedtuple('FluxLine','center size direction weight name')


def FluxLine(x0,y0,length,dir,name=None):
    return FluxLineT( center=mp.Vector3(x0,y0),
                      size=length*(xHat if dir is mp.Y else yHat),
                      direction=dir, weight=1.0, name=name)

######################################################################
# DFTCell is an improved data structure for working with frequency-domain
# field components in MEEP calculations. It consolidates and replaces the
# zoo of 8 quasi-redundant DFT-related data structures in core meep (namely:
# DftFlux, DftFields, DftNear2Far, DftForce, DftLdos, DftObj, FluxRegion,
# FieldsRegion) and has a different relationship to individual MEEP
# simulations described by instances of mp.simulation().
#
# In a nutshell, the core constituents of a DFT cell are the three metadata
# fields that define the set of frequency-domain field amplitudes tabulated
# by the cell: a grid subvolume (including the associated 'xyzw' metadata),
# a set of field components, and a set of frequencies. These fields are passed
# to the DFTCell constructor and do not change over the lifetime of the DFTCell,
# which will generally encompass the lifetimes of several mp.simulation()
# instances.
#
# On the other hand, the frequency-domain field-amplitude data arrays produced by
# DFT calculations are considered less intrinsic: the DFT cell may have *no* such
# data (as when a calculation is first initiated), or may have multiple sets of
# data arrays resulting from multiple different timestepping runs. These multiple
# data sets may correspond, for example, to timestepping runs excited by different
# sources (such as the forward and adjoint sources in an adjoint-based value-and-gradient
# calculation) and/or to runs in the full and 'vacuum' versions of a geometry,
# where the latter is the 'bare' version of a geometry (with scatterers and obstacles
# removed) that one runs to tabulate incident-field profiles.
#
# Although not strictly related to DFT calculations, dft_cells describing flux-monitor
# regions also know how to compute and cache arrays of eigenmode field amplitudes.
#
# We use the following semantic conventions for choosing names for entities in the
# hierarchy of
#
#   -- At the most granular level, we have 1D, 2D, or 3D arrays (slices)
#      of frequency-domain amplitudes for a single field component at a
#      single frequency in a single calculation. For our purposes these
#      only ever arise as loop variables, which we typically call 'F' for field.
#
#   -- We often have occasion to refer to the full set of such slices
#      for all field components stored in the DFTCell, again at a single
#      frequency in a single simulation. Since this set typically includes
#      data for both E and H field components, we call it simply EH.
#      Thus, in general, EH[2] = data array for field component #2,
#      at a single frequency in a single simulation.
#
#   -- EHData refers to a collection (one-dimensional list) of EH arrays, one
#      for each component in the cell. Thus e.g. EHData[2] = array slice of
#      amplitudes for field component components[2], all at a single frequency.
#
#   -- EHCatalog refers to a collection (one-dimensional list) of EHData
#      entities, one for each frequency in the DFTCell---that is, a 2D matrix,
#      with rows corresponding to frequencies and columns corresponding to
#      components of EH arrays. Thus e.g. EHCatalog[3][2] = array slice of
#      amplitudes for component components[2] at frequency #3.
#
#   -- Arrays of eigenmode field amplitudes are named similarly with the
#      substitution "EH" -> "eh"
#
# Note: for now, all arrays are stored in memory. For large calculations with
# many DFT frequencies this may become impractical. TODO: disk caching.
######################################################################
global_dft_cell_names=[]
class DFTCell(object):

#    Why does this not work?! I would have thought DFTCell.DFTCellNames would
#    be globally accessible like a static member of a C++ class, and in fact it
#    **is** accessible from most code, but seems to go out of scope when we're
#    in the body of a method of some other class.
#    DFTCellNames=[]
#    @classmethod
#    def get_index(cls, region_name):
#        cell_name = region_name + '_flux'
#        if cell_name in cls.DFTCellNames:
#            return cls.DFTCellNames.index(cell_name)
#        raise ValueError("reference to nonexistent DFT cell {}",region_name)
    @classmethod
    def get_index(cls, region_name):
        cell_name = region_name + '_flux'
        if cell_name in global_dft_cell_names:
            return global_dft_cell_names.index(cell_name)
        raise ValueError("reference to nonexistent DFT cell {}",region_name)


    ######################################################################
    ######################################################################
    ######################################################################
    def __init__(self, grid_info=None, region=None, center=origin, size=None,
                 components=None, fcen=None, df=0, nfreq=1, name=None):

        if region is not None:
            self.center, self.size, self.region = region.center, region.size, region
        elif size is not None:
            self.center, self.size, self.region = center, size, mp.Volume(center=center, size=size)
        else:
            self.center, self.size, self.region = origin, grid_info.size, mp.Volume(center=center, size=size)

        self.nHat       = region.direction if hasattr(region,'direction') else None
        self.celltype   = 'flux' if self.nHat is not None else 'fields'  # TODO extend to other cases
        self.components =       components if components is not None                    \
                           else EHTransverse[self.nHat] if self.celltype is 'flux'      \
                           else EHxyz
        self.fcen       = fcen
        self.df         = df if nfreq>1 else 0.0
        self.nfreq      = nfreq
        self.freqs      = [fcen] if nfreq==1 else np.linspace(fcen-0.5*df, fcen+0.5*df, nfreq)

        self.sim        = None  # mp.simulation for current simulation
        self.dft_obj    = None  # meep DFT object for current simulation

        self.EH_cache   = {}    # cache of frequency-domain field data computed in previous simulations
        self.eigencache = {}    # cache of eigenmode field data to avoid redundant recalculationsq

        self.name = name
        if self.name is None:
            if hasattr(self.region,'name'):
                self.name = '{}_{}'.format(self.region.name, self.celltype)
            elif grid_info and self.size==grid_info.size:
                self.name = 'fullgrid_{}'.format(self.celltype)
            else:
                 self.name = '{}_{}'.format(self.celltype, len(global_dft_cell_names))
#                self.name = '{}_{}'.format(self.celltype, len(self.DFTCellNames))
#        self.DFTCellNames.append(self.name)
        global_dft_cell_names.append(self.name)

        # FIXME At present the 'xyzw' metadata cannot be computed until a mp.simulation / meep::fields
        #       object has been created, but in fact the metadata only depend on the GridInfo
        #       (resolution and extents of the computational lattice) and are independent
        #       of the specific material geometry and source configuration of any particular
        #       'fields' instance or simulation. In keeping with the spirit of 'DFTCell' it should
        #       be possible to compute the metadata once and for all right here before any mp.simulation()
        #       or meep::fields is created, but that will require some code refactoring. For the time being
        #       we punt on this until later, after a fields object has been created.
        self.xyzw = self.slice_dims = None

    ######################################################################
    # 'register' the cell with a MEEP timestepping simulation to request
    # computation of frequency-domain fields
    ######################################################################
    def register(self, sim):
        self.sim     = sim
        self.dft_obj =      sim.add_flux(self.fcen,self.df,self.nfreq,self.region) if self.celltype=='flux'   \
                       else sim.add_dft_fields(self.components, self.freqs[0], self.freqs[-1], self.nfreq, where=self.region)

        # take the opportunity to fill in the metadata if not done yet; #FIXME to be removed as discussed above
        if self.xyzw is None:
            self.xyzw  = sim.get_dft_array_metadata(center=self.center, size=self.size)
            self.slice_dims = np.shape(self.xyzw[3])


    ######################################################################
    # Compute an array of frequency-domain field amplitudes, i.e. a
    # frequency-domain array slice, for a single field component at a
    # single frequency in the current simulation. This is like
    # mp.get_dft_array(), but 'zero-padded:' when the low-level DFT object
    # does not have data for the requested component (perhaps because it vanishes
    # identically by symmetry), this routine returns an array of the expected
    # dimensions with all zero entries, instead of a rank-0 array that prints
    # out as a single preposterously large or small floating-point number,
    # which is the not-very-user-friendly behavior of mp.get_dft_array().
    ######################################################################
    def get_EH_slice(self, c, nf=0):
        EH = self.sim.get_dft_array(self.dft_obj, c, nf)
        return EH if np.ndim(EH)>0 else 0.0j*np.zeros(self.slice_dims)

    ######################################################################
    # Return a 1D array (list) of arrays of frequency-domain field amplitudes,
    # one for each component in this DFTCell, at a single frequency in a
    # single MEEP simulation. The simulation in question may be the present,
    # ongoing simulation (if label==None), in which case the array slices are
    # read directly from the currently active meep DFT object; or it may be a
    # previous simulation (identified by label) for which DFTCell::save_fields
    # was called at the end of timestepping.
    ######################################################################
    def get_EH_slices(self, nf=0, label=None):
        if label is None:
            return [ self.get_EH_slice(c, nf=nf) for c in self.components ]
        elif label in self.EH_cache:
            return self.EH_cache[label][nf]
        raise ValueError("DFTCell {} has no saved data for label '{}'".format(self.name, label))

    ######################################################################
    # substract incident from total fields to yield scattered fields
    ######################################################################
    def subtract_incident_fields(self, EHT, nf=0):
        EHI = self.get_EH_slices(nf=nf, label='incident')
        for nc, c in enumerate(self.components):
            EHT[nc] -= EHI[nc]

    ####################################################################
    # This routine tells the DFTCell to create and save an archive of
    # the frequency-domain array slices for the present simulation---i.e.
    # to copy the frequency-domain field data out of the sim.dft_obj
    # structure and into an appropriate data buffer in the DFTCell,
    # before the sim.dft_obj data vanish when sim is deleted and replaced
    # by a new simulation. This routine should be called after timestepping
    # is complete. The given label is used to identify the stored data
    # for purposes of future retrieval.
    ######################################################################
    def save_fields(self, label):
        if label in self.EH_cache:
            raise ValueError("DFTCell {}: data for label {} has already been saved in cache".format(self.name,label))
        self.EH_cache[label] = [self.get_EH_slices(nf=nf) for nf in range(len(self.freqs))]

    ######################################################################
    # Return a 1D array (list) of arrays of field amplitudes for all
    # tangential E,H components at a single frequency---just like
    # get_EH_slices()---except that the sliced E and H fields are the
    # fields of eigenmode #mode.
    ######################################################################
    def get_eigenfield_slices(self, mode, nf=0):

        # look for data in cache
        tag='M{}.F{}'.format(mode,nf)
        log('DFTCell {}: Getting eigenfields for tag {}...'.format(self.name,tag))
        if self.eigencache and tag in self.eigencache:
            log("...found in cache")
            return self.eigencache[tag]

        # data not in cache; compute eigenmode and populate slice arrays
        freq=self.freqs[nf]
        dir=self.nHat
        vol=mp.Volume(self.region.center,self.region.size)
        k_initial=mp.Vector3()
        eigenmode=self.sim.get_eigenmode(freq, dir, vol, mode, k_initial)

        def get_eigenslice(eigenmode, xyzw, c):
            slice=[eigenmode.amplitude(mp.Vector3(x,y,z), c)            \
                    for x in xyzw[0] for y in xyzw[1] for z in xyzw[2]
                  ]
            return np.reshape(slice,self.slice_dims)

        eh_slices=[get_eigenslice(eigenmode,self.xyzw,c) for c in self.components]

        # store in cache before returning
        if self.eigencache is not None:
            log("Adding eigenfields for tag {}".format(tag))
            self.eigencache[tag]=eh_slices

        return eh_slices

    ##################################################
    # compute an objective quantity, i.e. an eigenmode
    # coefficient or a scattered or total power.
    ##################################################
    def eval_quantity(self, qcode, mode, nf=0):

        w  = self.xyzw[3]
        EH = self.get_EH_slices(nf)
        if qcode.islower():
             self.subtract_incident_fields(EH,nf)

        if qcode in 'sS':
            return 0.25*np.real(np.sum(w*( np.conj(EH[0])*EH[3] - np.conj(EH[1])*EH[2]) ))
        elif qcode in 'PM':
            eh = self.get_eigenfield_slices(mode, nf)  # EHList of eigenmode fields
            eH = np.sum( w*(np.conj(eh[0])*EH[3] - np.conj(eh[1])*EH[2]) )
            hE = np.sum( w*(np.conj(eh[3])*EH[0] - np.conj(eh[2])*EH[1]) )
            sign=1.0 if qcode=='P' else -1.0
            return (eH + sign*hE)/8.0
        else: # TODO: support other types of objectives quantities?
            ValueError('DFTCell {}: unsupported quantity type {}'.format(self.name,qcode))

#########################################################
# a 'qrule' is a specification for how to evaluate an
# objective quantity: which DFT cell, which physical
# quantity (power flux, eigenmode coefficient, etc,
# encoded in 'code') and (if necessary) which eigenmode.
# qrules are constructed from the string names of
# objective variables like 'P2_3' or 'M1_north' or 's_0'.
#########################################################
qrule = namedtuple('qrule', 'code mode ncell')

def qname_to_qrule(qname):

    pieces=qname.split('_')
    codemode, cellstr = pieces[0], '_'.join(pieces[1:])
    ncell=int(cellstr) if cellstr.isdigit() else DFTCell.get_index(cellstr)
    if codemode.upper()=='S':
        return qrule(codemode, 0, ncell)
    elif codemode[0] in 'PM':
        if codemode[1:].isdigit() and int(codemode[1:])>0:
            return qrule(codemode[0], int(codemode[1:]) , ncell)
        raise ValueError("Objective quantity {}: invalid mode index {}".format(qname,codemode[1:]))
    raise ValueError("Objective quantity {}: unknown quantity code {}".format(qname,codemode[0]))

#########################################################
# ObjectiveFunction is a simple class for evaluating
# a scalar function f of multiple inputs {q_i}, where
# the q_i are complex-valued in general and f may be
# real- or complex-valued. The fstr input to the class
# constructor is an string expression for the function.
# Class instances store the following data:
#  -- fexpr: sympy expression constructed from fstr
#  -- qsyms: array of sympy symbols identified by sympy
#            as the objective quantities, i.e. the inputs
#            on which f depends
#  -- qnames: stringified names of the qsyms
#  -- qrules: array of 'qrule' structures encoding how
#             the objective quantities are to be computed
#             from MEEP data (specifically, from frequency-
#             domain field data stored in DFTCells)
#  -- qvals: numpy array storing most recent updates of
#            objective-quantity values
#  -- riqsyms, riqvals: the same data content as
#            qsyms and qvalues, but with each complex-valued
#            'q' quantity split up into real-valued 'r' and 'i'
#            components. We do this to facilitate symbolic
#            differentiation of non-analytic functions
#            of the objective quantities such as |q_i|^2.
#########################################################
class ObjectiveFunction(object):

    ######################################################################
    # try to create a sympy expression from the given string and determine
    # names for all input variables (objective quantities) needed to
    # evaluate it
    ######################################################################
    def __init__(self, fstr='S_0'):

        # try to parse the function string to yield a sympy expression
        try:
            fexpr = sympy.sympify(fstr)
        except:
            raise ValueError("failed to parse function {}".format(fstr))

        # qnames = names of all objective quantities (i.e. all symbols
        # qnames = names of all objective quantities (i.e. all symbols
        #          identified by sympy as quantities on which fexpr depends)
        fprime = sympy.sympify(fstr.replace('0.0','1.0'))
        self.qsyms  = sorted(fprime.free_symbols, key=str)
        self.qnames = [str(s) for s in self.qsyms]

        # qrules = 'qrules' for all objective quantities, where a 'qrule'
        #          is metadata defining how a quantity is computed
        self.qrules = [qname_to_qrule(qn) for qn in self.qnames]

        # qvals = cached values of objective quantities
        self.qvals = 0.0j*np.zeros(len(self.qnames))

        # for each (generally complex-valued) objective quantity,
        # we now introduce two real-valued symbols for the real and
        # imaginary parts, stored in riqsymbols. q2ri is a dict of
        # sympy substitutions q -> qr + I*qi that we use below to
        # recast fexpr as a function of the ri quantities. riqvals
        # is a dict of numerical values of the ri quantities used
        # later to evaluate f and its partial derivatives.
        self.riqsymbols, self.riqvals, q2ri = [], {}, {}
        for nq,(qn,qs) in enumerate(zip(self.qnames,self.qsyms)):
            rqn, iqn = 'r'+qn, 'i'+qn
            rqs, iqs = sympy.symbols( [rqn, iqn], real=True)
            q2ri[qs] = rqs + iqs*sympy.I
            self.riqvals[rqn]=self.riqvals[iqn]=0.0
            self.riqsymbols += [rqs, iqs]

        self.fexpr = fexpr.subs(q2ri)

        # expressions for partial derivatives, dfexpr[q] = \partial f / \partial q_n
        self.dfexpr=[]
        for nq in range(len(self.qnames)):
            df_drqn = sympy.diff(self.fexpr,self.riqsymbols[2*nq+0])
            df_diqn = sympy.diff(self.fexpr,self.riqsymbols[2*nq+1])
            self.dfexpr.append( df_drqn - sympy.I*df_diqn )


  ######################################################################
  ######################################################################
  ######################################################################
    def get_fq(self, DFTCells, nf=0):

        # fetch updated values for all objective quantities
        for nq, qr in enumerate(self.qrules):
            self.qvals[nq] = DFTCells[qr.ncell].eval_quantity(qr.code,qr.mode,nf)
            self.riqvals[self.riqsymbols[2*nq+0]]=np.real(self.qvals[nq])
            self.riqvals[self.riqsymbols[2*nq+1]]=np.imag(self.qvals[nq])

        # plug in objective-quantity values to get value of objective function
        fval=self.fexpr.evalf(subs=self.riqvals)
        fval=complex(fval) if fval.is_complex else float(fval)

        return np.array( [fval] + list(self.qvals) )

    # compute values of all partial derivatives \partial f / \partial q
    def get_partials(self):
        return np.array( [ df.evalf(subs=self.riqvals) for df in self.dfexpr ] )

#########################################################
# end of ObjectiveFunction class definition
#########################################################

#########################################################
#########################################################
#########################################################
class AdjointSolver(object):

    #########################################################
    #########################################################
    #########################################################
    def __init__(self, obj_func, dft_cells, basis, sim, vis=None):

        self.obj_func    = obj_func
        self.dft_cells   = dft_cells
        self.basis       = basis
        self.sim         = sim
        self.vis         = vis

        # prefetch names of outputs computed by forward and
        # adjoint solves, for use in writing log files
        self.fqnames     = ['f'] + obj_func.qnames
        self.bnames      = basis.names

    #########################################################
    #########################################################
    #########################################################
    def eval_fq(self, nf=0):
        return self.obj_func.get_fq(self.dft_cells,nf=nf)

    def eval_gradf(self, nf=0):
        gradient=0.0j*np.zeros(len(self.basis()))
        cell=self.dft_cells[-1] # design cell
        EH_forward=cell.get_EH_slices(nf,label='forward')
        EH_adjoint=cell.get_EH_slices(nf) # no label->current simulation
        f=np.sum( [EH_forward[nc]*EH_adjoint[nc]
                      for nc,c in enumerate(cell.components) if c in Exyz], 0 )
#        f=0.0j*np.ones(cell.slice_dims)
#        for nc,c in enumerate(cell.components):
#            if c not in Exyz:
#                continue
#            f+=EH_forward[nc]*EH_adjoint[nc]
        return self.basis.project(f,cell.xyzw)

    #########################################################
    #########################################################
    #########################################################
    def run_until_converged(self, case='forward'):

        last_source_time = self.sim.sources[0].src.swigobj.last_time_max()
        verbose          = (adjoint_options['verbosity'] == 'verbose')
        reltol           = adjoint_options['dft_reltol']
        max_time         = adjoint_options['dft_timeout']*last_source_time
        check_interval   = adjoint_options['dft_interval']*last_source_time
        names            = self.bnames if case=='adjoint' else self.fqnames

        # register DFT cells
        self.sim.init_sim()
        [cell.register(self.sim) for cell in self.dft_cells]

        #update_plot(self.sim, which='Geometry')
        if self.vis:
            self.vis.update(self.sim,'Geometry')

        # construct field-animation step function if requested
        step_funcs = []
        acs = adjoint_options['animate_components']
        if acs is not None:
            ivl=adjoint_options['animate_interval']
            step_funcs = [ plot_field_components(self.sim, components=acs, interval=ivl) ]

        log("Beginning {} timestepping run...".format(case))
        self.sim.run(*step_funcs, until=self.sim.fields.last_source_time())

        if self.vis:
            self.vis.update(self.sim,'ForwardFD')

        # now continue timestepping with periodic convergence checks until
        # we converge or timeout
        last_vals = np.inf*np.ones(len(names))
        max_rel_delta=np.inf
        next_check_time=self.sim.round_time()
        while True:
            self.sim.run(*step_funcs, until=next_check_time)
            next_check_time = min(next_check_time + check_interval, max_time)

            vals = self.eval_gradf() if case=='adjoint' else self.eval_fq()

            rel_delta = [rel_diff(v,last_v) for v,last_v in zip(vals,last_vals)]
            last_vals = vals
            max_rel_delta=max(rel_delta)

            log(' ** t={} MRD={} ** '.format(self.sim.round_time(),max_rel_delta))
            if verbose:
                [log('{:10s}: {:+.4e}({:.1e})'.format(n,v,e)) for n,v,e in zip(names,vals,rel_delta)]
                [log('')]

            if max_rel_delta<=reltol or self.sim.round_time()>=max_time:
                return vals


    ##############################################################
    ##############################################################
    ##############################################################
    def place_adjoint_sources(self, qweights):

        # extract temporal envelope of forward sources from existing simulation
        envelope = self.sim.sources[0].src
        freq     = envelope.frequency
        omega    = 2.0*np.pi*freq
        factor   = 2.0j*omega
        if callable(getattr(envelope, "fourier_transform", None)):
            factor /= envelope.fourier_transform(freq)

        ##################################################
        # loop over all objective quantities, adding
        # appropriately-weighted adjoint sources for each
        # quantity
        ##################################################
        self.sim.reset_meep()
        self.sim.change_sources([])
        nf=0
        for qr, qw in zip(self.obj_func.qrules, qweights):

            if qw==0.0:
                continue

            code, mode, cell=qr.code, qr.mode, self.dft_cells[qr.ncell]
            EH =      cell.get_EH_slices(nf=nf, label='forward')  if mode==0 \
                 else cell.get_eigenfield_slices(mode=mode, nf=0)

            components  = cell.components
            x, y, z, w  = cell.xyzw[0], cell.xyzw[1], cell.xyzw[2], cell.xyzw[3]
            shape       = [np.shape(q)[0] for q in [x,y,z]]

            if code in 'PM':
                sign = 1.0 if code=='P' else -1.0
                signs=[+1.0,-1.0,+1.0*sign,-1.0*sign]
                self.sim.sources+=[mp.Source(envelope, cell.components[3-nc],
                                             cell.center, cell.size,
                                             amplitude=signs[nc]*factor*qw,
                                             amp_data=np.reshape(np.conj(EH[nc]),shape)
                                            ) for nc in range(len(components))
                                  ]

        self.sim.force_complex_fields=True

    ########################################################
    ########################################################
    ########################################################
    def solve(self, need_gradient=False):

        # forward run
        fq=self.run_until_converged(case='forward')

        if self.vis:
            self.vis.update(self.sim,'ForwardFD')

        if need_gradient==False:
            return fq[0], 0

        for cell in self.dft_cells:
            cell.save_fields('forward')

        # adjoint run(s)
        qweights=self.obj_func.get_partials()
        self.place_adjoint_sources(qweights)
        gradf=self.run_until_converged(case='adjoint')
        if self.vis:
            self.vis.update(self.sim,'AdjointFD')

        return fq[0], gradf

    ########################################################
    ########################################################
    ########################################################
    def get_gradq(self, qname):
        qweights=[1.0 if qn==qname else 0.0 for qn in self.obj_func.qnames]
        self.place_adjoint_sources(qweights)
        return self.run_until_converged(case='adjoint')
