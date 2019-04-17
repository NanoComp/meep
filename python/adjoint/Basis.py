######################################################################
# Basis.py: definition of the Basis abstract base class and implementations
#           of some derived classes
######################################################################
from __future__ import division
from numbers import Number
import sys
from sympy import lambdify, Symbol
import numpy as np
import meep as mp

from . import Grid
#from Objective import Grid

######################################################################
# try to load dolfin (FENICS) module, but hold off on complaining if
# unsuccessful until someone actually tries to do something that requires it
######################################################################
try:
    import dolfin as df
except ImportError:
    pass

######################################################################
# ultimately it might be nice to use a python class with a __call__
# method as a meep epsilon_function, but the typemaps in
# typemap_utils.cpp don't seem to allow it, so the 'function that
# returns a function' approach will do for now.
######################################################################
#class ParameterizedDielectric(object):
#
#    def __init__(self, basis, beta_vector, center=mp.Vector3()):
#        self.center, self.basis = center, basis
#        self.basis.set_coefficients(beta_vector)
#
#    def __call__(p):
#       return self.basis.eval_expansion( p-self.center )

def ParameterizedDielectric(basis, beta_vector=None, center=mp.Vector3()):
    if beta_vector is not None:
        basis.set_coefficients(beta_vector)
    def eps_func(p):
       return 1.0 + basis.eval(p-center)
    return eps_func

#######################################################################
# Given (a) a scalar function f whose description may take any of several
# possible forms and (b) a grid of points in space, return a function
# GridFunc(n) that maps integer-valued grid-point indices
# to function values.
#######################################################################
class GridFunc(object):

    def __init__(self,f,grid):
        self.p=grid.points
        self.fm=self.fv=self.ff=None
        if isinstance(f,np.ndarray) and f.shape==grid.shape:
            self.fm = f.flatten()
        elif isinstance(f,Number):
            self.fv = f
        elif callable(f):
            self.ff = lambda n: f(self.p[n])
        elif isinstance(f,str):
            ffunc=lambdify( [Symbol(v) for v in 'xyz'],f)
            self.ff = lambda n:ffunc(self.p[n].x,self.p[n].y,self.p[n].z)
        else:
            raise ValueError("GridFunc: failed to construct function")

    def __call__(self, n):
        return self.fm[n] if self.fm is not None else self.fv if self.fv is not None else self.ff(n)

######################################################################
# invoke python's 'abstract base class' formalism in a version-agnostic way
######################################################################
from abc import ABCMeta, abstractmethod
ABC = ABCMeta('ABC', (object,), {'__slots__': ()}) # compatible with Python 2 and 3

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Basis is the abstract base class from which classes describing specific
# basis sets should inherit.
#----------------------------------------------------------------------
#----------------------------------------------------------------------
class Basis(ABC):

    def __init__(self, dim):
        self.dim=dim
        self.beta_vector=np.zeros(self.dim)

    @property
    def dimension(self):
        return self.dim

    ######################################################################
    # derived classes must override __call__
    ######################################################################
    @abstractmethod
    def __call__(self, p):
        raise NotImplementedError("derived class must implement __call__() method")

    ######################################################################
    # functions with boilerplate built-in implementations that may optionally
    # be overriden for speed.
    #
    # set_coefficients updates the internally-cached vector of expansion coefficients
    #
    # bvector returns the full vector of basis functions evaluated at a given point
    #
    # project computes and returns the coefficients in the basis-set projection
    # of the given function (which may be any of the objects from which we
    # know how to construct a GridFunc(), i.e. a function, a number, a string,
    # a matrix), also optionally caching them internally as in
    # set_coefficients. subtract_one may be used to ensure that the
    # function expanded in the basis is the given function minus 1.
    # (This is useful because the permittivity returned by a
    # ParameterizedDielectricFunction is defined to be 1 + the basis
    # expansion.)
    #
    # by default, the integrals needed to evaluate the expansion of an
    # arbitrary function are computed via brute-force numerical quadrature
    # in overlap() and gram_matrix(). derived classes should override these
    # methods with more efficient alternatives for particular basis sets.
    ######################################################################
    def set_coefficients(self, beta_vector, p=None):
        self.beta_vector[:]=beta_vector[:]
        return 0.0 if p is None else self(p)


    def get_bvector(self, p):
        old_beta_vector=self.beta_vector
        bvector=[self.set_coefficients(unit_vector(d,self.dim),p=p) for d in range(self.dim)]
        self.set_coefficients(old_beta_vector)
        return bvector


    def project(self,f,grid,subtract_one=False,cache=False):
        beta_vector=np.linalg.solve(self.gram_matrix(grid),
                                    self.overlap(f,grid,subtract_one=subtract_one))
        if cache:
            self.set_coefficients(beta_vector)
        return beta_vector


    ######################################################################
    # get the vector of inner products of all basis functions with an
    # arbitrary function f, i.e. the vector with components v_n \equiv <f,b_n>.
    # f may be a function that inputs a mp.Vector3 (coordinates of
    # eval point) and returns a floating-point number, or a
    # matrix of the same dimension as w, or a constant, or a string expression
    # defining a function of x,y, and z.
    ######################################################################
    def overlap(self,f,grid,subtract_one=False):
        fn = GridFunc(f,grid)
        f_dot_b = 0.0*fn(0)*np.zeros(self.dim)
        offset=1.0 if subtract_one else 0.0
        for n,(pp,ww) in enumerate(zip(grid.points,grid.weights)):
            f_dot_b += ww*(fn(n)-offset)*self.get_bvector(pp)
        return f_dot_b

    ##########################################################
    # return the matrix of basis-function inner products,
    #  gram_{ij} = <b_i | b_j>, by brute-force numerical quadrature.
    ##########################################################
    def gram_matrix_bf(self,grid):
        gm=np.zeros([self.dim,self.dim])
        for pp,ww in zip(grid.points, grid.weights):
            bvector = self.get_bvector(pp)
            for nr,br in enumerate(bvector):
                for nc,bc in enumerate(bvector):
                    gram[nr,nc]+=ww*br*bc
        return gram

    ##########################################################
    # default implementation of gram_matrix is to use
    # numerical quadrature, but subclasses should override this
    # with more efficient methods.
    ##########################################################
    def gram_matrix(self,grid=None):
        return gram_matrix_bf(self,grid)

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# end of Basis base class
#----------------------------------------------------------------------
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#----------------------------------------------------------------------
# A FiniteElementBasis is basis of 2D or 3D finite-element functions
# over a rectangle or parallelepiped, implemented using the FENICS
# finite-element package.
#----------------------------------------------------------------------
#----------------------------------------------------------------------
class FiniteElementBasis(Basis):

    ############################################################
    ############################################################
    ############################################################
    def __init__(self, size, center=mp.Vector3(),
                       num_elements=None, element_size=0.25,
                       element_type='Lagrange', element_order=1):

        if 'dolfin' not in sys.modules:
            raise ImportError('failed to load dolfin (FENICS) module, needed for FiniteElementBasis')

        vmin, vmax = center-0.5*size, center+0.5*size
        pmin, pmax = [ df.Point(v.__array__()) for v in [vmin,vmax] ]
        nn = num_elements if num_elements else [ int(np.ceil(s/element_size)) for s in size ]

        mesh    = df.RectangleMesh(pmin,pmax,nn[0],nn[1],diagonal='left') if size.z==0 else df.BoxMesh(pmin,pmax,nn[0],nn[1],nn[2])
        fs      = df.FunctionSpace(mesh,element_type,element_order)
        self.f  = df.Function(fs)

        # precompute and save gram matrix
        u       = df.TrialFunction(fs)
        v       = df.TestFunction(fs)
        A       = df.assemble(u*v*df.dx)
        self.gm = A.array()

        super().__init__( fs.dim() )


    ############################################################
    ############################################################
    ############################################################
    def __call__(self, p):
        return self.f(p)

    def set_coefficients(self, beta_vector):
        self.f.vector().set_local(beta_vector)

    def gram_matrix(self,grid=None):
        return self.gm

    ############################################################
    ############################################################
    ############################################################
    def get_bvector(self, p):
        fs      = self.f.function_space()
        mesh    = fs.mesh()
        dm      = fs.dofmap()
        el      = fs.element()

        # determine the mesh cell in which p lies
        pnt     = df.Point(p.__array__())
        cidx    = mesh.bounding_box_tree().compute_first_entity_collision(pnt)
        cell    = df.Cell(mesh,cidx)
        cdofs   = cell.get_vertex_coordinates()

        # get global indices and values (at p) of this cell's basis functions
        bfidxs  = dm.cell_dofs(cidx)
        bfvals  = el.evaluate_basis_all(p.__array__(),cdofs,cell.orientation())
        bvec    = np.zeros(fs.dim())
        for i,v in zip(bfidxs, bfvals):
            bvec[i]=v
        return bvec

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# end of FiniteElementBasis class
#----------------------------------------------------------------------
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Poor man's version of FiniteElementBasis to enable 2D adjoint
# modeling even if the dolfin / FENICS module is not available.
# Equivalent to FiniteElementBasis for element_type='Lagrange' and
# element_order=1.
#
# How it works: The sides of an (lx x ly) rectangle are subdivided into
# num_elements[0] x num_elements[1] segments (or into segments of
# length element_size). To each node in the resulting (NX+1)*(NY+1)
# grid of nodes we assign a single basis function. The basis function
# associated with node #n is only supported on those triangles
# having node #n as a vertex; it takes the value 1 at node #n and
# falls linearly to zero at the neighboring nodes. Indexing: The
# node with grid coordinates (nx,ny) (0 \le ni \le N_i) is assigned '
# index nx*(NY+1) + ny, following the conventional MEEP scheme for grid indexing.
#----------------------------------------------------------------------
#----------------------------------------------------------------------
class SimpleFiniteElementBasis(Basis):

    def __init__(self, size=mp.Vector3(), center=mp.Vector3(),
                       num_elements=None, element_size=0.25):
        if size.z>0.0:
            mp.abort('SimpleFiniteElementBasis not implemented for 3D problems')
        self.p0 = center
        self.l  = [size.x,size.y]
        self.N  = num_elements if num_elements is not None else [int(np.ceil(ll/element_size)) for ll in self.l]
        self.d  = [ll/(1.0*NN) for ll,NN in zip(self.l,self.N)]
        super().__init__( (self.N[0]+1)*(self.N[1]+1) )

    def in_grid(self, n, ofs=[0,0]):
        return np.all([nn+oo in range(0,NN+1) for(nn,oo,NN) in zip(n,ofs,self.N)])

    # scalar index of basis function associated with node n + optional offset
    def bindex(self, n, ofs=[0,0]):
        return -1 if not self.in_grid(n,ofs) else (n[0]+ofs[0])*(self.N[1]+1) + (n[1]+ofs[1])

    @property
    def shape(self):
        return (self.N[0]+1,self.N[1]+1)

    ##############################################################
    # on input, p[0,1] are the x,y coordinates of an evaluation
    # point in the grid. The return value is a list of
    # (bindex,bvalue) pairs giving the index and value of all
    # basis functions supported at p.
    ##############################################################
    def contributors(self, p):
        p            = p - self.p0
        pshift       = [ pp + 0.5*ll for pp,ll in zip(p,self.l) ]
        node         = [ int(np.floor(pp/dd)) for pp,dd in zip(pshift,self.d) ]
        xi           = [ (pp-nn*dd)/(1.0*dd) for (pp,nn,dd) in zip(pshift,node,self.d) ]
        xisum, lower = xi[0]+xi[1], xi[0] <= (1.0-xi[1])
        indices      = [self.bindex(node,ofs) for ofs in [(1,0), (0,1), (0,0) if lower else (1,1)]]
        vals         = [xi[0], xi[1], 1.0-xisum] if lower else [1.0-xi[1],1.0-xi[0],xisum-1.0]
        return [ (i,v) for i,v in zip(indices,vals) if i!= -1 ]

    def __call__(self, p):
        return sum( [ self.beta_vector[idx]*val for idx,val in self.contributors(p) ] )

    def get_bvector(self, p):
        bvector=np.zeros(self.dim)
        for idx, val in self.contributors(p):
            bvector[idx]=val
        return bvector

    def gram_matrix(self, grid=None):
        diag,off_diag=np.array([1.0/2.0,1.0/12.0]) * (self.d[0]*self.d[1])
        gm = diag*np.identity(self.dim)
        offsets=[(dx,dy) for dx in [-1,0,1] for dy in [-1,0,1] if dx!=dy ]
        for (i,n) in enumerate([(nx,ny) for nx in range(0,self.N[0]+1) for ny in range(0,self.N[1]+1)]):
            for j in [ self.bindex(n,ofs) for ofs in offsets if self.in_grid(n,ofs) ]:
                gm[i,j]=gm[j,i]=(diag if i==j else off_diag)
        return gm

######################################################################
# The remainder of this file provides implementations of spectral
# basis sets, e.g. plane waves for rectilinear domains.
######################################################################

##################################################
# Fourier basis for a periodic function of one variable
# u (-0.5 < u < 0.5).
# For a given max frequency kmax, there are 2*kmax+1 basis functions,
# [1, sin(2*pi*u), cos(2*pi*u), sin(4*pi*u), ..., cos(2*kmax*pi*u)]
##################################################
def sinusoid(k,u):
    arg = 2.0*np.pi*np.floor((k+1)/2)*u
    return np.sin(arg) if (k%2) else np.cos(arg)

def sinusoid_names(arg,kmax,tex=False):
    c='\\' if tex else ''
    SC=[c+'sin',c+'cos']
    return  ['1'] + [ SC[p]+'({}{})'.format('' if k==1 else k,arg) for k in range(1,kmax+1) for p in [0,1]]

def product_name(f1,f2,tex=False):
    pname = f1 if f2=='1' else f2 if f1=='1' else f1 + '*' + f2
    return pname if not tex else '$' + pname.replace('*','') + '$'

##################################################
# Plane-wave basis for a rectangular region.
#
# f_{nx,ny} = S_{nx}(x/lx) * S_{ny}(y/ly)
#
# S_{0}(u)    = 1
# S_{2k-1}(u) = sin(2*pi*k*u)
# S_{2k}  (u) = cos(2*pi*k*u)
#
# The size of the basis is (2*kx_max+1)*(2*ky_max+1)
##################################################
class PlaneWaveBasis(Basis):

    def __init__(self, lx, ly, center=mp.Vector3(), kx_max=0, ky_max=0):
        self.p0         = np.array([center.x,center.y])
        self.l          = [lx,ly]
        self.kmax       = [kx_max, ky_max]
        self.nn         = range(2*kx_max+1), range(2*ky_max+1)
        fxnames,fynames = [sinusoid_names(arg,kmax) for (arg,kmax) in zip('xy',self.kmax)]
        self.fnames     = [product_name(fxn,fyn) for fxn in fxnames for fyn in fynames]
        fxnames,fynames = [sinusoid_names(arg,kmax,tex=True) for (arg,kmax) in zip(['\\overline{x}','\\overline{y}'],self.kmax)]
        self.tex_fnames = [product_name(fx,fy,tex=True) for fx in fxnames for fy in fynames]
        super().__init__(len(self.fnames))

    @property
    def shape(self):
        return ( len(self.nn[0]) , len(self.nn[1]) )

    @property
    def names(self):
        return self.fnames

    @property
    def tex_names(self):
        return self.tex_fnames

    def __call__(self, p=None):
        p = np.array([0.0,0.0]) if p is None else p-self.p0
        u = [pi/li for pi,li in zip(p,self.l)]
        return np.array([ sinusoid(nx,u[0])*sinusoid(ny,u[1]) for nx in self.nn[0] for ny in self.nn[1] ])

######################################################################
## basis of expansion functions f{m,n}(r,phi) for
## a disc or annulus.
##
## f(r,phi) = (legendre polynomial in ur) * (sinusoid in phi)
##
## f_{nr, 0}     = L_{nr}(ur)
## f_{nr, 2*k-1} = L_{nr}(ur) * sin(k*phi)
## f_{nr, 2*k}   = L_{nr}(ur) * cos(k*phi)
##
## for nr=[0,...,nr_max], k=[0,1,...,kphi_max]
##
## Here ur is a rescaled version of r that runs over [-1:1]
## as rho runs over [inner_radius:outer_radius].
##
## The size of the basis is (nr_max+1)*(2*kphiMax+1)
######################################################################
class FourierLegendreBasis(Basis):

    def __init__(self, center=mp.Vector3(), radius=None, outer_radius=None, inner_radius=0.0,
                 nr_max=0, kphi_max=0):

        self.p0           = np.array([center.x,center.y])
        self.nrphi        = (nr_max+1, 2*kphi_max+1)
        self.outer_radius = radius if outer_radius is None else outer_radius
        self.inner_radius = 0.0 if inner_radius is None else inner_radius
        self.radial_span  = self.outer_radius - self.inner_radius
        frnames           = ['1'] + ['P_{}(u)'.format(n) for n in range(1,nr_max+1)]
        fpnames           = sinusoid_names('theta',kphi_max)
        self.fnames       = [product_name(fr,fp) for fr in frnames for fp in fpnames]
        frnames           = ['1'] + ['P_{'+str(n)+'}(\overline{r})' for n in range(1,nr_max+1)]
        fpnames           = sinusoid_names('\\theta',kphi_max,tex=True)
        self.tex_fnames   = [product_name(fr,fp,tex=True) for fr in frnames for fp in fpnames]
        super().__init__(len(self.fnames))

    @property
    def shape(self):
        return (self.nrphi[0],self.nrphi[1])

    @property
    def tex_names(self):
        return self.tex_fnames

    @property
    def names(self):
        return self.fnames


    def __call__(self, p=None):
        p = np.zeros(2) if p is None else p-self.p0
        b = np.zeros( self.nrphi[0] * self.nrphi[1] )
        r = np.sqrt(p[0]*p[0] + p[1]*p[1])
        if r<self.inner_radius or r>self.outer_radius:
            b[0]=1.0;
            return b
        ur = -1.0 + 2.0*(r-self.inner_radius)/self.radial_span
        uphi = np.arctan2(p[1],p[0])/(2.0*np.pi)
        (Lm1, L)=(0.0,1.0) # initialize Legendre recursion
        for nr in range(self.nrphi[0]):
            for nphi in range(self.nrphi[1]):
                b[ nr*self.nrphi[1] + nphi ] = L*sinusoid(nphi,uphi)
            (Lm2,Lm1)=(Lm1,L)
            L = ( (2*nr+1)*ur*Lm1 - nr*Lm2 )/(nr+1)  # Legendre recursion
        return b
