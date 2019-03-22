######################################################################
#             functions described by expansions in user-defined sets
#             of basis functions, plus predefined implementations of
#             some simple basis sets for common cases
######################################################################
from __future__ import division
import numpy as np
import meep as mp
import fenics

#################################################
# ideally the following should be passable to meep
# as an epsilon_function, but the typemaps in
# typemap_utils.cpp don't seem to allow it
##################################################
#class ParameterizedDielectric(object):
#
#    def __init__(self, center, basis, beta_vector):
#        self.center, self.basis = center, basis
#        self.basis.set_coefficients(beta_vector)
#
#    def __call__(p):
#       return self.basis.eval_expansion( p-self.center )

def ParameterizedDielectric(center, basis, beta_vector):
    basis.set_coefficients(beta_vector)
    def eps_func(p):
       return basis.eval_expansion( p-center )
    return eps_func


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
    return  ['1'] + [ SC[p]+'({}{})'.format('' if k is 1 else k,arg) for k in range(1,kmax+1) for p in [0,1]]

def product_name(f1,f2,tex=False):
    pname = f1 if f2=='1' else f2 if f1=='1' else f1 + '*' + f2
    return pname if not tex else '$' + pname.replace('*','') + '$'

######################################################################
# invoke python's 'abstract base class' formalism in a version-agnostic way
######################################################################
from abc import ABCMeta, abstractmethod
ABC = ABCMeta('ABC', (object,), {'__slots__': ()}) # compatible with Python 2 *and* 3:

######################################################################
# Basis is the abstract base class from which classes describing specific
# basis sets should inherit
######################################################################
class Basis(ABC):

    def __init__(self, dim):
        self.dim=dim
        self.beta_vector=np.zeros(dim)

    @property
    def dimension(self):
        return self.dim

    # derived classes must override __call__, which returns the full
    # vector of basis-function values at a given evaluation point
    @abstractmethod
    def __call__(self, p=[0.0,0.0]):
        raise NotImplementedError("derived class must implement __call__() method")

    def set_coefficients(self,beta_vector):
        self.beta_vector[:]=beta_vector[:]

    def eval_expansion(self, p=[0.0,0.0]):
        return np.dot(self.beta_vector, self.__call__(p) )

    # derived classes may override shape(), names(), texnames()
    @property
    def shape(self):
        return (self.dim,1)

    @property
    def names(self):
        return ['b{}'.format(d) for d in range(self.dim)]

    @property
    def texnames(self):
        return [r'$b_{}$'.format(d) for d in range(self.dim)]

    #########################################################j
    # project an arbitrary function f(x) onto the basis, i.e.
    # return the vector of normalized overlap integrals
    # <f|b_n> / <b_n|b_n>
    ##########################################################
    def project(self,f,xyzw):
        (x,y,z,w)=xyzw[0],xyzw[1],xyzw[2],xyzw[3]
        ffunc   = f if callable(f) else None
        fmatrix = f if not callable(f) and np.shape(f)==np.shape(w) else None
        if ffunc is None and fmatrix is None:
            raise ValueError("invalid function specification in {}.project()".format(self.__class__.__name__))
        f_dot_b, b_dot_b = 1.0j*np.zeros(self.dim), np.zeros(self.dim)
        xyz=[mp.Vector3(xx,yy,zz) for xx in x for yy in y for zz in z]
        it=np.nditer(w,flags=['f_index','multi_index'])
        while not it.finished:
            n, nn       = it.index, it.multi_index
            ww, bb, ff  = w[nn], self(xyz[n]), ffunc(xyz[n]) if ffunc else fmatrix[nn]
            f_dot_b    += (ww*ff)*bb
            b_dot_b    += ww*bb*bb
            it.iternext()
        return np.array( [ fdb/bdb for (fdb,bdb) in zip(f_dot_b,b_dot_b)] )

    #########################################################j
    ##########################################################
    def gram_matrix(self,xyzw):
        (x,y,z,w)=xyzw[0],xyzw[1],xyzw[2],xyzw[3]
        xyz=[mp.Vector3(xx,yy,zz) for xx in x for yy in y for zz in z]
        it=np.nditer(w,flags=['f_index','multi_index'])
        dim=self.dim
        gram=0.0*np.zeros([dim,dim])
        while not it.finished:
            n, nn       = it.index, it.multi_index
            bvec        = self(xyz[n])
            ww          = w[nn]
            for dr in range(dim):
                for dc in range(dim):
                    gram[dr,dc]+=ww*bvec[dr]*bvec[dc]
            it.iternext()
        return gram

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

    def __init__(self, lx, ly, kx_max=0, ky_max=0):
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

    def __call__(self, p=[0.0,0.0]):
        u=[pi/li for pi,li in zip(p,self.l)]
        return np.array([ sinusoid(nx,u[0])*sinusoid(ny,u[1]) for nx in self.nn[0] for ny in self.nn[1] ])

###################################################
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
###################################################
class FourierLegendreBasis(Basis):

    def __init__(self, radius=None, outer_radius=None, inner_radius=0.0,
                 nr_max=0, kphi_max=0):

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


    def __call__(self, p=[0.0,0.0]):
        b=np.zeros( self.nrphi[0] * self.nrphi[1] )
        r=np.sqrt(p[0]*p[0] + p[1]*p[1])
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

###################################################
## basis of 2D finite-element functions over a rectangle or disc
###################################################
class FiniteElementBasis(Basis):

    def __init__(self, N=4, lx=None, ly=None, radius=None,
                       fe_type='Lagrange', fe_order=1):

        try:
            import fenics
        except ImportError:
            raise ImportError('failed to load fenics module (needed for FiniteElementBasis)')

        if lx is not None and ly is not None:
            mesh=fenics.RectangleMesh(fenics.Point(-0.5*lx,-0.5*ly),
                                      fenics.Point(+0.5*lx,+0.5*ly),
                                      int(np.ceil(N*lx)),int(np.ceil(N*ly)))
        elif radius is not None:
            mesh=fenics.SphereMesh(fenics.Point(0.0,0.0),radius,1.0/N)
        else:
            raise ValueError('invalid parameters in FiniteElementBasis')

        self.fs  = fenics.FunctionSpace(mesh,fe_type,fe_order)
        self.f   = fenics.Function(self.fs)
        super().__init__( self.fs.dim() )

        # FIXME
        self.bfs = []
        for n in range(self.fs.dim()):
            bf=fenics.Function(self.fs)
            bf.vector()[n]=1.0
            self.bfs.append(bf)

    def __call__(self, p=[0.0,0.0]):
        return np.array( [bf(p[0],p[1]) for bf in self.bfs] )

    def set_coefficients(self, beta_vector):
        self.f.vector().set_local(beta_vector)

    def eval_expansion(self, p=[0.0,0.0]):
        return 1.0 + self.f(p[0],p[1])
