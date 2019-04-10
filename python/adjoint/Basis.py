######################################################################
# Basis.py:
#             of basis functions, plus predefined implementations of
#             some simple basis sets for common cases
######################################################################
from __future__ import division
from numbers import Number
from sympy import lambdify, Symbol
import numpy as np
import meep as mp
#import fenics

######################################################################
# ideally it should be possible to use a python class with a __call__
# method as a meep epsilon_function, but the typemaps in
# typemap_utils.cpp don't seem to allow it, so the 'function that
# returns a function' approach will do for now.
######################################################################
#class ParameterizedDielectric(object):
#
#    def __init__(self, center, basis, beta_vector):
#        self.center, self.basis = center, basis
#        self.basis.set_coefficients(beta_vector)
#
#    def __call__(p):
#       return self.basis.eval_expansion( p-self.center )

def ParameterizedDielectric(center, basis, beta_vector=None):
    if beta_vector is not None:
        basis.set_coefficients(beta_vector)
    def eps_func(p):
       return 1.0 + basis.eval( p-center )
    return eps_func

#######################################################################
# Given (1) a scalar function f whose description may take any of several
# possible forms and (2) a grid of points {p}, return a function
# GridFunc(n) that inputs an integer n the function of f at the nth grid point.
#######################################################################
class GridFunc(object):

    def __init__(self,f,p,shape):
        self.p = p
        self.fm=self.fv=self.ff=None
        if isinstance(f,np.ndarray) and np.shape(f)==shape:
            self.fm = f.flatten()
        elif isinstance(f,Number):
            self.fv = f
        elif callable(f):
            self.ff=lambda n: f(p[n])
        elif isinstance(f,str):
            ffunc=lambdify( [Symbol(v) for v in 'xyz'],f)
            self.ff = lambda n:ffunc(p[n][0],p[n][1],p[n][2])
        else:
            raise ValueError("GridFunc: failed to construct function")

    def __call__(self, n):
        return self.fm[n] if self.fm is not None else self.fv if self.fv is not None else self.ff(n)

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

######################################################################
# invoke python's 'abstract base class' formalism in a version-agnostic way
######################################################################
from abc import ABCMeta, abstractmethod
ABC = ABCMeta('ABC', (object,), {'__slots__': ()}) # compatible with Python 2 *and* 3:

######################################################################
# Basis is the abstract base class from which classes describing specific
# basis sets should inherit.
######################################################################
class Basis(ABC):

    def __init__(self, dim):
        self.dim=dim
        self.beta_vector=np.zeros(self.dim)

    @property
    def dimension(self):
        return self.dim

    ######################################################################
    # derived classes must override __call__, which returns the full
    # vector of basis-function values at a given evaluation point
    ######################################################################
    @abstractmethod
    def __call__(self, p=[0.0,0.0]):
        raise NotImplementedError("derived class must implement __call__() method")

    ######################################################################
    # functions with boilerplate built-in implementations that may be
    # overriden for speed.
    #
    # eval() returns the basis expansion (\sum \beta_n b_n(x)) for a given
    # vector of expansion coefficients at a given point.
    #
    # tabulate returns a numpy array giving the basis-function expansion
    # at each point in a given grid of points.
    #
    # set_coefficients caches the given vector of expansion coefficients
    # for use on subsequent calls to eval().
    #
    # expand computes and returns the coefficients in the basis-set expansion
    # of the given function (which may be any of the objects from which we
    # know how to construct a GridFunc(), i.e. a function, a number, a string,
    # a matrix), also optionally caching them internally as in
    # set_coefficients. subtract_one may be used to ensure that the
    # function expanded in the basis is the given function minus 1.
    # (this is useful because the permittivity returned by a
    # ParameterizedDielectricFunction is defined to be 1 + the basis expansion.
    # else the internally cached coefficients.
    #
    # by default, the integrals needed to evaluate the expansion of an
    # arbitrary function are computed via brute-force numerical quadrature
    # in overlap() and gram_matrix(). derived classes should override these
    # methods with more efficient alternatives for particular basis sets.
    ######################################################################
    def eval(self, p=[0.0,0.0], beta_vector=None):
        return np.dot(self.beta_vector if beta_vector is None else beta_vector, self.__call__(p) )

    def tabulate(self, xyzw, beta_vector=None):
        grid   = np.zeros(np.shape(xyzw[3]))
        garray = grid.flatten()
        parray = [mp.Vector3(xx,yy,zz) for xx in xyzw[0] for yy in xyzw[1] for zz in xyzw[2]]
        for n,p in enumerate (parray):
            garray[n]
        return np.reshape([self.eval(p,beta_vector) for p in parray],np,np.shape(xyzw[3]))

    def set_coefficients(self,beta_vector):
        self.beta_vector[:]=beta_vector[:]

    def expand(self,f,xyzw,cache=False,subtract_one=False):
        beta_vector=np.linalg.solve(self.gram_matrix(xyzw),self.overlap(f,xyzw,subtract_one))
        if cache:
            self.beta_vector=beta_vector
        return beta_vector

    ######################################################################
    # get the vector of inner products of all basis functions with an
    # arbitrary function f, i.e. the vector with components v_n \equiv <f,b_n>.
    # xyzw = array metadata for the region over which to integrate.
    # f may be a function that inputs a mp.Vector3 (coordinates of
    # eval point) and returns a floating-point number, or a
    # matrix of the same dimension as w, or a constant, or a string expression
    # defining a function of x,y, and z.
    ######################################################################
    def overlap(self,f,xyzw, subtract_one=False):
        p,w = [mp.Vector3(xx,yy,zz) for xx in xyzw[0] for yy in xyzw[1] for zz in xyzw[2]], xyzw[3].flatten()
        fn  = GridFunc(f,p,np.shape(xyzw[-1]))
        f_dot_b = 0.0*fn(0)*np.zeros(self.dim)
        offset=1.0 if subtract_one else 0.0
        for n,(pp,ww) in enumerate(zip(p,w)):
            f_dot_b += ww*self(pp)*(fn(n) - offset)
        return f_dot_b

    ##########################################################
    # return the matrix of basis-function inner products,
    #  gram_{ij} = <b_i | b_j>, by brute-force numerical quadrature.
    ##########################################################
    def gram_matrix_bf(self,xyzw):
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

    ##########################################################
    # default implementation of gram_matrix is to use
    # numerical quadrature, but subclasses should override this
    # with more efficient methods.
    ##########################################################
    def gram_matrix(self,xyzw):
        return self.gram_matrix_bf(xyzw)

    ######################################################################
    # derived classes may override shape(), names(), texnames(), which are
    # used for plotting the basis
    ######################################################################
    @property
    def shape(self):
        return (self.dim,1)

    @property
    def names(self):
        return ['b{}'.format(d) for d in range(self.dim)]

    @property
    def tex_names(self):
        return [r'$b_{}$'.format(d) for d in range(self.dim)]

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
## basis of 2D finite-element functions over a rectangle or disc,
## implemented using the FENICS finite-element package; this is
## temporarily shelved pending an assessment of how limiting it
## may be to have a FENICS dependency. Instead, the routine
## below implements a baby version supporting only first-order
## Lagrange functions on a rectangular domain.
###################################################
#class FiniteElementBasis(Basis):
#
#    def __init__(self, N=4, lx=None, ly=None, radius=None,
#                       fe_type='Lagrange', fe_order=1):
#
#        try:
#            import fenics
#        except ImportError:
#            raise ImportError('failed to load fenics module (needed for FiniteElementBasis)')
#
#        if lx is not None and ly is not None:
#            mesh=fenics.RectangleMesh(fenics.Point(-0.5*lx,-0.5*ly),
#                                      fenics.Point(+0.5*lx,+0.5*ly),
#                                      int(np.ceil(N*lx)),int(np.ceil(N*ly)))
#        elif radius is not None:
#            mesh=fenics.SphereMesh(fenics.Point(0.0,0.0),radius,1.0/N)
#        else:
#            raise ValueError('invalid parameters in FiniteElementBasis')
#
#        self.fs  = fenics.FunctionSpace(mesh,fe_type,fe_order)
#        self.f   = fenics.Function(self.fs)
#        super().__init__( self.fs.dim() )
#
#        # FIXME
#        self.bfs = []
#        for n in range(self.fs.dim()):
#            bf=fenics.Function(self.fs)
#            bf.vector()[n]=1.0
#            self.bfs.append(bf)
#
#    def __call__(self, p=[0.0,0.0]):
#        if self.fs.mesh().bounding_box_tree().compute_collisions(fenics.Point(p[0],p[1])) == []:
#            return 0.0
#        return np.array( [bf(p[0],p[1]) for bf in self.bfs] )
#
#    def set_coefficients(self, beta_vector):
#        self.f.vector().set_local(beta_vector)
#
#    def eval_expansion(self, p=[0.0,0.0]):
#        if self.fs.mesh().bounding_box_tree().compute_collisions(fenics.Point(p[0],p[1])) == []:
#            return 0.0
#        return 1.0 + self.f(p[0],p[1])

######################################################################
## finite-element basis for a rectangle.
## density is the number of basis functions per unit length.
##
## (This is a somewhat simplistic implementation that supports only
## first-order elements over triangles, possibly to be replaced
## eventually by something more sophisticated from the many available
## choices of python FEM codes.)
##
## How it works: The sides of an (lx x ly) rectangle are divided in to
## segments at a resolution of `density` functions per unit length.
## To each node in the resulting (NX+1)*(NY+1) grid of nodes we assign
## a single basis function. The basis function associated with node #n
## is only supported on those triangles naving node #n as a vertex;
## it takes the value 1 at node #n and falls linearly to zero at
## the neighboring nodes. Indexing: The node with grid coordinates
## (nx,ny) (0 \le ni \le N_i) is assigned index nx*(NY+1) + ny,
## following the convention MEEP scheme for grid indexing.
######################################################################
class FiniteElementBasis(Basis):

    def __init__(self, lx, ly, density=4):
        self.l=[lx,ly]
        self.N=[int(np.ceil(density*ll)) for ll in self.l]
        self.d=[ll/(1.0*NN) for ll,NN in zip(self.l,self.N)]
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
        pshift       = [ pp + 0.5*ll for pp,ll in zip(p,self.l) ]
        node         = [ int(np.floor(pp/dd)) for pp,dd in zip(pshift,self.d) ]
        xi           = [ (pp-nn*dd)/(1.0*dd) for (pp,nn,dd) in zip(pshift,node,self.d) ]
        xisum, lower = xi[0]+xi[1], xi[0] <= (1.0-xi[1])
        indices      = [self.bindex(node,ofs) for ofs in [(1,0), (0,1), (0,0) if lower else (1,1)]]
        vals         = [xi[0], xi[1], 1.0-xisum] if lower else [1.0-xi[1],1.0-xi[0],xisum-1.0]
        return [ (i,v) for i,v in zip(indices,vals) if i!= -1 ]

    def __call__(self, p=[0.0,0.0]):
        b=np.zeros(self.dim)
        for idx,val in self.contributors(p):
           b[idx]=val
        return b

    def eval(self, p=[0.0,0.0], beta_vector=None):
        bv = beta_vector if beta_vector else self.beta_vector
        return sum( [ bv[idx]*val for idx,val in self.contributors(p) ] )

    def gram_matrix(self, xyzw):
        diag,off_diag=np.array([1.0/2.0,1.0/12.0]) * (self.d[0]*self.d[1])
        gm = diag*np.identity(self.dim)
        offsets=[(dx,dy) for dx in [-1,0,1] for dy in [-1,0,1] if dx!=dy ]
        for (i,n) in enumerate([(nx,ny) for nx in range(0,self.N[0]+1) for ny in range(0,self.N[1]+1)]):
            for j in [ self.bindex(n,ofs) for ofs in offsets if self.in_grid(n,ofs) ]:
                gm[i,j]=gm[j,i]=(diag if i==j else off_diag)
        return gm
