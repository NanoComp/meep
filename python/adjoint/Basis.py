######################################################################
# Basis.py -- general support for spatially-varying permittivity
#             functions described by expansions in user-defined sets
#             of basis functions, plus predefined implementations of
#             some simple basis sets for common cases
######################################################################
import numpy as np
import meep as mp

#################################################
#################################################
##################################################
def parameterized_dielectric(center, basis, beta):
    def _eps_func(p):
       return np.dot(beta,basis( (p-center).__array__()))
    return _eps_func


##################################################
# Fourier basis for a periodic function of one variable
# u (-0.5 < u < 0.5).
# For a given max frequency kmax, there are 2*kmax+1 basis functions,
# [1, sin(2*pi*u), cos(2*pi*u), sin(4*pi*u), ..., cos(2*kmax*pi*u)]
##################################################
def sinusoid(k,u):
    arg = 2.0*np.pi*np.floor((k+1)/2)*u
    return np.sin(arg) if (k%2) else np.cos(arg)

def sinusoid_names(arg,kmax):
    sinusoid_names=['1']
    for nu in range(1,kmax+1):
        sinusoid_names.append('sin({}{})'.format('' if nu==1 else str(nu),arg))
        sinusoid_names.append('cos({}{})'.format('' if nu==1 else str(nu),arg))
    return sinusoid_names


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
def plane_wave_basis(lx, ly, kx_max=0, ky_max=0):

    def _get_bf_vector(p=[0.0,0.0]):
        b=np.zeros( (2*kx_max+1)*(2*ky_max+1) )
        for kx in range(0,2*kx_max+1):
            for ky in range(0,2*ky_max+1):
                b[ kx*(2*ky_max+1) + ky ] = sinusoid(kx,p[0]/lx)*sinusoid(ky,p[1]/ly)
        return b

    return _get_bf_vector


##################################################
# basis of expansion functions f{m,n}(rho,phi) for
# a disc or annulus.
#
# f(rho,phi) = (legendre polynomial in rho) * (sinusoid in phi)
#
# f_{nr, 0}     = L_{nr}(u)
# f_{nr, 2*k-1} = L_{nr}(u) * sin(k*phi)
# f_{nr, 2*k}   = L_{nr}(u) * cos(k*phi)
#
# for nr=[0,...,nr_max], k=[0,1,...,kphi_max]
#
# Here u is a rescaled version of rho that runs over [-1:1]
# as rho runs over [inner_radius:outer_radius].
#
# The size of the basis is (nr_max+1)*(2*kphiMax+1)
##################################################
def fourier_legendre_basis(radius=None, outer_radius=None, inner_radius=0.0,
                           nr_max=0, kphi_max=0):

    num_radial=nr_max+1
    num_angular=(2*kphi_max+1)
    outer_radius=radius if outer_radius is None else outer_radius

    def _get_bf_vector(p=[0.0,0.0]):
        b=np.zeros( num_radial*num_angular )
        r=np.sqrt(p[0]*p[0] + p[1]*p[1])
        if r<inner_radius or r>outer_radius:
            b[0]=1.0;
            return b
        u = -1.0 + 2.0*(r-inner_radius)/(outer_radius-inner_radius)
        phi=np.arctan2(p[1],p[0])
        (Lm1, L)=(0.0,1.0)
        for nr in range(0,num_radial):
            for nphi in range(0,num_angular):
                b[ nr*num_angular + nphi ] = L*sinusoid(nphi,phi/(2.0*np.pi))
            (Lm2,Lm1)=(Lm1,L)
            L = ( (2*nr+1)*u*Lm1 - nr*Lm2 )/(nr+1)  # Legendre recursion
        return b

    return _get_bf_vector

def fl_basis_names(nr_max, kphi_max):
    snames=sinusoid_names('phi',kphi_max)
    flb_names=snames[:]
    for nr in range(1,nr_max+1):
        flb_names += ['L{}(r)'.format(nr) + ('' if s=='1' else s) for s in snames]
    return flb_names

#########################################################
# utility routine to project an arbitrary function f(x)
# onto a basis set {b_n(x)}, i.e. computes the N-vector
# of overlap integrals <f(x) | b_n(x) >
#########################################################
def project_basis(basis,xyzw,f):
    projection=1.0j*np.zeros(len(basis()))
    (x,y,z,w)=xyzw[0],xyzw[1],xyzw[2],xyzw[3]
    xyz=[mp.Vector3(xx,yy,zz) for xx in x for yy in y for zz in z]
    it=np.nditer(w,flags=['f_index','multi_index'])
    while not it.finished:
        n, nn = it.index, it.multi_index
        projection+=w[nn]*basis(xyz[n])*f[nn]
        it.iternext()
    return projection
