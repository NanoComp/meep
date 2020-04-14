import autograd.numpy as npa
import numpy as np
from scipy import special, signal

def smoothing_filter(x,sigma,delta,Nx,Ny):
    '''
    sigma ............. smoothing radius
    delta ............. perturbation parameter
    beta .............. thresholding parameter
    Nx ................ number of points in filter kernel in x direction
    Ny ................ number of points in filter kernel in y direction
    '''
    eta = 0.5 -  special.erf(delta/sigma)
    kernel = np.outer(signal.gaussian(Nx, sigma), signal.gaussian(Ny, sigma)) # Gaussian filter kernel
    kernel_fft = np.fft.fft2(kernel / np.sum(kernel.flatten())) # Normalize response and get freq response
    return npa.real(npa.fft.ifft2(npa.fft.fft2(x.reshape((Nx,Ny)))*kernel_fft).flatten())

def projection_filter(x,sigma,delta,beta):
    '''
    sigma ............. smoothing radius
    delta ............. perturbation parameter
    beta .............. thresholding parameter
    '''
    eta = 0.5 -  special.erf(delta/sigma)
    case1 = eta*npa.exp(-beta*(eta-x)/eta) - (eta-x)*npa.exp(-beta)
    case2 = 1 - (1-eta)*npa.exp(-beta*(x-eta)/(1-eta)) - (eta-x)*npa.exp(-beta)
    return npa.where(x < eta,case1,case2)