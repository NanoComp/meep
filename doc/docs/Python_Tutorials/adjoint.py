##################################################
# adjoint.py --- a collection of pymeep routines
# for adjoint-based sensitivity analysis of
# material geometries
##################################################
import datetime;

e_components=[mp.Ex, mp.Ey, mp.Ez]
eh_components=[mp.Ex, mp.Ey, mp.Ez, mp.Hx, mp.Hy, mp.Hz]

##################################################
# the four field components stored in a dft_flux
# object with the given normal direction
##################################################
def get_flux_components(dir):
    if dir==mp.X:
      return [mp.Ey, mp.Ez, mp.Hz, mp.Hy]
    elif dir==mp.Y:
      return [mp.Ez, mp.Ex, mp.Hx, mp.Hz]
    else: # dir=mp.Z
      return [mp.Ex, mp.Ey, mp.Hy, mp.Hx]

def get_dft_flux_fields(sim,flux):
    components=get_flux_components(flux.normal_direction)
    fields=[ sim.get_dft_array(flux, c, 0) for c in components ]
    f1=fields[1];
    for nf in range(0,len(fields)):
        if len(np.shape(fields[nf]))==0:
            fields[nf]=np.zeros(np.shape(f1))
    return fields

#########################################################
# phase of complex number in degrees
#########################################################
def phase(z):
    return np.angle(z)*180.0/np.pi

#########################################################
# Given a list of dft_flux objects, continue timestepping
# until all fluxes have stopped changing (to the specified
# absolute tolerance) or the time limit is reached.
#########################################################
def run_until_flux_converged(sim, dft_fluxes, check_interval=20.0,
                             tol=1.0e-6, max_time=2000):

    sim.run(until=sim.fields.last_source_time())
    next_check_time = sim.round_time() + check_interval
    last_S = np.array([mp.get_fluxes(flux)[0] for flux in dft_fluxes])
    done=False;
    while not done:
        sim.fields.step()
        if sim.round_time() < next_check_time:
            continue
        next_check_time += check_interval

        S = np.array([mp.get_fluxes(flux)[0] for flux in dft_fluxes])
        delta=np.abs(S-last_S)
        rel_delta=[delta[n]/(1 if S[n]==0 else np.abs(S[n])) for n in range(0,len(S))]
        last_S = S
        
        if mp.am_master():
          dt=datetime.datetime.now().strftime("%D::%T: ");
          print dt, "{}: ".format(sim.round_time()), "S=",S,", dS=", rel_delta
        done = ( (max(rel_delta)<tol) or (sim.round_time() > max_time) )

#########################################################
#########################################################
#########################################################
def compute_gradient(sim, region, forward_fields, adjoint_fields, basis):
    gradient=0.0j*np.zeros(len(basis(mp.Vector3())))
    xyzw=sim.get_dft_array_metadata(region)
    for nc in range(0,3):
        (E0,EA)=(forward_fields[nc],adjoint_fields[nc])
        if len(np.shape(E0))==0 or (np.shape(EA)!=np.shape(E0)):
            continue
        it=np.nditer(E0,flags=['multi_index'])
        while not it.finished:
            ii=it.multi_index
            (p,w)=(mp.Vector3(*tuple(xyzw[ii][0:3])), xyzw[ii][3])
            gradient += w*E0[ii]*EA[ii]*basis(p)
            it.iternext()
    return gradient

##################################################
##################################################
##################################################
def run_until_gradient_converged(sim, fcen, design_region,
                                 forward_fields, basis,
                                 check_interval=20.0, tol=1.0e-4,
                                 max_time=2000):

    adjoint_dfts=sim.add_dft_fields(e_components, fcen, fcen, 1, where=design_region)
    sim.run(until=sim.fields.last_source_time())
    adjoint_fields=[sim.get_dft_array(adjoint_dfts,c,0) for c in e_components]
    last_grad=compute_gradient(sim, design_region, forward_fields, adjoint_fields, basis)
    next_check_time=sim.fields.last_source_time() + check_interval
    done=False
    while not done:
        sim.fields.step()
        if sim.round_time() < next_check_time:
            continue
        next_check_time += check_interval

        adjoint_fields=[sim.get_dft_array(adjoint_dfts,c,0) for c in e_components]
        grad=compute_gradient(sim, design_region, forward_fields, adjoint_fields, basis)
        norm_grad  = np.linalg.norm(grad)
        norm_delta = np.linalg.norm(grad-last_grad)/norm_grad
        last_grad= grad
        if mp.am_master():
           stamp=datetime.datetime.now().strftime("%D::%T: ");
           print stamp, "{:.3e}: |grad|={:.3e}, delta={:.0e}".format(sim.round_time(),norm_grad,norm_delta)
        done = ( (norm_delta<tol) or (sim.round_time() > max_time) )
    
    adjoint_fields = [sim.get_dft_array(adjoint_dfts, c, 0) for c in e_components]
    return grad, adjoint_fields

#########################################################
#########################################################
#########################################################
def place_adjoint_sources(sim, fcen, df, region,
                          objective_fields, mode=0, coeff=1,
                          terms=[1,1,1,1]):
    
    sim.init_fields()
    src_t=mp.GaussianSource(fcen,fwidth=df)
    xyzw=sim.get_dft_array_metadata(center=region.center,size=region.size)
    sim.sources = []
    components = get_flux_components(region.direction)

    peak_time = 5.0*df;
    factor    = np.exp(2.0j*np.pi*fcen*peak_time)/(2.0*df*2.0j*np.pi*fcen)

    for nc in range(0,4):
        if terms[nc]==0:
            continue
        it=np.nditer(objective_fields[1],flags=['multi_index'])
        while not it.finished:
            ii=it.multi_index
            (p,w)=(mp.Vector3(*tuple(xyzw[ii][0:3])), xyzw[ii][3])
            sign=1.0 if (nc%2)==1 else -1.0
            ncbar=(nc+2)%4
            c=components[nc]
            cbar=components[ncbar]
            amp=mode.amplitude(p,cbar) if mode!=0 else objective_fields[ncbar][ii]
            sim.sources += [ mp.Source( src_t, center=p, component=c,
                                        amplitude=factor*w*np.conj(coeff*amp)
                                      )
                           ]
            it.iternext()

#########################################################
# second timestepping run with adjoint sources to compute
# derivatives
#########################################################
def run_adjoint(sim, fcen, df, basis,
                objective_region, forward_objective_fields,
                design_region,    forward_design_fields,
                mode=0, coeff=1, terms=[0,1,0,0]):

    place_adjoint_sources(sim,fcen,df,
                          objective_region, forward_objective_fields,
                          mode=mode, coeff=coeff, terms=terms)
    sim.init_sim()
                          
    gradient,adjoint_design_fields \
     =run_until_gradient_converged(sim, fcen,
                                   design_region, forward_design_fields,
                                   basis,
                                   check_interval=20.0, tol=1.0e-4,
                                   max_time=2000)

    return gradient, adjoint_design_fields

#########################################################
#########################################################
#########################################################
def get_objective_and_gradient(sim, forward_sources,
                               objective_regions, design_region, basis,
                               objective):

    src_time=forward_sources[0].src
    fcen=src.frequency
    df=1/src.width

    sim.restart_fields()
    sim.change_sources(forward_sources)
    sim.init_sim()
    objective_fluxes=[sim.add_flux(fcen,fcen,1,region) for region in objective_regions]
    run_until_flux_converged(sim, objective_fluxes)

    sim.restart_fields()
    objective_fluxes=[sim.add_flux(fcen,fcen,1,region) for region in objective_regions]
    place_adjoint_sources(sim,fcen,df,
                          objective_region, forward_objective_fields,
                          mode=mode, coeff=coeff, terms=terms)
                          
    sim.init_sim()
    gradf,adjoint_design_fields \
     =run_until_gradient_converged(sim, fcen,
                                   design_region, forward_design_fields,
                                   basis,
                                   check_interval=20.0, tol=1.0e-4,
                                   max_time=2000)

    return f, gradf

#########################################################
# the remainder of this file is a set of predefined basis
# functions
#########################################################
def sinusoid(k,u):
    arg = 2.0*np.pi*np.floor((k+1)/2)*u
    return np.sin(arg) if (k%2) else np.cos(arg)

##################################################
# Plane-wave basis for a rectangular region.
#                                                   
# f_{nx,ny} = S_{nx}(x) * S_{ny}(y)
# 
# S_{0}(u)    = 1
# S_{2k-1}(u) = sin(2*pi*k*u)
# S_{2k}  (u) = cos(2*pi*k*u)
#
# The size of the basis is (2*kxMax+1)*(2*kyMax+1)
##################################################
def fourier_basis(kxMax, kyMax, add_offset=False):

    def _get_bf_vector(pbar):
        b=np.zeros( (2*kxMax+1)*(2*kyMax+1) )
        for kx in range(0,2*kxMax+1):
            for ky in range(0,2*kyMax+1):
                b[ kx*(2*kyMax+1) + ky ] = sinusoid(kx,pbar[0])*sinusoid(ky,pbar[1])
        if add_offset:
            for n in range(1,len(b)):
                b[n]+=1.0
        return b
            
    return _get_bf_vector

##################################################
# basis of expansion functions f_{m,n}(rho,phi) for
# a disc or annulus.
#
# f_{nr,np}(rho,phi) = (legendre polynomial in rho) * (sinusoid in phi)
#
# f_{nr, 0 }    = L_nr(u)
# f_{nr, 2*k-1) = L_n(u) * sin(k*phi)
# f_{nr, 2*k)   = L_n(u) * cos(k*phi)
#
# for nr=[0,...,NR], k=[0,1,...,kMax]
#
# Here u is a rescaled version of rho that runs over [-1:1]
# as rho runs over [rho_min:rho_max].
#
# The size of the basis is (NR+1)*(2*kMax+1)
##################################################
def disc_basis(NR=2, kMax=2, rho_max=1.0, rho_min=0.0):

    def _get_bf_vector(pbar):
        b=np.zeros( (NR+1)*(2*kMax+1) )
        (x,y)=pbar[0],pbar[1]
        rho=np.sqrt( x*x + y*y )
        if rho<rho_min or rho>rho_max:
            return b
        phi=np.arctan2(y,x)
        u = -1.0 + 2.0*(rho-rho_min)/(rho_max-rho_min)
        (Lm1, L)=(0.0,1.0)
        for nr in range(0,NR+1):
            for nk in range(0,2*kMax+1):
                b[ nr*(2*kMax+1) + nk ] = L*sinusoid(nk,phi/(2.0*np.pi))
            (Lm2,Lm1)=(Lm1,L)
            L = ( (2*nr+1)*u*Lm1 - nr*Lm2 )/(nr+1)  # Legendre recursion
        return b
            
    return _get_bf_vector

##################################################
# dielectric function for a design region
##################################################
def custom_dielectric(center, size, basis, coeffs, return_medium=True):

 def _eps_value(p):
    # call user's function with scaled/shifted coordinate to get vector of
    # basis-function values, then dot-product with vector of basis-function coefficients
    pbar=[ ((p[n]-center[n]) / (1 if size[n]==0 else size[n])) for n in range(0,3)]
    return np.dot(coeffs, basis(pbar))

 def _eps_medium(p):
    eps=_eps_value(p)
    return mp.Medium( epsilon = 12 if eps<=0 else eps )

 return _eps_medium if return_medium else _eps_value
