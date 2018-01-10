import meep as mp
import numpy as np
import matplotlib.pyplot as plt

##################################################
# x-dependent width of waveguide
##################################################
def h_func(x, L, p, hA, hB):
  x0=x/L
  if (x0 < -0.5):
    return hA;
  if (x0 > +0.5):
    return hB;
  if (p==0):
    return 0.5*(hA+hB) + (hB-hA)*x0;
  else: # if (p==1):
    return 0.5*(hA+hB) + (hB-hA)*x0*(1.5 - 2.0*x0*x0);

##################################################
# user-defined function for spatially-varying epsilon
##################################################
def my_eps_func(loc, L, p, hA, hB, eps_out, eps_in):

    if ( abs(loc.y) > 0.5*h_func(loc.x, L, p, hA, hB) ):
     return eps_out;    # outside waveguide
    else:
     return eps_in;     # inside waveguide

##################################################
##################################################
##################################################
class wvg_taper:

    ##################################################
    # constructor
    ##################################################
    def __init__(self,
                 hA=1.0, hB=3.0,      # smaller, larger waveguide thickness
                 L=2.0, p=0,          # taper length and smoothness index
                 eps_waveguide=11.7,  # permittivity inside waveguide
                 eps_ambient=1.0,     # permittivity of medium
                 LX=5.0, LY=3.0,      # half-lengths of computational cell
                 DPML=0.5,            # PML thickness
                 resolution=25.0,     # grid points per unit length
               ): 

        eps_func = lambda loc: my_eps_func(loc, L, p, hA, hB,
                                                eps_ambient, eps_waveguide)
                                          
        self.sim=mp.Simulation( cell_size=mp.Vector3(2*LX, 2*LY),
                                resolution=resolution,
                                boundary_layers=[mp.PML(DPML)],
                                epsilon_func = eps_func
                              )

        self.sim.init_fields()

    ##################################################
    # plot permittivity over the computational grid
    ##################################################
    def plot_eps(self):

     eps=self.sim.get_array(center    = mp.Vector3(0,0),
                            size      = self.sim.cell_size,
                            component = mp.Dielectric)

     interp='gaussian'
     cmap='coolwarm'
     LX=0.5*self.sim.cell_size.x;
     LY=0.5*self.sim.cell_size.y;
     xtnt=[-LX,LX,-LY,LY];
     plt.figure()
     plt.imshow(eps.transpose(), interpolation=interp, cmap=cmap, extent=xtnt)
     plt.xlabel("x")
     plt.ylabel("y")
     plt.colorbar()
     plt.show(block=False)

    ##################################################
    # plot eigenmode profiles
    ##################################################
    def plot_modes(self):
       
       res=1.0*self.sim.resolution;
       LX=0.5*self.sim.cell_size.x;
       LY=0.5*self.sim.cell_size.y;
       xA=-0.5*LX;
       xB=+0.5*LX;
       vA=mp.volume( mp.vec(xA, -LY), mp.vec(xA,+LY) )
       kpoint=mp.vec(0.426302,0);
       omega=0.15;
       band_num=1;
       parity=0;
       match_frequency=True;
       tol=1.0e-4;
# FIXME why does python complain about this function call
# having incorrect arguments?
       vedata = self.sim.fields.get_eigenmode(omega, mp.X, vA, vA,
                                              band_num, kpoint,
                                              match_frequency,
                                              parity, res, tol);
       vgrp = get_group_velocity(vedata);
       print("vgrp: {}".format(vgrp))

    ##################################################
    # add an eigenmode-source excitation for the #band_numth mode
    # of the smaller waveguide, then timestep to accumulate DFT
    # flux in the larger waveguide.
    # if frame_interval>0, a movie is created showing
    # the fields on the xy plane with one frame
    # every frame_interval time units (in meep time)
    ##################################################
    def get_flux(self, fcen=0.15, df=0.075, nfreq=1, band_num=1,
                 frame_interval=0):
       
       #--------------------------------------------------
       # add eigenmode source at midpoint of smaller waveguide
       #--------------------------------------------------
       f=self.sim.fields;
       res=1.0*self.sim.resolution;
       LX=0.5*self.sim.cell_size.x;
       LY=0.5*self.sim.cell_size.y;
       xA=-0.5*LX;
       xB=+0.5*LX;
       vA=mp.volume( mp.vec(xA, -LY), mp.vec(xA,+LY) )
       vB=mp.volume( mp.vec(xB, -LY), mp.vec(xB,+LY) )
       vC=mp.volume( mp.vec(-LX, -LY), mp.vec(LX,LY) )
       src=mp.GaussianSource(fcen, fwidth=df);
       kpoint=mp.vec(0.426302,0);
       parity=0;
       match_frequency=True;
       tol=1.0e-4;
       amp=1.0;
       f.add_eigenmode_source(mp.Dielectric, src, mp.X, vA, vA, 
                              band_num, kpoint, match_frequency,
                              parity, res, tol, amp);

       #--------------------------------------------------
       # add DFT flux region at midpoint of larger waveguide
       #--------------------------------------------------
       fluxB=f.add_dft_flux_plane(vB, fcen-0.5*df, fcen+0.5*df, nfreq);

       #--------------------------------------------------
       # for DFT flux region for moviemaking if requested
       #--------------------------------------------------
       fluxC=0
       if frame_interval>0:
         fluxC=f.add_dft_flux_plane(vC, fcen-0.5*df, fcen+0.5*df, nfreq);

       #--------------------------------------------------
       # timestep until Poynting flux through larger waveguide has 
       # decayed to 0.1% its max value
       #--------------------------------------------------
       pvInterval=1.0; # check PV decay every 1.0 meep time units
       nextPVTime=f.round_time() + pvInterval;
       nextFrameTime=f.round_time();
       MaxPV=0.0;
       Stop=False;
       while Stop==False:

         f.step();

         # check for poynting-flux decay at regular intervals
         FieldsDecayed=False;
         if f.round_time() > nextPVTime:
             nextPVTime += pvInterval;
             ThisPV=f.flux_in_box(mp.X,vB)
             if (ThisPV > MaxPV):
                MaxPV = ThisPV;
             elif (ThisPV < 0.001*MaxPV):
                FieldsDecayed=True;

         # output movie frames at regular intervals if requested
         # TODO implement me

         SourcesFinished = f.round_time() > f.last_source_time();
         Stop = (SourcesFinished and FieldsDecayed);
         
       print("finished timestepping at {}".format(f.round_time()))
       return fluxB

    ##################################################
    # postprocessing: write DFT fields in larger waveguide to HDF5
    ##################################################
    def flux2hdf5(self, flux, vol, filename):
      self.sim.fields.output_flux_fields(flux, vol, filename);

    ##################################################
    # postprocessing: compute coefficients in normal-mode
    # expansion of DFT fields in larger waveguide
    ##################################################
    def get_eigenmode_coefficients(self, flux, d, vol,
                                   bands, k_guess, k_guess_data):
      f=self.sim.fields
      return f.get_eigenmode_coefficients(flux, d, vol, bands,
                                          k_guess, k_guess_data)

##################################################
##################################################
##################################################
wt=wvg_taper();
wt.plot_eps();
#wt.plot_modes();
fluxB = wt.get_flux();

LX=0.5*wt.sim.cell_size.x;
LY=0.5*wt.sim.cell_size.y;
xB=+0.5*LX;
vB=mp.volume( mp.vec(xB, -LY), mp.vec(xB,+LY) )
wt.flux2hdf5(fluxB, vB, "vBFields")
