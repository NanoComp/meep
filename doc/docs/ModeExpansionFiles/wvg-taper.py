import meep as mp
import numpy as np
import matplotlib.pyplot as plt

##################################################
# x-dependent width of waveguide
##################################################
def w_func(x, L, p, wA, wB):
  x0=x/L
  if (x0 < -0.5):
    return wA;
  if (x0 > +0.5):
    return wB;
  if (p==0):
    return 0.5*(wA+wB) + (wB-wA)*x0;
  else: # if (p==1):
    return 0.5*(wA+wB) + (wB-wA)*x0*(1.5 - 2.0*x0*x0);

##################################################
# user-defined function for spatially-varying epsilon
##################################################
def my_eps_func(loc, L, p, wA, wB, eps_out, eps_in):

    if ( abs(loc.y) > 0.5*w_func(loc.x, L, p, wA, wB) ):
     return eps_out;    # outside waveguide
    else:
     return eps_in;     # inside waveguide

##################################################
# user-provided functions to estimate k-vectors for waveguide eigenmodes
##################################################
def equal_float(a,b,Tol=1.0e-6):
    if ( abs(a-b) < Tol*max(abs(a),abs(b)) ):
        return True;
    else:
        return False;

def k_guess(freq, band_num, w):
  
    # hard-coded dispersion relations for waveguides of given size
    if ( equal_float(w,1.0) and equal_float(freq, 0.15) ):
        if (band_num>=1): return mp.vec(0.419984,0)

    if ( equal_float(w,3.0) and equal_float(freq, 0.15) ):
        if (band_num==1):
            return mp.vec(0.494476,0)
        if (band_num==2):
            return mp.vec(0.486399,0)
        if (band_num==3):
            return mp.vec(0.435861,0)
        if (band_num==4):
            return mp.vec(0.397068,0)
        if (band_num==5):
            return mp.vec(0.322812,0)
        if (band_num>=6):
            return mp.vec(0.211186,0)

    return mp.vec(0.419987,0)

##################################################
##################################################
##################################################
class wvg_taper:

    ##################################################
    # constructor
    ##################################################
    def __init__(self,
                 wA=1.0, wB=3.0,      # smaller, larger waveguide thickness
                 L=2.0, p=0,          # taper length and smoothness index
                 eps_waveguide=11.7,  # permittivity inside waveguide
                 eps_ambient=1.0,     # permittivity of medium
                 LX=5.0, LY=3.0,      # half-lengths of computational cell
                 DPML=0.5,            # PML thickness
                 resolution=25.0,     # grid points per unit length
               ): 

        self.force_complex_fields = True;

        #--------------------------------------------------------------------
        #- initialize structure and fields with user-defined epsilon function
        #--------------------------------------------------------------------
   
        eps_func = lambda loc: my_eps_func(loc, L, p, wA, wB,
                                           eps_ambient, eps_waveguide)
                                          
        self.sim=mp.Simulation( cell_size=mp.Vector3(2*LX, 2*LY),
                                resolution=resolution,
                                boundary_layers=[mp.PML(DPML)],
                                force_complex_fields=True,
                                epsilon_func = eps_func
                              )

        self.sim.run(mp.at_beginning(mp.output_epsilon), until=1.0);
        f=self.sim.fields;

        #--------------------------------------------------
        # add DFT flux regions at midpoints of smaller and larger waveguides
        #--------------------------------------------------
        xA=-0.5*LX;
        xB=+0.5*LX;
        #LYP=LY-DPML;
        LYP=LY;
        self.wA=wA;
        self.wB=wB;
        self.vA=mp.volume( mp.vec(xA, -LYP), mp.vec(xA,+LYP) )
        self.vB=mp.volume( mp.vec(xB, -LYP), mp.vec(xB,+LYP) )
        fcen=0.15;
        df=0.5*fcen;
        nf=1;
        print("Adding flux region A...");
        self.fluxA=f.add_dft_flux_plane(self.vA, fcen-0.5*df, fcen+0.5*df, nf);
        print("Adding flux region B...");
        self.fluxB=f.add_dft_flux_plane(self.vB, fcen-0.5*df, fcen+0.5*df, nf);

    ##################################################
    # plot permittivity over the computational grid
    ##################################################
    def plot_eps(self):

     eps=self.sim.get_array(center    = mp.Vector3(0,0),
                            size      = self.sim.cell_size,
                            component = mp.Dielectric)
     print("eps shape={}".format(eps.shape))

     interp='gaussian'
     cmap='coolwarm'
     LX=0.5*self.sim.cell_size.x;
     LY=0.5*self.sim.cell_size.y;
     extent=[-LX,LX,-LY,LY];
     plt.figure()
     plt.imshow(eps.transpose(), interpolation=interp, cmap=cmap, extent=extent)
     plt.xlabel("x")
     plt.ylabel("y")
     plt.colorbar()
     plt.show(block=False)

    ##################################################
    # generate plots of eigenmode profiles
    ##################################################
    def plot_modes(self):
       
       f=self.sim.fields;
       vA=self.vA;
       vB=self.vB;
       wA=self.wA;
       wB=self.wB;
       fluxA=self.fluxA;
       fluxB=self.fluxB;
       res=1.0*self.sim.resolution;
       freq=0.15;
       parity=0;
       match_freq=True;
       tol=1.0e-4;

       # eigenmode #1 in narrower waveguide
       modeA1 = f.get_eigenmode(freq, mp.X, vA, vA,
                                1, k_guess(freq, 1, wA),
                                match_freq, parity, res, tol);
       f.output_mode_fields(modeA1, fluxA, vA, "modeA1");
       if mp.am_master():
           print("vgrp(A1)={}".format(mp.get_group_velocity(modeA1)));

       # eigenmodes #1-4 in wider waveguide
       modeB1 = f.get_eigenmode(freq, mp.X, vB, vB,
                                1, k_guess(freq,1,wB),
                                match_freq, parity, res, tol);
       f.output_mode_fields(modeB1, fluxB, vB, "modeB1");
       if mp.am_master():
           print("vgrp(B1)={}".format(mp.get_group_velocity(modeB1)));

       modeB2 = f.get_eigenmode(freq, mp.X, vB, vB,
                                2, k_guess(freq,2,wB),
                                match_freq, parity, res, tol);
       f.output_mode_fields(modeB2, fluxB, vB, "modeB2");
       if mp.am_master():
           print("vgrp(B2)={}".format(mp.get_group_velocity(modeB1)));

       modeB3 = f.get_eigenmode(freq, mp.X, vB, vB,
                                3, k_guess(freq,3,wB),
                                match_freq, parity, res, tol);
       f.output_mode_fields(modeB3, fluxB, vB, "modeB3");
       if mp.am_master():
           print("vgrp(B3)={}".format(mp.get_group_velocity(modeB3)));

       modeB4 = f.get_eigenmode(freq, mp.X, vB, vB,
                                4, k_guess(freq,4,wB),
                                match_freq, parity, res, tol);
       f.output_mode_fields(modeB4, fluxB, vB, "modeB4");
       if mp.am_master():
           print("vgrp(B4)={}".format(mp.get_group_velocity(modeB4)));

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
       vA=self.vA;
       vB=self.vB;
       res=1.0*self.sim.resolution;
       src=mp.GaussianSource(fcen, fwidth=df);
       kpoint=k_guess(fcen,band_num,self.wA)
       parity=0;
       match_freq=True;
       tol=1.0e-4;
       amp=1.0;
       f.add_eigenmode_source(mp.Dielectric, src, mp.X, vA, vA, band_num,
                              kpoint, match_freq, parity, res, tol, amp);

       #--------------------------------------------------
       # add DFT flux region for moviemaking if requested
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
##################################################
         frameInterval=5.0;
         if f.round_time() > nextFrameTime:
             nextFrameTime += frameInterval;
             eps=self.sim.get_array(center    = mp.Vector3(0,0),
                                    size      = self.sim.cell_size,
                                    component = mp.Sx)
             print("Sx shape={}".format(eps.shape))
             plt.clf()
             plt.imshow(eps.transpose(), interpolation='gaussian', cmap='coolwarm')
             plt.show()
##################################################

         SourcesFinished = f.round_time() > f.last_source_time();
         Stop = (SourcesFinished and FieldsDecayed);
         
       print("finished timestepping at {}".format(f.round_time()))
       f.output_flux_fields(wt.fluxA, wt.vA, "fluxA");
       f.output_flux_fields(wt.fluxB, wt.vB, "fluxB");

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
#
##################################################
##################################################
##################################################
wt=wvg_taper();
wt.plot_eps();

wt.get_flux();
#wt.plot_modes();
