import meep as mp
import numpy as np
import h5py as h5py
import matplotlib.pyplot as plt
import matplotlib.animation as animation

##################################################
# x-dependent width of waveguide
##################################################
def w_func(x, L, p, wA, wB):
  if L==0:
    return wA if x<0 else wB
  x0=x/L
  if (x0 < -0.5):
    return wA;
  elif (x0 > +0.5):
    return wB;
  elif p==2:
    return 0.5*(wA+wB) + (wB-wA)*x0*(15.0 + x0*x0*(-40.0 + x0*x0*48.0))/8.0;
  elif p==1:
    return 0.5*(wA+wB) + (wB-wA)*x0*(1.5 - 2.0*x0*x0);
  else: # default t p==0, simple linear taper
    return 0.5*(wA+wB) + (wB-wA)*x0;

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
            return mp.vec(0.494499,0)
        if (band_num==2):
            return mp.vec(0.486657,0)
        if (band_num==3):
            return mp.vec(0.434539,0)
        if (band_num==4):
            return mp.vec(0.397068,0)
        if (band_num==5):
            return mp.vec(0.322812,0)
        if (band_num>=6):
            return mp.vec(0.211186,0)

    return mp.vec(0.0,0.0)

##################################################
##################################################
##################################################
class wvg_taper:

    ##################################################
    # constructor
    ##################################################
    def __init__(self,
                 wA=1.0, wB=3.0,        # smaller, larger waveguide thickness
                 LWaveguide=3.0,        # length of each waveguide section
                 LTaper=3.0, pTaper=0,  # taper length and smoothness index
                 eps_waveguide=11.7,    # permittivity inside waveguide
                 eps_ambient=1.0,       # permittivity of medium
                 LY=6.0,                # width of computational cell
                 DPML=0.5,              # PML thickness
                 fcen=0.15, df=0.075,   # center frequency / width
                 band_num=1,            # index of eigenmode source
                 resolution=25.0,       # grid points per unit length
               ): 

        #--------------------------------------------------------------------
        #- user-defined epsilon function
        #--------------------------------------------------------------------
        eps_func = lambda loc: my_eps_func(loc, LTaper, pTaper, wA, wB,
                                           eps_ambient, eps_waveguide)

        #--------------------------------------------------------------------
        #- eigenmode source at midpoint of smaller waveguide
        #--------------------------------------------------------------------
        LX = 2.0*(DPML + LWaveguide) + LTaper;
        xA = -0.5*LX + DPML + 0.5*LWaveguide;
        xB = +0.5*LX - DPML - 0.5*LWaveguide;
        sources = [ mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=df),
                                       center=mp.Vector3(xA,0.0),
                                       size=mp.Vector3(0.0,LY),
                                       eig_band=band_num
                                      )
                  ]
                                          
        self.sim=mp.Simulation( cell_size=mp.Vector3(LX, LY),
                                resolution=resolution,
                                boundary_layers=[mp.PML(DPML)],
                                force_complex_fields=True,
                                epsilon_func = eps_func,
                                sources=sources
                              )
       
        self.sim.run(mp.at_beginning(mp.output_epsilon), until=1.0)
        f=self.sim.fields;

        #--------------------------------------------------
        # add DFT flux regions at midpoints of smaller and larger waveguides
        #--------------------------------------------------
        YP            = 0.5*LY - DPML;
        self.vA       = mp.volume( mp.vec(xA, -YP), mp.vec(xA,+YP) )
        self.vB       = mp.volume( mp.vec(xB, -YP), mp.vec(xB,+YP) )
        nf=1;
        self.fluxA=f.add_dft_flux_plane(self.vA, fcen-0.5*df, fcen+0.5*df, nf);
        self.fluxB=f.add_dft_flux_plane(self.vB, fcen-0.5*df, fcen+0.5*df, nf);

        #--------------------------------------------------
        # save some other fields in the wvg_taper class for later use
        #--------------------------------------------------
        self.xA       = xA;
        self.xB       = xB;
        self.wA       = wA;
        self.wB       = wB;
        self.LTaper   = LTaper;
        self.pTaper   = pTaper;
        self.fcen     = fcen;
        self.df       = df;
        self.band_num = band_num;

    ##################################################
    # plot permittivity over the computational grid
    ##################################################
    def plot_eps(self):

     eps=self.sim.get_array(center    = mp.Vector3(0,0),
                            size      = self.sim.cell_size,
                            component = mp.Dielectric)

     interp='gaussian'
     cmap='coolwarm'
     LX=self.sim.cell_size.x;
     LY=self.sim.cell_size.y;
     extent=[-0.5*LX,0.5*LX,-0.5*LY,0.5*LY];
     plt.figure()
 
     # plot epsilon grid
     plt.rc('xtick', labelsize=15) 
     plt.rc('ytick', labelsize=15) 
     plt.imshow(eps.transpose(), interpolation=interp, cmap=cmap, extent=extent)
     plt.xlabel(r"$x$",fontsize=40)
     plt.ylabel(r"$y$",fontsize=40,rotation=0)
     plt.colorbar()

     # add lines and annotations to indicate locations of flux planes
     plt.plot(self.xA*np.ones(100), np.linspace(-0.5*LY, 0.5*LY, 100), "--", lw=5, color='white')
     plt.text(self.xA, -0.5*LY-0.5, r'$x_A$', fontsize=40, horizontalalignment='center')
     plt.plot(self.xB*np.ones(100), np.linspace(-0.5*LY, 0.5*LY, 100), "--", lw=5, color='white');
     plt.text(self.xB, -0.5*LY-0.5, r'$x_B$', fontsize=40, horizontalalignment='center')

     plt.show(block=False)

    ##################################################
    # subroutine invoked by plot_modes and plot_flux
    # to plot field components stored in an hdf5 file
    ##################################################
    def plot_fields(self, num_rows, plot_row, file_name, row_name, is_flux):

      # read (tranverse) field components on cross-sectional plane from file
      file = h5py.File(file_name + ".h5", 'r')
      suffix = "_0" if is_flux else "";
      ey = file["ey" + suffix + ".r"][()] + 1j*file["ey" + suffix + ".i"][()];
      ez = file["ez" + suffix + ".r"][()] + 1j*file["ez" + suffix + ".i"][()];
      hy = file["hy" + suffix + ".r"][()] + 1j*file["hy" + suffix + ".i"][()];
      hz = file["hz" + suffix + ".r"][()] + 1j*file["ez" + suffix + ".i"][()];
      if np.size(ey.shape)>1:
        ey = np.real(ey[1,:]);
        ez = np.real(ez[1,:]);
        hy = np.real(hy[1,:]);
        hz = np.real(hz[1,:]);

      # plot ey,ez,hy,hz
      eMax=max(np.absolute(ey).max(), np.absolute(ez).max());
      hMax=max(np.absolute(hy).max(), np.absolute(hz).max());
      NY=ey.shape[0];
      yMax=self.vA.get_min_corner().y();
      yRange=np.linspace(-yMax,yMax,NY);
      r0=4*(plot_row-1);

      # plot labels
      plt.subplot(num_rows,4,r0+1); plt.plot(yRange,ey); plt.ylim(-eMax,eMax);
      plt.subplot(num_rows,4,r0+2); plt.plot(yRange,ez); plt.ylim(-eMax,eMax);
      plt.subplot(num_rows,4,r0+3); plt.plot(yRange,hy); plt.ylim(-hMax,hMax);
      plt.subplot(num_rows,4,r0+4); plt.plot(yRange,hz); plt.ylim(-hMax,hMax);
      plt.subplot(num_rows,4,1);    plt.title(r're($E_y$)',fontsize=40)
      plt.subplot(num_rows,4,2);    plt.title(r're($E_z$)',fontsize=40)
      plt.subplot(num_rows,4,3);    plt.title(r're($H_y$)',fontsize=40)
      plt.subplot(num_rows,4,4);    plt.title(r're($H_z$)',fontsize=40)
      plt.subplot(num_rows,4,r0+1); plt.ylabel(row_name,fontsize=32, rotation=0,
                                               horizontalalignment='right')
      plt.show(block=False);
      
    ##################################################
    # generate plots of eigenmode profiles
    ##################################################
    def plot_modes(self, nbMax=7):
       
       ##################################################
       # calculate the eigenmodes and write field components
       # on cross-sectional planes to HDF5 files
       ##################################################
       f     = self.sim.fields;
       vA    = self.vA;
       vB    = self.vB;
       fcen  = self.fcen;

       plt.clf();

       # eigenmode #1 in narrower waveguide
       nb=self.band_num;
       k0=k_guess(fcen,nb,self.wA);
       modeA = f.get_eigenmode(fcen, mp.X, vA, vA, nb, k0, True, 0, 0.0, 1.0e-4)
       f.output_mode_fields(modeA, self.fluxA, vA, "modeA");
       self.plot_fields(nbMax, 1, "modeA", "Mode A1", False);

       # eigenmodes #1-6 in wider waveguide
       for nb in range(1,nbMax):
         k0=k_guess(fcen,nb,self.wB);
         modeB = f.get_eigenmode(fcen, mp.X, vB, vB, nb, k0, True, 0, 0.0, 1.0e-4)
         file_name="modeB" + str(nb);
         row_name="Mode B" + str(nb);
         f.output_mode_fields(modeB, self.fluxB, vB, file_name);
         self.plot_fields(nbMax, nb+1, file_name, row_name, False);

       for np in range(25,29):
           plt.subplot(nbMax,4,np); plt.xlabel(r'$y$',fontsize=32)

    ##################################################
    # add an eigenmode-source excitation for the #band_numth mode
    # of the smaller waveguide, then timestep to accumulate DFT
    # flux in the larger waveguide.
    # if frame_interval>0, a movie is created showing
    # the fields on the xy plane with one frame
    # every frame_interval time units (in meep time)
    ##################################################
    def get_flux(self, frame_interval=0):

       #--------------------------------------------------
       #--------------------------------------------------
       #--------------------------------------------------
       nextFrameTime = 1e100;
       f=self.sim.fields;
       vA=self.vA;
       vB=self.vB;
       LTaper=self.LTaper;
       pTaper=self.pTaper;
       origin=mp.Vector3(0,0,0)
       LX=self.sim.cell_size.x;
       LY=self.sim.cell_size.y;
       extent=[-0.5*LX,0.5*LX,-0.5*LY,0.5*LY];
       full_grid=mp.Vector3(LX, LY, 0)
       frames=[]
       fig=0
       plt.ion()
       num_frames=0
       if frame_interval>0:
         fig=plt.figure();
         nextFrameTime=f.round_time();

       #--------------------------------------------------
       # timestep until Poynting flux through larger waveguide has 
       # decayed to 0.1% its max value
       #--------------------------------------------------
       pvInterval=1.0; # check PV decay every 1.0 meep time units
       nextPVTime=f.round_time() + pvInterval;
       MaxPV=0.0;
       Stop=False;
       SMax=0.0; # max poynting flux at any point
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

         # write movie-frame files at regular intervals if requested
         if f.round_time() > nextFrameTime:
             nextFrameTime += frame_interval;
             eps=self.sim.get_array(center=origin,size=full_grid,component=mp.Dielectric)
             Sx=self.sim.get_array(center=origin,size=full_grid,component=mp.Sx)
             SMax = max( abs(Sx.max()), SMax )
             plt.clf();
             eps=-eps*SMax/11.7;
             plt.imshow(eps.transpose(), extent=extent,
                        interpolation='gaussian', cmap='coolwarm')
             plt.imshow(Sx.transpose(), alpha=0.5, extent=extent,
                        interpolation='gaussian', cmap='coolwarm')
             plt.savefig('L{}_p{}_f{}.png'.format(LTaper,pTaper,num_frames));
             num_frames+=1;

         SourcesFinished = f.round_time() > f.last_source_time();
         Stop = (SourcesFinished and FieldsDecayed);
         
       print("finished timestepping at {}".format(f.round_time()))
       file_name="fluxB_p" + str(self.pTaper) + "_L" + str(self.LTaper);
       f.output_flux_fields(wt.fluxB, wt.vB, file_name);
       if self.LTaper==0:
           row_name="Flux B\n(L=0)"
       else:
           row_name="Flux B\n(L=" + str(self.LTaper) + ", p=" + str(self.pTaper) + ")"
       nbMax=7;
       self.plot_fields(nbMax,1,file_name,row_name,True)
       plt.subplot(nbMax,4,1);
       for np in range(1,5):
           plt.subplot(nbMax,4,np); plt.xlabel(r'$y$',fontsize=32)
       
       # postprocessing to generate movie:
       # ffmpeg -i frame%d.png output.mpg

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
wt=wvg_taper(LTaper=0.0,pTaper=0,resolution=25)
#wt.plot_modes()
wt.get_flux(frame_interval=1);
wt=wvg_taper(LTaper=1.0,pTaper=0,resolution=50)
wt.get_flux(frame_interval=1);
wt=wvg_taper(LTaper=2.0,pTaper=0,resolution=50)
wt.get_flux(frame_interval=1);
wt=wvg_taper(LTaper=3.0,pTaper=0,resolution=50)
wt.get_flux(frame_interval=1);
wt=wvg_taper(LTaper=3.0,pTaper=0,resolution=50)
wt.get_flux(frame_interval=1);
wt=wvg_taper(LTaper=3.0,pTaper=1,resolution=50)
wt.get_flux(frame_interval=1);
wt=wvg_taper(LTaper=3.0,pTaper=2,resolution=50)
wt.get_flux(frame_interval=1);
