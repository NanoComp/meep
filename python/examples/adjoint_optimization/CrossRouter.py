import numpy as np
import meep as mp

from meep.adjoint import (OptimizationProblem, FluxLine,
                          xHat, yHat, zHat, origin,
                          ParameterizedDielectric, FiniteElementBasis)

##################################################
##################################################
##################################################
class CrossRouter(OptimizationProblem):

    ##################################################
    ##################################################
    ##################################################
    def add_args(self, parser):

        # add new problem-specific arguments
        parser.add_argument('--wh',       type=float, default=1.5,  help='width of horizontal waveguide')
        parser.add_argument('--wv',       type=float, default=1.5,  help='width of vertical waveguide')
        parser.add_argument('--l_stub',   type=float, default=3.0,  help='waveguide input/output stub length')
        parser.add_argument('--eps',      type=float, default=6.0,  help='waveguide permittivity')
        parser.add_argument('--r_design', type=float, default=0.0,  help='design region radius')
        parser.add_argument('--l_design', type=float, default=4.0,  help='design region side length')
        parser.add_argument('--nfe',      type=int,   default=2,    help='number of finite elements per unit length')

        # the following three quantities are weighting prefactors assigned to
        # the north, south, and east power fluxes in the objective function
        parser.add_argument('--n_weight', type=float, default=1.00, help='')
        parser.add_argument('--s_weight', type=float, default=0.00, help='')
        parser.add_argument('--e_weight', type=float, default=0.00, help='')

        # set problem-specific defaults for existing (general) arguments
        parser.set_defaults(fcen=0.5)
        parser.set_defaults(df=0.2)
        parser.set_defaults(dpml=1.0)
        parser.set_defaults(epsilon_design=6.0)

    ##################################################
    ##################################################
    ##################################################
    def init_problem(self, args):

        #----------------------------------------
        # size of computational cell
        #----------------------------------------
        lcen          = 1.0/args.fcen
        dpml          = 0.5*lcen if args.dpml == -1.0 else args.dpml
        design_length = 2.0*args.r_design if args.r_design > 0.0 else args.l_design
        sx = sy       = dpml + args.l_stub + design_length + args.l_stub + dpml
        cell_size     = mp.Vector3(sx, sy, 0.0)

        #----------------------------------------
        #- design region bounding box
        #----------------------------------------
        design_center = origin
        design_size   = mp.Vector3(design_length, design_length)
        design_region = mp.Volume(center=design_center, size=design_size)

        #----------------------------------------
        #- objective and source regions
        #----------------------------------------
        gap            =  args.l_stub/6.0                    # gap between source region and flux monitor
        d_flux         =  0.5*(design_length + args.l_stub)  # distance from origin to NSEW flux monitors
        d_source       =  d_flux + gap                       # distance from origin to source
        d_flx2         =  d_flux + 2.0*gap
        l_flux_NS      =  2.0*args.wv
        l_flux_EW      =  2.0*args.wh
        north          =  FluxLine(0.0, +d_flux, l_flux_NS, mp.Y, 'north')
        south          =  FluxLine(0.0, -d_flux, l_flux_NS, mp.Y, 'south')
        east           =  FluxLine(+d_flux, 0.0, l_flux_EW, mp.X, 'east')
        west1          =  FluxLine(-d_flux, 0.0, l_flux_EW, mp.X, 'west1')
        west2          =  FluxLine(-d_flx2, 0.0, l_flux_EW, mp.X, 'west2')

        objective_regions  = [north, south, east, west1, west2]

        source_center  =  mp.Vector3(-d_source, 0.0)
        source_size    =  mp.Vector3(0.0,l_flux_EW)

        #----------------------------------------
        #- optional extra regions for visualization
        #----------------------------------------
        extra_regions  = [mp.Volume(center=origin, size=cell_size)] if args.full_dfts else []

        #----------------------------------------
        # basis set
        #----------------------------------------
        basis = FiniteElementBasis(lx=args.l_design, ly=args.l_design, density=args.nfe)

        #----------------------------------------
        #- objective function
        #----------------------------------------
        fstr=(   '   {:s}*Abs(P1_north)**2'.format('0.0' if args.n_weight==0.0 else '{}'.format(args.n_weight))
               + ' + {:s}*Abs(M1_south)**2'.format('0.0' if args.s_weight==0.0 else '{}'.format(args.s_weight))
               + ' + {:s}*Abs(P1_east)**2'.format('0.0'  if args.e_weight==0.0 else '{}'.format(args.e_weight))
               + ' + 0.0*(P1_north + M1_south + P1_east + P1_west1 + P1_west2)'
               + ' + 0.0*(M1_north + M1_south + M1_east + M1_west1 + M1_west2)'
               + ' + 0.0*(S_north + S_south + S_east + S_west1 + S_west2)'
             )

        #----------------------------------------
        #- internal storage for variables needed later
        #----------------------------------------
        self.args            = args
        self.dpml            = dpml
        self.cell_size       = cell_size
        self.basis           = basis
        self.design_center   = design_center
        self.design_size     = design_size
        self.source_center   = source_center
        self.source_size     = source_size

        if args.eps_design is None:
            args.eps_design = args.eps

        return fstr, objective_regions, extra_regions, design_region, basis

    ##############################################################
    ##############################################################
    ##############################################################
    def create_sim(self, beta_vector, vacuum=False):

        args=self.args

        hwvg=mp.Block(center=origin, material=mp.Medium(epsilon=args.eps),
                      size=mp.Vector3(self.cell_size.x,args.wh))
        vwvg=mp.Block(center=origin, material=mp.Medium(epsilon=args.eps),
                      size=mp.Vector3(args.wv,self.cell_size.y))

        if args.r_design>0.0:
            router=mp.Cylinder(center=self.design_center, radius=args.r_design,
                               epsilon_func=ParameterizedDielectric(self.design_center,
                                                                    self.basis,
                                                                    beta_vector))
        else:
            router=mp.Block(center=self.design_center, size=self.design_size,
                            epsilon_func=ParameterizedDielectric(self.design_center,
                                                                 self.basis,
                                                                 beta_vector))
        geometry=[hwvg, vwvg, router]

        envelope = mp.GaussianSource(args.fcen,fwidth=args.df)
        amp=1.0
        if callable(getattr(envelope, "fourier_transform", None)):
            amp /= envelope.fourier_transform(args.fcen)
        sources=[mp.EigenModeSource(src=envelope,
                                    center=self.source_center,
                                    size=self.source_size,
                                    eig_band=args.source_mode,
                                    amplitude=amp
                                   )
                ]

        sim=mp.Simulation(resolution=args.res, cell_size=self.cell_size,
                          boundary_layers=[mp.PML(self.dpml)], geometry=geometry,
                          sources=sources)

        if args.complex_fields:
            sim.force_complex_fields=True

        return sim

######################################################################
# if executed as a script, we look at our own filename to figure out
# the name of the class above, create an instance of this class called
# op, and call its run() method.
######################################################################
if __name__ == '__main__':
    op=globals()[__file__.split('/')[-1].split('.')[0]]()
    op.run()
