import sys
import argparse
import numpy as np
import meep as mp

from meep.adjoint import (OptimizationProblem, DFTCell, adjoint_options,
                          xHat, yHat, zHat, origin, FluxLine,
                          ParameterizedDielectric, FourierLegendreBasis)

##################################################
##################################################
##################################################
class HoleCloak(OptimizationProblem):

    ##################################################
    ##################################################
    ##################################################
    def add_args(self, parser):

        # add new problem-specific arguments
        parser.add_argument('--dair',        type=float, default=-1.0, help='')
        parser.add_argument('--w_wvg',       type=float, default=4.0,  help='')
        parser.add_argument('--eps_wvg',     type=float, default=6.0,  help='')
        parser.add_argument('--r_disc',      type=float, default=0.5,  help='')
        parser.add_argument('--r_cloak',     type=float, default=1.5,  help='')
        parser.add_argument('--nr_max',      type=int,   default=3,    help='')
        parser.add_argument('--kphi_max',    type=int,   default=2,    help='')
        parser.add_argument('--eps_disc',    type=float, default=1.0,  help='permittivity in hole region (0.0 for PEC)')

        # set problem-specific defaults for existing (general) arguments
        parser.set_defaults(fcen=0.5)
        parser.set_defaults(df=0.2)
        parser.set_defaults(dpml=1.0)

    ##################################################
    ##################################################
    ##################################################
    def init_problem(self, args):

        #----------------------------------------
        # size of computational cell
        #----------------------------------------
        lcen       = 1.0/args.fcen
        dpml       = 0.5*lcen if args.dpml==-1.0 else args.dpml
        dair       = 0.5*args.w_wvg if args.dair==-1.0 else args.dair
        L          = 3.0*lcen
        Lmin       = 6.0*dpml + 2.0*args.r_cloak
        L          = max(L,Lmin)
        sx         = dpml+L+dpml
        sy         = dpml+dair+args.w_wvg+dair+dpml
        cell_size  = mp.Vector3(sx, sy, 0.0)

        #----------------------------------------
        #- design region
        #----------------------------------------
        design_center = origin
        design_size   = mp.Vector3(2.0*args.r_cloak, 2.0*args.r_cloak)
        design_region = mp.Volume(center=design_center, size=design_size)

        #----------------------------------------
        #- objective regions
        #----------------------------------------
        fluxW_center  =  (+args.r_cloak+ dpml)*xHat
        fluxE_center  =  (-args.r_cloak- dpml)*xHat
        flux_size     =  2.0*args.w_wvg*yHat

        #fluxW_region  = mp.FluxRegion(center=fluxW_center, size=flux_size, direction=mp.X)
        #fluxE_region  = mp.FluxRegion(center=fluxE_center, size=flux_size, direction=mp.X)
        x0_east       =  args.r_cloak + dpml
        x0_west       = -args.r_cloak - dpml
        y0            = 0.0
        flux_length   = 2.0*args.w_wvg
        east          = FluxLine(x0_east,y0,flux_length,mp.X,'east')
        west          = FluxLine(x0_west,y0,flux_length,mp.X,'west')

        objective_regions  = [east, west]

        #----------------------------------------
        #- optional extra regions for visualization
        #----------------------------------------
        extra_regions      = [mp.Volume(center=origin, size=cell_size)] if args.full_dfts else []

        #----------------------------------------
        # basis set
        #----------------------------------------
        basis = FourierLegendreBasis(outer_radius=args.r_cloak, inner_radius=args.r_disc,
                                     nr_max=args.nr_max, kphi_max=args.kphi_max)

        #----------------------------------------
        #- source location
        #----------------------------------------
        source_center    = (x0_west-dpml)*xHat
        source_size      = flux_length*yHat

        #----------------------------------------
        #- objective function
        #----------------------------------------
        fstr='Abs(P1_east)**2+0.0*(P1_west + M1_east + M1_west + S_west + S_east)'

        #----------------------------------------
        #- internal storage for variables needed later
        #----------------------------------------
        self.args            = args
        self.dpml            = dpml
        self.cell_size       = cell_size
        self.basis           = basis
        self.design_center   = design_center
        self.source_center   = source_center
        self.source_size     = source_size

        return fstr, objective_regions, extra_regions, design_region, basis

    ##############################################################
    ##############################################################
    ##############################################################
    def create_sim(self, beta_vector, vacuum=False):

        args=self.args
        sx=self.cell_size.x

        wvg=mp.Block(center=origin, material=mp.Medium(epsilon=args.eps_wvg),
                     size=mp.Vector3(self.cell_size.x,args.w_wvg))
        cloak=mp.Cylinder(center=self.design_center, radius=args.r_cloak,
                          epsilon_func=ParameterizedDielectric(self.design_center,
                                                               self.basis,
                                                               beta_vector))
        disc=mp.Cylinder(center=self.design_center, radius=args.r_disc,
                         material=(mp.metal if args.eps_disc==0 else
                                   mp.Medium(epsilon=args.eps_disc)))

        geometry=[wvg] if vacuum else [wvg, cloak, disc]

        envelope = mp.GaussianSource(args.fcen,fwidth=args.df)
        amp=1.0
        if callable(getattr(envelope, "fourier_transform", None)):
            amp /= envelope.fourier_transform(args.fcen)
        sources=[mp.EigenModeSource(src=envelope,
                                    center=self.source_center,
                                    size=self.source_size,
                                    eig_band=self.args.source_mode,
                                    amplitude=amp
                                   )
                ]

        sim=mp.Simulation(resolution=args.res, cell_size=self.cell_size,
                          boundary_layers=[mp.PML(args.dpml)], geometry=geometry,
                          sources=sources)

        if args.complex_fields:
            sim.force_complex_fields=True

        return sim

######################################################################
# if executed as a script, we look at our own filename to figure out
# the name of the class above, create an instance of this class called
# opt_prob, and call its run() method.
######################################################################
if __name__ == '__main__':
    opt_prob=globals()[__file__.split('/')[-1].split('.')[0]]()
    opt_prob.run()
