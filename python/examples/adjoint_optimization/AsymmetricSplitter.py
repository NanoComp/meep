import sys
import argparse
import numpy as np
import meep as mp

from meep.adjoint import (OptimizationProblem, DFTCell, adjoint_options,
                          xHat, yHat, zHat, origin, FluxLine,
                          ParameterizedDielectric, FiniteElementBasis)

##################################################
##################################################
##################################################
class Splitter12(OptimizationProblem):

    ##################################################
    ##################################################
    ##################################################
    def add_args(self, parser):

        # add new problem-specific arguments
        parser.add_argument('--dair',        type=float, default=-1.0, help='')
        parser.add_argument('--w_in',        type=float, default=1.0,  help='width of input waveguide')
        parser.add_argument('--w_out1',      type=float, default=0.5,  help='width of output waveguide 1')
        parser.add_argument('--w_out2',      type=float, default=0.5,  help='width of output waveguide 2')
        parser.add_argument('--l_stub',      type=float, default=3.0,  help='length of waveguide input/output stub')
        parser.add_argument('--l_design',    type=float, default=2.0,  help='length of design region')
        parser.add_argument('--h_design',    type=float, default=6.0,  help='height of design region')
        parser.add_argument('--eps_in',      type=float, default=6.0,  help='input waveguide permittivity')
        parser.add_argument('--eps_out1',    type=float, default=2.0,  help='output waveguide 1 permittivity')
        parser.add_argument('--eps_out2',    type=float, default=12.0, help='output waveguide 2 permittivity')
        parser.add_argument('--nfe',         type=int,   default=2,    help='number of finite elements per unit length')

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
        dair       = 0.5*args.w_in if args.dair==-1.0 else args.dair
        sx         = dpml + args.l_stub + args.l_design + args.l_stub + dpml
        sy         = dpml + dair + args.h_design + dair + dpml
        cell_size  = mp.Vector3(sx, sy, 0.0)

        #----------------------------------------
        #- design region
        #----------------------------------------
        design_center = origin
        design_size   = mp.Vector3(args.l_design, args.h_design, 0.0)
        design_region = mp.Volume(center=design_center, size=design_size)

        #----------------------------------------
        #- objective regions
        #----------------------------------------
        x_in          =  -0.5*(args.l_design + args.l_stub)
        x_out         =  +0.5*(args.l_design + args.l_stub)
        y_out1        =  +0.25*args.h_design
        y_out2        =  -0.25*args.h_design

        flux_in       =  FluxLine(x_in,     0.0, 2.0*args.w_in,   mp.X, 'in')
        flux_out1     =  FluxLine(x_out, y_out1, 2.0*args.w_out1, mp.X, 'out1')
        flux_out2     =  FluxLine(x_out, y_out2, 2.0*args.w_out2, mp.X, 'out2')

        objective_regions  = [flux_in, flux_out1, flux_out2]

        #----------------------------------------
        #- optional extra regions for visualization if the --full-dfts options is present.
        #----------------------------------------
        extra_regions      = [mp.Volume(center=origin, size=cell_size)] if args.full_dfts else []

        #----------------------------------------
        # forward source region
        #----------------------------------------
        source_center    =  (x_in - 0.25*args.l_stub)*xHat
        source_size      =  2.0*args.w_in*yHat

        #----------------------------------------
        # basis set
        #----------------------------------------
        basis = FiniteElementBasis(args.l_design, args.h_design, args.nfe)

        #----------------------------------------
        #- objective function
        #----------------------------------------
        fstr = (   'Abs(P1_out1)**2'
                 + '+0.0*(P1_out1 + M1_out1)'
                 + '+0.0*(P1_out2 + M1_out2)'
                 + '+0.0*(P1_in   + M1_in + S_out1 + S_out2 + S_in)'
               )

        #----------------------------------------
        #- internal storage for variables needed later
        #----------------------------------------
        self.args            = args
        self.dpml            = dpml
        self.cell_size       = cell_size
        self.basis           = basis
        self.design_center   = origin
        self.design_size     = design_size
        self.source_center   = source_center
        self.source_size     = source_size

        return fstr, objective_regions, extra_regions, design_region, basis

    ##############################################################
    ##############################################################
    ##############################################################
    def create_sim(self, beta_vector, vacuum=False):

        args=self.args
        sx=self.cell_size.x

        x_in   = -0.5*(args.l_design + args.l_stub)
        x_out  = +0.5*(args.l_design + args.l_stub)
        y_out1 = +0.25*args.h_design
        y_out2 = -0.25*args.h_design

        wvg_in = mp.Block( center=mp.Vector3(x_in,0.0),
                           size=mp.Vector3(args.l_stub,args.w_in),
                           material=mp.Medium(epsilon=args.eps_in))
        wvg_out1 = mp.Block( center=mp.Vector3(x_out,y_out1),
                             size=mp.Vector3(args.l_stub,args.w_out1),
                             material=mp.Medium(epsilon=args.eps_out1))
        wvg_out2 = mp.Block( center=mp.Vector3(x_out,y_out2),
                             size=mp.Vector3(args.l_stub,args.w_out2),
                             material=mp.Medium(epsilon=args.eps_out2))
        design   = mp.Block( center=origin,
                             size=mp.Vector3(args.l_design,args.h_design),
                             epsilon_func=ParameterizedDielectric(self.design_center,
                                                                  self.basis,
                                                                  beta_vector)
                           )

        geometry=[wvg_in, wvg_out1, wvg_out2, design]

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
