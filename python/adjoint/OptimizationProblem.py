import os
import sys
import argparse
import pickle
from datetime import datetime as dt2
from collections import namedtuple

import numpy as np

from .Objective import (adjoint_options, unit_vector, log,
                        DFTCell, ObjectiveFunction, AdjointSolver)

from .Visualization import (AdjointVisualizer, set_plot_default, def_plot_options)

######################################################################
# invoke python's 'abstract base class' formalism in a version-agnostic way
######################################################################
from abc import ABCMeta, abstractmethod
ABC = ABCMeta('ABC', (object,), {'__slots__': ()}) # compatible with Python 2 *and* 3:

######################################################################
# OptimizationProblem is an abstract base class from which classes
# describing specific optimization problems should inherit
######################################################################
class OptimizationProblem(ABC):

    ######################################################################
    # pure virtual methods that must be overridden by child
    ######################################################################
    @abstractmethod
    def init_problem(self, args):
        raise NotImplementedError("derived class must implement init_problem() method")

    @abstractmethod
    def create_sim(self, betavector, vacuum=False):
        raise NotImplementedError("derived class must implement create_sim() method")

    ######################################################################
    # virtual methods that *may* optionally be overridden by child
    ######################################################################
    def add_args(self, args):
        pass

    ######################################################################
    ######################################################################
    ######################################################################
    def __init__(self, cmdline=None):


        # define and parse command-line arguments
        parser = self.init_args()       # general args
        self.add_args(parser)           # problem-specific arguments

        # set options from cmdline parameter (if present) or sys.argv
        argv =        [ a for a in cmdline.split(' ') if a ] if cmdline \
                 else sys.argv[1:] if len(sys.argv)>1 else []
        self.cmdline = ' '.join(argv)
        self.args = args = parser.parse_args(argv)

        # call subclass for problem-specific initialization first...
        fstr, objective_regions, extra_regions, design_region, self.basis \
            = self.init_problem(self.args)

        # and now do some general initialization of basic data fields:
        #    DFT cells
        fcen, df, nfreq         = args.fcen, args.df, args.nfreq
        DFTCell.reset_cell_names()
        self.objective_cells    = [ DFTCell(region=v, fcen=fcen, df=df, nfreq=nfreq) for v in objective_regions ]
        self.extra_cells        = [ DFTCell(region=v, fcen=fcen, df=df, nfreq=nfreq) for v in extra_regions ] if args.full_dfts else []
        self.design_cell        = DFTCell(region=design_region, fcen=fcen, df=df, nfreq=nfreq, name='design_fields')
        self.dft_cells          = self.objective_cells + self.extra_cells + [self.design_cell]

        # objective function
        self.obj_func           = ObjectiveFunction(fstr=fstr)
        self.fqnames            = ['f'] + self.obj_func.qnames

        # design variables
        self.dim                = self.basis.dim
        self.beta_vector        = self.init_beta_vector()

        # options affecting meep calculations
        adjoint_options['dft_reltol']   = args.dft_reltol
        adjoint_options['dft_timeout']  = args.dft_timeout
        adjoint_options['dft_interval'] = args.dft_interval
        adjoint_options['verbosity']    = 'verbose' if args.verbose         \
                                          else 'concise' if args.concise    \
                                          else adjoint_options['verbosity']
        adjoint_options['logfile']      = args.logfile
        self.filebase                   = args.filebase if args.filebase else self.__class__.__name__

        # options controlling the optimizer
        adjoint_options['logfile']      = args.logfile

        # miscellaneous general options affecting visualization

        if args.label_source_regions:
            set_plot_default('fontsize',def_plot_options['fontsize'], 'src')

        adjoint_options['animate_components'] = args.animate_component
        adjoint_options['animate_interval']   = args.animate_interval

        # other data structures that are initialized on a just-in-time basis
        self.solver = self.sim = self.vis = self.dfdEps = None

    ######################################################################
    # constructor helper method that initializes the command-line parser
    #  with general-purpose (problem-independent) arguments
    ######################################################################
    def init_args(self):

        parser = argparse.ArgumentParser()

        #--------------------------------------------------
        # parameters affecting meep computations
        #--------------------------------------------------
        parser.add_argument('--res',          type=float, default=20,      help='resolution')
        parser.add_argument('--dpml',         type=float, default=-1.0,    help='PML thickness (-1 --> autodetermined)')
        parser.add_argument('--fcen',         type=float, default=0.5,     help='center frequency')
        parser.add_argument('--df',           type=float, default=0.25,    help='frequency width')
        parser.add_argument('--source_mode',  type=int,   default=1,       help='mode index of eigenmode source')
        parser.add_argument('--dft_reltol',   type=float, default=adjoint_options['dft_reltol'],   help='convergence threshold for end of timestepping')
        parser.add_argument('--dft_timeout',  type=float, default=adjoint_options['dft_timeout'],  help='max runtime in units of last_source_time')
        parser.add_argument('--dft_interval', type=float, default=adjoint_options['dft_interval'], help='meep time DFT convergence checks in units of last_source_time')

        #--------------------------------------------------
        # flags affecting outputs from meep computations
        #--------------------------------------------------
        parser.add_argument('--nfreq',          type=int,              default=1,           help='number of output frequencies')
        parser.add_argument('--full_dfts',      dest='full_dfts',      action='store_true', help='compute DFT fields over full volume')
        parser.add_argument('--complex_fields', dest='complex_fields', action='store_true', help='force complex fields')
        parser.add_argument('--filebase',       type=str,              default=self.__class__.__name__, help='base name of output files')

        #--------------------------------------------------
        # initial values for basis-function coefficients
        #--------------------------------------------------
        parser.add_argument('--betafile',  type=str,   default='',    help='file of expansion coefficients')
        parser.add_argument('--beta',       nargs=2,   default=[],    action='append',  help='set value of expansion coefficient')
        parser.add_argument('--eps_design', type=str,  default=None,  help='functional expression for initial design permittivity')

        #--------------------------------------------------
        # options describing the calculation to be done
        #--------------------------------------------------
        # do a single calculation of the objective function, optionally with gradient
        parser.add_argument('--eval_objective',  dest='eval_objective',  action='store_true', help='evaluate objective function value')
        parser.add_argument('--eval_gradient',   dest='eval_gradient',   action='store_true', help='evaluate objective function value and gradient')
        parser.add_argument('--gradient_qname',  type=str, default=None, help='name of objective quantity to differentiate via adjoint method')

        # compute finite-difference approximation to derivative for test purposes
        parser.add_argument('--fd_order',      type=int,   default=0,        help='finite-difference order (0,1,2)')
        parser.add_argument('--fd_index',      default=[], action='append',  help='index of differentiation variable')
        parser.add_argument('--fd_rel_delta',  type=float, default=1.0e-2,   help='relative finite-difference delta')

        # run the full iterative optimization
        parser.add_argument('--optimize',       dest='optimize', action='store_true', help='perform automated design optimization')

        #--------------------------------------------------
        #- options affecting optimization ---------------
        #--------------------------------------------------
        parser.add_argument('--alpha',          type=float,     default=1.0, help='gradient descent relaxation parameter')
        parser.add_argument('--min_alpha',      default=1.0e-3, help='minimum value of alpha')
        parser.add_argument('--max_alpha',      default=10.0,   help='maximum value of alpha')
        parser.add_argument('--boldness',       default=1.25,   help='sometimes you just gotta live a little')
        parser.add_argument('--timidity',       default=0.75,   help='can\'t be too careful in this dangerous world')
        parser.add_argument('--max_iters',      type=int, default=100, help='max number of optimization iterations')
        parser.add_argument('--overlap_dfdEps', dest='overlap_dfdEps', action='store_true')

        #--------------------------------------------------
        # flags configuring adjoint-solver options
        #--------------------------------------------------l
        parser.add_argument('--verbose',      dest='verbose',   action='store_true', help='produce more output')
        parser.add_argument('--concise',      dest='concise',   action='store_true', help='produce less output')
        parser.add_argument('--visualize',    dest='visualize', action='store_true', help='produce visualization graphics')
        parser.add_argument('--label_source_regions', dest='label_source_regions', action='store_true', help='label source regions in visualization plots')
        parser.add_argument('--logfile',        type=str,   default=None,    help='log file name')
        parser.add_argument('--pickle_data',  dest='pickle_data', action='store_true', help='save state to binary data file')
        parser.add_argument('--animate_component', action='append', help='plot time-domain field component')
        parser.add_argument('--animate_interval', type=float, default=1.0, help='meep time between animation frames')

        return parser

    ######################################################################
    # constructor helper method to process command-line arguments for
    # initializing the vector of design variables
    ######################################################################
    def init_beta_vector(self):

        beta_vector=np.zeros(self.dim)

        #######################################################################
        # if a --betafile was specified, try to parse it in the form
        #    beta_0 \n beta_1 \n ... (if single column)
        # or
        #    i1	 beta_i1 \n i2 beta_i2 \n  (if two columns)
        #######################################################################
        if self.args.betafile:
            fb = np.loadtxt(self.args.betafile)   # 'file beta'
            if np.ndim(fb)==1:
                indices,vals = range(len(fb)), fb
            elif np.ndim(fb)==2 and np.shape(fb)[1] == 2:
                indices,vals = fb[:,0], fb[:,1]
            else:
                raise ValueError("{}: invalid file format".format(self.args.betafile))
            for i,v in zip(indices,vals): beta_vector[i]=v

        #######################################################################
        # parse arguments of the form --beta index value
        #######################################################################
        for ivpair in self.args.beta:     # loop over (index,value) pairs
            beta_vector[int(ivpair[0])]=ivpair[1]

        #######################################################################
        # if a functional form for --eps_design was specified, project that
        # function onto the basis and use this as the initial design point
        #######################################################################
        if np.count_nonzero(beta_vector)==0:
            eps_design=self.args.eps_design if self.args.eps_design else 1.0
            sim = self.create_sim(beta_vector)
            sim.init_sim()
            xyzw=sim.get_dft_array_metadata(center=self.dft_cells[-1].center, size=self.dft_cells[-1].size)
            beta_vector=self.basis.expand(eps_design,xyzw,cache=True,subtract_one=True)

        return beta_vector

    ######################################################################
    # terminate script without exiting the (i)python console or notebook
    ######################################################################
    def terminate(self, msg):
        raise Exception(msg)

    ######################################################################
    ######################################################################
    ######################################################################
    def plot_geometry(self):
        sim = self.create_sim(self.beta_vector)
        sim.init_sim()
        [cell.register(sim) for cell in self.dft_cells]
        vis = AdjointVisualizer(cases=['Geometry'])
        vis.update(sim,'Geometry')
        self.sim,self.vis = sim,vis # save copies for debugging; not strictly necessary
        return

    ######################################################################
    ######################################################################
    ######################################################################
    def update_design_variables(self, beta_vector):
        if self.solver==None:
            self.solver=AdjointSolver(self.obj_func, self.dft_cells, self.basis)
        self.solver.sim=self.create_sim(beta_vector)

    ######################################################################
    # compute the objective function, plus possibly its gradient (via
    # adjoints), at a single point in design space
    ######################################################################
    def eval_objective(self, beta_vector, need_gradient=False, vis=None):
        self.update_design_variables(beta_vector)
        self.solver.vis=vis
        return self.solver.solve(need_gradient=need_gradient)

    ######################################################################
    # an OptState stores the current state of an optimization problem.
    ######################################################################
    OptState = namedtuple('OptState', 'n alpha beta fq gradf dfdEps')

    ######################################################################
    ######################################################################
    ######################################################################
    def log_state(self, state, substate=None):

        ts=dt2.now().strftime('%T')
        with open(self.iterfile,'a') as f:
            if substate:
                f.write('{}:  Subiter {}.{}: '.format(ts,state.n,substate.n))
                f.write('f={}, alpha={}\n'.format(substate.fq[0],substate.alpha))
                return

            dfdEps_avg=np.sum(self.dft_cells[-1].xyzw[3] * state.dfdEps)
            f.write('\n\n{}: Iter {}: f={}, alpha={} dfdeAve={}\n'
                    .format(ts,state.n,state.fq[0],state.alpha,dfdEps_avg))
            [f.write('#{} {} = {}\n'.format(state.n,nn,qq)) for nn,qq in zip(self.fqnames[1:], state.fq[1:]) ]
            [f.write('#{} b{} = {}\n'.format(state.n,n,b)) for n,b in enumerate(state.beta)]
            f.write('\n\n')
            self.annals.append(state)

    ######################################################################
    ######################################################################
    ######################################################################
    def line_search(self,state):

        self.log_state(state)
        cease_file = '/tmp/terminate.{}'.format(os.getpid())

        bs, xyzw, ovrlp = self.basis, self.dft_cells[-1].xyzw, self.args.overlap_dfdEps
        dbeta = bs.overlap(state.dfdEps, xyzw) if ovrlp else bs.expand(state.dfdEps, xyzw)

        alpha, iter, subiter = state.alpha, state.n, 0
        while alpha>self.args.min_alpha:

            beta = state.beta + alpha*np.real(dbeta)
            for n in range(len(beta)):
                beta[n] = max(0.0, beta[n])
            self.update_design_variables(beta)
            fq   = self.solver.forward_solve()

            substate = self.OptState(subiter,alpha,beta,fq,0,0)
            self.log_state(state,substate)

            if fq[0] > state.fq[0]:    # found a new optimum, declare victory and a new iteration
                gradf = self.solver.adjoint_solve()
                alpha = min(alpha*self.args.boldness,self.args.max_alpha)
                return self.OptState(iter+1, alpha, beta, fq, gradf, self.solver.dfdEps)

            if os.path.isfile(cease_file):  # premature termination requested by user
                os.remove(cease_file)
                return None

            alpha*=self.args.timidity

        return None   # unable to improve objective by proceeding any distance in given direction

    ######################################################################
    ######################################################################
    ######################################################################
    def optimize(self):

        ss = int( dt2.now().strftime("%s") ) - 1553607629
        self.iterfile = '{}.{}.iters'.format(self.filebase,ss)
        with open(self.iterfile,'w') as f:
            f.write('#{} ran'.format(self.__class__.__name__))
            f.write(dt2.now().strftime(' %D::%T\n'))
            f.write('# with args {}\n\n'.format(self.cmdline))
        self.annals=[]

        ######################################################################
        # initialize AdjointSolver and get objective function value and gradient
        # at the initial design point
        ######################################################################
        alpha        = self.args.alpha
        beta         = self.beta_vector
        self.update_design_variables(beta)
        fq,gradf     = self.solver.solve(need_gradient=True)
        dfdEps       = self.solver.dfdEps

        state        = self.OptState(1,alpha,beta,fq,gradf,dfdEps)

        while state and state.n<self.args.max_iters:
            state = self.line_search(state)

        return self.annals

    ######################################################################
    # the 'run' class method of OptimizationProblem, which is where control
    # passes after __init__ if we were executed as a script, looks at
    # command-line flags to determine whether to (a) launch a full
    # iterative optimization of the design variables or (b) study
    # the geometry for just the one given set of design variables
    ######################################################################
    def run(self):

        #--------------------------------------------------------------
        # if --optimize is present, run iterative design optimization
        #--------------------------------------------------------------
        args=self.args
        if args.optimize:
            self.optimize()
            return

        #--------------------------------------------------------------
        # if no computation was requested, just plot the geometry
        #--------------------------------------------------------------
        if not(args.eval_gradient or args.eval_objective or args.fd_order>0):
            self.plot_geometry()
            return

        #--------------------------------------------------------------
        #--------------------------------------------------------------
        #--------------------------------------------------------------
        vis = None
        if args.visualize:
            cases=['Geometry','Forward'] + (['Adjoint'] if args.eval_gradient else [])
            vis=AdjointVisualizer(cases=cases)

        fq,gradf=self.eval_objective(self.beta_vector, need_gradient=args.eval_gradient, vis=vis)

        self.output([],[],actions=['begin'],msg='{} {}'.format(self.filebase,self.cmdline))
        self.output(['res','fcen','df','source_mode'], [args.res,args.fcen,args.df,args.source_mode])
        [self.output(['beta'+str(n)],[beta]) for (n,beta) in enumerate(self.beta_vector)]
        [self.output([fqn], [fqv]) for (fqn,fqv) in zip(self.fqnames,fq)]
        if args.eval_gradient:
            [self.output( ['df/db'+str(n)], [g]) for n,g in enumerate(gradf)]

        #--------------------------------------------------------------
        # calculate the objective function value, optionally its gradient,
        # and optionally the finite-difference approximation to the
        # derivative of a single objective quantity
        #--------------------------------------------------------------
        indices = [] if args.fd_order==0 else [int(index) for index in args.fd_index]
        for index in indices:
            dbeta    = args.fd_rel_delta*(1.0 if self.beta_vector[index]==0 else np.abs(self.beta_vector[index]))
            beta_hat = unit_vector(index,self.dim)
            fqp, _ = self.eval_objective(self.beta_vector + dbeta*beta_hat)
            d1fq   = (fqp - fq) / dbeta
            if args.fd_order==2 and self.beta_vector[index]!=0.0:
                fqm, _ = self.eval_objective(self.beta_vector - dbeta*beta_hat)
                d2fq   = (fqp - fqm) / (2.0*dbeta)
            elif args.fd_order==2 and self.beta_vector[index]==0.0:
                fqpp, _ = self.eval_objective(self.beta_vector + 2*dbeta*beta_hat)
                d2fq   = (4*fqp - fqpp - 3.0*fq) / (2.0*dbeta)
            for i,fqn in enumerate(self.fqnames):
                nlist = ['d{}/db{}__O{}__'.format(fqn,index,ord+1) for ord in range(args.fd_order)]
                vlist = [ d1fq[i] ] + ([d2fq[i]] if args.fd_order==2 else[])
                [self.output(nlist, vlist)]

        self.output([],[],actions=['end'])

        if args.pickle_data:
            f = open(self.filebase + '.pickle', 'wb')
            pickle.dump(self.solver.dfdEps,f)
            f.close()

        return

    ######################################################################
    # 'lcdoi' = legend, console, digest, output, iterations
    ######################################################################
    def output(self, names, values, actions=[], msg=None, files='lcdoi'):

        if 'begin' in actions:
            self.nout = 1 # running index of output quantity
            self.legend  = open(self.filebase + '.legend','w') if 'l' in files else None
            self.digest  = open(self.filebase + '.digest','a') if 'd' in files else None
            self.outfile = open(self.filebase + '.out','a')    if 'o' in files else None
            self.iterlog = open(self.filebase + '.iters','a')  if 'i' in files else None

        msgfiles  = [sys.stdout]  if 'c' in files else []
        msgfiles += [self.digest] if 'd' in files else []
        if msg is not None:
            tm=dt2.now().strftime('%D::%T ')
            for f in msgfiles:
                f.write('\n\n** {} {}\n**'.format(tm,msg))

        #--------------------------------------------------------------
        #- inline utility functions for real/complex numerical output
        #--------------------------------------------------------------
        def fw(f,s=None):
            if f is not None:
                f.write(s) if s is not None else f.close()

        def myangle(z):
            theta=np.angle(z,deg=True)
            return theta if theta>=0.0 else theta+360.0

        def pretty_print(val,for_file=False,polar=False):
            fmt='{:+.8e} ' if for_file else '{:+.4e}'
            if not np.iscomplex(val):
                return fmt.format(np.real(val))
            (v1,v2)=(np.abs(val),myangle(val)) if polar else (np.real(val), np.imag(val))
            if for_file:
                return (fmt + fmt).format(v1,v2)
            elif polar:
                return (fmt + '@' + '{:.0e}').format(v1,v2)
            else:
                return ('(' + fmt + ',' + fmt + ')').format(v1,v2)

        namestrs, valstrs = '', ''
        for name,val in zip(names, values):

            namestrs += (' ,' if namestrs else '') + name
            valstrs  += (' ,' if valstrs  else '') + pretty_print(val)
            if np.iscomplex(val):
                valstrs += '(' + pretty_print(val,polar=True) + ' '

            fw(self.outfile,pretty_print(val, for_file=True))

            idxstr=str(self.nout) + (','+str(self.nout+1) if np.iscomplex(val) else '')
            fw(self.legend,'{}: {}\n'.format(idxstr,name))
            self.nout += (2 if np.iscomplex(val) else 1) if self.legend else 0

        for f in msgfiles:
            f.write('{:30s}:  {}\n'.format(namestrs,valstrs))

        if 'end' in actions:
            fw(self.outfile,'\n')
            fw(self.outfile)
            fw(self.digest)
            fw(self.legend)
