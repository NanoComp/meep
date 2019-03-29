######################################################################
# ParallelDesignTester is a simple facility for running large numbers
# of single-design-point calculations on a meep.adjoint optimization
# problem, using a multi-process server pool to allow multiple jobs
# to be run simultaneously.
#
# ** Usage as a command-line script:
#
# % python ParallelDesignTester.py [ --nproc XX ]  --name MyProb.py --casefile CaseFile
#
#     where:
#
#       -- MyProb.py is a python script describing a meep.adjoint
#          optimization problem, which implements a class MyProb
#          that inherits from meep.adjoint.OptimizationProblem
#
#       -- CaseFile is a text file in which each line specifies
#          a separate set of command-line arguments for MyProb.py,
#          describing an objective-function evaluation at a single
#          point in design space
#
#       -- The optional --nproc argument specifies the number of
#          server processes created, i.e. the number of simultaneous
#          jobs that are run. The default value is half the number
#          of CPU cores available on the machine.
#
# ** Usage from other scripts or a console/notebook interface:
#
# from meep.adjoint import ParallelDesignTester
#
# ParallelDesignTester('MyProb',command_lines).launch(nproc=nproc)
#
#  (where command_lines=['command_line_1','command_line_2',...])
#
#
######################################################################
import os
from os import environ as env
import argparse
import multiprocessing
from multiprocessing import cpu_count as ncpus
from datetime import datetime as dt
import importlib
from importlib import util
from importlib import import_module as imp_mod

######################################################################
######################################################################
######################################################################
class ParallelDesignTester(object):

    def __init__(self, name, cmdlines):
        self.name, self.cmdlines = name, cmdlines

    def __call__(self, n):
        getattr(imp_mod(self.name),self.name)(self.cmdlines[n]).run()

    def launch(self, nproc=None):

        # create and cd to working directory, then launch servers and run jobs
        wdir = '{}_{}'.format(self.name,dt.now().strftime("%m%d.%H%M%S"))
        os.mkdir(wdir)
        os.chdir(wdir)
        nproc=nproc if nproc else max(ncpus()//2,1)
        multiprocessing.Pool(nproc).map(self,range(len(self.cmdlines)))

######################################################################
######################################################################
######################################################################
if __name__ == '__main__':
    parser, nproc = argparse.ArgumentParser(), multiprocessing.cpu_count()//2
    parser.add_argument('--name',     type=str, default=None,  help='name of OptimizationProblem subclass')
    parser.add_argument('--casefile', type=str, default=None,  help='list of command-line strings')
    parser.add_argument('--nproc',    type=int, default=nproc, help='number of server processes to launch')
    args=parser.parse_args()

    if args.name is None or not importlib.util.find_spec(args.name):
       raise ValueError("missing or invalid class/module --name " + args.name if args.name else '')
    if args.casefile is None or not os.path.isfile(args.casefile):
        raise ValueError("missing or invalid list of command-line strings --casefile " + args.casefile if args.casefile else '')
    with open(args.casefile) as f:
        cmdlines=[l for l in [ln.strip(' \t\n') for ln in f.readlines()] if l and l[0]]
    ParallelDesignTester(args.name,cmdlines).launch(args.nproc)
