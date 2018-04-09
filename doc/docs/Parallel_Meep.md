---
# Parallel Meep
---

Meep supports distributed-memory parallelism using [MPI](https://en.wikipedia.org/wiki/MPI). This allows it to scale up from small dual-processor machines to supercomputers, and to work on very large problems that may not fit into the memory of one machine. Meep simulations can use hundreds of processors. Of course, your problem must be sufficiently large in order to [benefit from many processors](FAQ/#should-i-expect-linear-speedup-from-the-parallel-meep).

Installing Parallel Meep
------------------------

To install the parallel version of Meep, you must have a version of MPI installed on your system. See [Installation](Installation/#mpi).

We also strongly recommend installing the [HDF5 package](Installation/#hdf5) with parallel I/O support if you are going to run with more than a few processors. HDF5 needs to be configured with the flag `--enable-parallel`. You may also have to set the `CC` environment variable to `mpicc`. Unfortunately, the parallel HDF5 library then does not work with serial code, so you have may have to choose to install either the serial or the parallel Meep, but not both.

If you don't install HDF5 with parallel I/O support, you can still do I/O from MPI &mdash; Meep has some hacks to let it write HDF5 files using serial I/O from multiple processes, one at a time. However, this does not scale very well to many processors. We've also observed some MPI implementations to freeze under the strain of trying to write from many processes at once.

Then you just configure Meep with the flag `--with-mpi`. If you run the resulting `meep` executable as usual, it runs on a single process; to run with multiple processors you should use `mpirun` as described below. Because you can run the parallel Meep in a single process like this, there is no need to separately compile and install the serial version of Meep.

Using Parallel Meep
-------------------

The parallel version of Meep is designed to operate completely transparently: you use the same Python/Scheme file as for the serial version, and the output is the same but it is just faster.

In order to run MPI programs, you typically have to use a command like `mpirun` with an argument to say how many processes you want to use. Consult your MPI documentation. For example, with many popular MPI implementations, to run with 4 processes you would use something like:

*Python*
```sh
mpirun -np 4 python foo.py > foo.out
```

*Scheme*
```sh
mpirun -np 4 meep foo.ctl > foo.out
```

There is one important requirement: every MPI process must be able to read the `foo.py`/`foo.ctl` input file or whatever your control file is called. On most systems, this is no problem, but if for some reason your MPI processes don't all have access to the local filesystem then you may need to make copies of your input file or something.

You cannot run Meep interactively on multiple processors.

### Different Forms of Parallelization

Parallel Meep works by taking your simulation and dividing the computational cell among the MPI processes. This is the only way of parallelizing a single simulation and enables simulating very large problems.

However, there is an alternative strategy for parallelization. If you have many smaller simulations that you want to run, say for many different values of some parameter, then you can just run these as separate jobs. Such parallelization is known as [embarrassingly parallel](https://en.wikipedia.org/wiki/Embarrassingly_parallel) because no communication is required. Meep provides no explicit support for this mode of operation, but of course it is quite easy to do yourself: just launch as many Meep jobs as you want, perhaps changing the parameters via the command-line using a shell script.

Technical Details
-----------------

When you run Meep under MPI, the following is a brief description of what is happening behind the scenes. For the most part, you shouldn't need to know this stuff. Just use the same Python/Scheme file exactly as you would for a uniprocessor simulation.

First, every MPI process executes the Python/Scheme file in parallel. The processes communicate however, to only perform one simulation in sync with one another. In particular, the computational cell is divided into "chunks", one per process, to roughly equally divide the work and the memory. For additional details, see Section 2.2 ("Grid chunks and owned points") of [Computer Physics Communications, Vol. 181, pp. 687-702, 2010](http://ab-initio.mit.edu/~oskooi/papers/Oskooi10.pdf).

When you time-step via Python's `meep.Simulation.run(until=...)` or Scheme's `run-until`, etc., the chunks are time-stepped in parallel, communicating the values of the pixels on their boundaries with one another. In general, any Meep function that performs some collective operation over the whole computational cell or a large portion thereof is parallelized, including: time-stepping, HDF5 I/O, accumulation of flux spectra, and field integration via `integrate_field_function` (Python) or `integrate-field-function` (Scheme), although the *results* are communicated to all processes.

Computations that only involve isolated points, such as `get_field_point` (Python) or `get-field-point` (Scheme), or `Harminv` (Python) or `harminv` (Scheme) analyses, are performed by all processes redundantly. In the case of `get_field_point` or `get-field-point`, Meep figures out which process "knows" the field at the given field, and then sends the field value from that process to all other processes. This is harmless because such computations are rarely performance bottlenecks.

Although all processes execute the Python/Scheme file in parallel, print statements are ignored for all process but one (process \#0). In this way, you only get one copy of the output.

If for some reason you need to distinguish different MPI processes in your Python/Scheme file, you can use the following two functions:

**`(meep-count-processors)`**  
**`meep.count_processors()`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Returns the number of processes that Meep is using in parallel.

**`(meep-my-rank)`**  
**`meep.my_rank()`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Returns the index of the process running the current file, from zero to (meep-count-processors)–1.

**`(meep-all-wait)`**  
**`meep.all_wait()`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Blocks until all processes execute this statment (MPI_Barrier).

**Warning**: do not attempt to perform different Meep commands in different processes by using the `(meep-my-rank)` or `meep.my_rank()`. All processes must for the most part execute the same Meep commands in the same sequence or they will deadlock, waiting forever for one another.

For large multicore jobs with I/O, it may be necessary to have `(meep-all-wait)` as the last line in the Scheme file to ensure that all processors terminate at the same point in the execution. Otherwise, one processor may finish and abruptly terminate the other processors.
