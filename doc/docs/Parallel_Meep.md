---
# Parallel Meep
---

Meep supports distributed-memory parallelism using [MPI](https://en.wikipedia.org/wiki/MPI). This allows it to scale up from small dual-processor machines to large parallel supercomputers, and to work on very large problems that may not even fit into the memory of one machine. We've run it using at least hundreds of processors. Of course, your problem must be sufficiently large in order to benefit from many processors.

Installing Parallel Meep
------------------------

To install the parallel version of Meep, you must have a version of MPI installed on your system. See the [Installation manual](Installation.md#mpi-parallel-machines).

We also *strongly* recommend installing the HDF5 with parallel I/O support if you are going to run with more than a few processors. (configure HDF5 with `--enable-parallel`; you may also have to set the `CC` environment variable to `mpicc`.) Unfortunately, the parallel HDF5 library then does not work with serial code, so you have may have to choose to install either the serial or the parallel Meep, but not both.

If you don't install HDF5 with parallel I/O support, you can still do I/O from MPI — Meep has some hacks to let it write HDF5 files using serial I/O from multiple processes, one at a time. However, this does not scale very well to many processors. (We've also observed some MPI implementations to freeze under the strain of trying to write from many processes at once. YMMV).

Then you just `configure` Meep `--with-mpi`. The `meep` executable is installed as `meep-mpi`, so that you can have both the serial and parallel versions installed at the same time.

Using Parallel Meep
-------------------

The parallel version of Meep is designed to operate completely transparently: you use the same `.ctl` file as for the serial version, and the output is the same, but it is just faster (hopefully).

In order to run MPI programs, you typically have to use a command like `mpirun` with an argument to say how many processes you want to use. (Consult your MPI documentation.) For example, with many popular MPI implementations, to run with 4 processes you would use something like:

```
mpirun -np 4 meep-mpi foo.ctl > foo.out
```


There is one important requirement: every MPI process must be able to read the `foo.ctl` input file (or whatever your control file is called). On most systems, this is no problem, but if for some reason your MPI processes don't all have access to the local filesystem then you may need to make copies of your input file or something.

(You cannot run Meep interactively on multiple processors.)

### Different forms of parallelization

Parallel Meep works by taking your simulation and dividing the computational cell among the MPI processes. This is the only way of parallelizing a single simulation, and allows you to attack very large problems.

However, there is an alternative strategy for parallelization. If you have many smaller simulations that you want to run, say for many different values of some parameter, then you can just run these as separate jobs. Such parallelization is known as "embarrassingly parallel" because no communication is required. Meep provides no explicit support for this mode of operation, but of course it is quite easy to do yourself: just launch as many Meep jobs as you want, perhaps changing the parameters via the command-line using a shell script.

Technical Details
-----------------

When you run Meep under MPI, the following is a brief description of what is happening behind the scenes. For the most part, you shouldn't *need* to know this stuff...just use the same `.ctl` exactly as you would for a uniprocessor simulation.

First, every MPI process executes the `.ctl` file in parallel. The processes communicate however, to only perform one simulation in sync with one another. In particular, the computational cell is divided into "chunks", one per process, to roughly equally divide the work and the memory.

When you time-step (via `run-until` or whatever), the chunks are time-stepped in parallel, communicating the values of the pixels on their boundaries with one another. In general, any Meep function that performs some collective operation over the whole computational cell or a large portion thereof is parallelized, including: time-stepping, HDF5 I/O, accumulation of flux spectra, and field integration (via `integrate-field-function`, although the *results* are communicated to all processes).

Computations that only involve isolated points, such as `get-field-point` or `harminv` analyses, are performed by all processes redundantly. (In the case of `get-field-point`, Meep figures out which process "knows" the field at the given field, and then sends the field value from that process to all other processes.) This is harmless because such computations are rarely performance bottlenecks.

Although all processes execute the `.ctl` in parallel, `print` statements are ignored for all process but one (process \#0). In this way, you only get one copy of the output.

If for some reason you need to distinguish different MPI processes in your `.ctl` file, you can use the following two functions:

(meep-count-processors)
Returns the number of processes that Meep is using in parallel.

(meep-my-rank)
Returns the index of the process running the current `.ctl` file, from zero to `(meep-count-processors)` – 1.

**Warning**: do not attempt to perform different Meep commands in different processes by using the `(meep-my-rank)`. All processes must (for the most part) execute the same Meep commands in the same sequence or they will deadlock, waiting forever for one another.

For large multicore jobs with I/O, it may be necessary to have `(meep-all-wait)` as the last line in the `.ctl` file to ensure that all processors terminate at the same point in the execution. Otherwise, one processor may finish and abruptly terminate the other processors.
