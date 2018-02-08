---
# C++ Tutorial
---

Instead of using the [Python interface](Python_User_Interface/), Meep is also callable as a C++ library by writing a C++ program that links to it. The C++ interface provides the most flexibility in setting up simulations. There are numerous examples in the `tests/` directory of the Meep codebase which cover a wide range of Meep's functionality. These are a good additional reference. Finally, we should also note that, while Meep is nominally in C++, it is perhaps better described as "C+". That is, most of the coding style is C-like with a few C++ features.

[TOC]

Differences from libctl
-----------------------

The C++ interface has several differences from the libctl interface besides the obvious difference in syntax.

The most notable difference is that, while the libctl interface puts the origin (0,0,0) at the *center* of the computational cell, the C++ interface by default puts the origin at the *corner* of the computational cell. That is, an $L\times L\times L$ cell goes from ($-L/2$,$-L/2$,$-L/2$) to ($L/2$,$L/2$,$L/2$) in libctl, but from (0,0,0) to ($L$,$L$,$L$) in C++. This can be changed by calling `grid_volume::shift_origin`.

Overview
--------

We begin with a brief outline of a Meep C++ program, with minimal explanations, leaving more details for the examples below and the [C++ interface](C++_User_Interface.md). The C++ program should begin with:

```c++
#include <meep.hpp>
using namespace meep;
```

to include the definitions of the Meep routines. Later, when you compile (see below), you must also link to the Meep libraries.

Your main program should then initialize Meep, and will generally then define a computational `grid_volume` and the associated `structure` describing the geometry and materials, initialize the `fields`, add sources, and then time-step. In short:

```c++
int main(int argc, char **argv) {
  initialize mpi(argc, argv); // do this even for non-MPI Meep
  double resolution = 20;     // pixels per distance
  grid_volume v = vol2d(5, 10, resolution); // 5x10 2d cell
  structure s(v, eps, pml(1.0));
  fields f(&s);

  f.output_hdf5(Dielectric, v.surroundings());

  double freq = 0.3, fwidth = 0.1;
  gaussian_src_time src(freq, fwidth);
  f.add_point_source(Ey, src, vec(1.1, 2.3));
  while (f.time() < f.last_source_time()) {
     f.step();
  }

  f.output_hdf5(Hz, v.surroundings());

  return 0;
}
```

This example doesn't do much &mdash; it just runs a Gaussian source and outputs the $H_z$ field at the end. The dielectric structure is determined by the user-defined function `eps`, which has the form:

```c++
double eps(const vec &p) {
  if (p.x() < 2 && p.y() < 3)
    return 12.0;
  return 1.0;
}
```

which returns the dielectric function $\varepsilon(\mathbf{x})$ which is just a 2$\times$3 rectangle of $\varepsilon$=12 in the upper-left corner. Unlike in the Scheme interface, by default the origin of the coordinate system is at the *corner* of the computational cell.

Now that you have the basic flavor, we can proceed to some more specific examples.

Computing the Quality Factor of a Resonator
-------------------------------------------

In this first tutorial, we will write the script to compute the quality factor of a 1d Fabry-Perot cavity. For a 1d system, Meep considers a computational cell along the $z$ coordinate. The control file will be a C++ file having extension \*.cpp. In order to use all the classes and subroutines available in Meep, the first two lines of any control file must be the following:

```c++
#include <meep.hpp>
usingnamespace meep;
```

The particular Fabry-Perot cavity we will investigate consists of an air region bounded by two distributed Bragg reflectors which are quarter-wave stacks of $\varepsilon$ of 12 and 1. We choose the size of the defect to be twice as large as the air thickness in the quarter-wave stack, so that a defect mode is found near midgap. This structure will include $N$ = 5 periods in the Bragg reflector on either side of the defect. We can use a larger $N$ but the quality factor may then be too large to compute. The parameters are set up as follows:

```c++
const double eps1 = 12.0; // epsilon of layer1
const double eps2 = 1.0;  // epsilon of layer 2
const double grating_periodicity = 1.0;
const double d1 = sqrt(eps2) / (sqrt(eps1)+sqrt(eps2)); // quarter wave stack dimensions
const double d2 = sqrt(eps1) / (sqrt(eps1)+sqrt(eps2));
const double half_cavity_width = d2;
const int N = 5;
```

Meep supports perfectly matching layers (PML) as absorbing boundary conditions. The PML begins at the edge of the computational volume and works inwards. Hence, we specify the size of the computational cell as follows: 

```c++
const double pml_thickness = 1.0;
const double z_center = half_cavity_width + N*grating_periodicity + pml_thickness;
```

Note that `z_center` is half the cell length. To specify a dielectric structure, we define a function that takes as input one parameter, a position vector, and returns the value of the dielectric at that position.

```c++
double eps(const vec &p) {
 vec r = p - vec(z_center);
 if (abs(r.z()) < half_cavity_width)
   return 1.0;
 else {
   double dz = abs(r.z()) - half_cavity_width;
   while (dz > grating_periodicity)
     dz -= grating_periodicity;
   return (dz < d1) ? eps1 : eps2;
 }
}
```

We are now ready to set up the computational cell, excite the sources, time step the fields, and compute the resulting quality factor. Here we set up the main part of the control file incorporating some of Meep's classes and sub routines. We will excite the mode at the midgap of the bandgap, which we expect to have the largest quality factor since it will have the largest exponential decay.

```c++
int main(int argc, char **argv) {
 initialize mpi(argc,argv);
 const double amicron = 10.0;
 const grid_volume vol = vol1d(2*z_center,amicron);
 structure s(vol,eps,pml(pml_thickness));
 fields f(&s);
```

Note the constructor for the `grid_volume` class in 1d which takes as parameters the size of the cell and the resolution. The sources in Meep are excitations of the polarization vector in $\mathbf{D}=\varepsilon\mathbf{E}+\mathbf{P}$. The polarization can be any one of the six cartesian or cylindrical fields. There are a variety of sources including dipole and current sources, gaussian pulses and a continuous wave sources. For now, we use a gaussian source centered at the midgap frequency with a narrow width, along with the start time and end time of the source.

```c++
 const double w_midgap = (sqrt(eps1)+sqrt(eps2))/(4*sqrt(eps1)*sqrt(eps2));
 gaussian_src_time src(w_midgap, 1.0, 0.0, 5.0);
 f.add_point_source(Ex, src, vol.center());
```

Here we make use of the built in function, `vol.center()`, that automatically computes the geometric center of the of the computational volume. To see what the unit cell of the dielectric function looks, output as an HDF5 file, we invoke the built in command:

```c++
f.output_hdf5(Dielectric, v.surroundings())
```

This versatile command can be used to visualize all of the field components, field energy when called at a particular time step. The resulting output file will be `eps-000000.00.h5` and to view it, we first need to convert it to the portable network graphics (PNG) format with the `h5topng` tool. We are now ready to begin the time stepping of the fields, and do so in a simple loop:

```c++
while (f.time() < f.last_source_time()) f.step();
```

To calculate the quality factor, we need to monitor a field beyond the time at which our point source has been turned off. We can do this in Meep by establishing an array where all fields at a particular point can be monitored. We set up the area as well as the output arrays we will need later as such:

```c++
int maxbands = 5;
int ttot = int(400/f.dt)+1;
complex`<double>` *p = new complex`<double>`[ttot];
complex`<double>` *amps = new complex`<double>`[maxbands];
double *freq_re = new double[maxbands];
double *freq_im = new double[maxbands];
```

The variable `maxbands` is the number of mode frequencies we will be looking for while `ttot` represents the total number of timesteps for which we will monitor the field component at. For this simulation given the large quality factor that we expect, we will time step for a long period of time which in this case is 400 unit seconds. Now we are ready to time step and collect the field information:

```c++
 int i=0;
 while (f.time() < f.last_source_time() + 400) {
   f.step();
   p[i++] = f.get_field(Ex,vol.center());
}
```

The `get_field` function does exactly that, obtains the value of the field at a given position within the computational cell using linear interpolation where necessary. Now we have all the information we need and must obtain the frequencies of the modes which we do by invoking a special tool for harmonic inversion, `harminv`, to extract the real and imaginary frequencies:

```c++
 int num = do_harminv(p, ttot, f.dt, 0.8*w_midgap, 1.2*w_midgap, maxbands, amps, freq_re, freq_im);
```

The integer returned by `harminv` is the number of mode frequencies obtained from the input data `p`. The particular call to `harminv` included passing the arrays by value and telling `harminv` to look for frequencies within 20% of the mid gap frequency up to a maximum of 5 bands. At this point, the necessary information to compute the quality has been stored in the `freq_re` and `freq_im` arrays and we compute the quality factor using the formula, $Q=-\omega_r/2\omega_i$. All that is left to do, is to print these results with the quality factor:

```c++
master_printf("frequency,amplitude,quality factor\n");
for (int i=0; i<num; ++i)
  master_printf("%0.6g%+0.6gi,%0.5g%+0.5gi,%g\n",freq_re[i],freq_im[i],real(amps[i]),imag(amps[i]),-0.5*freq_re[i]/freq_im[i]);
return 0;
}
```

That's it, we are done. The output for this program is `3.26087e+06`. The large quality factor is to be expected since our system includes a relatively high number of bragg reflectors ($N$=5) in each direction and we are exciting the mode at mid gap. Suppose we want to see what the resonant mode looks like, say over one period. The following block of code time steps the field for exactly one period while outputting the $E_x$ field for the computational cell at each time step:

```c++
double curr_time = f.time();
while (f.time() < curr_time + 1/w_midgap) {
 f.output_hdf5(Ex, vol.surroundings());
 f.step();
}
```

After we convert the hdf5 files to PNG format and superimpose the images on the dielectric background to produce an animated, we obtain the following:

<center>
![](images/Fabryperot.gif)
</center>

Compiling
---------

In order to compile your code and link against the Meep library, you must link in several additional libraries that Meep requires. There are three possible ways to do this:

(1) After compiling the Meep package following the instructions elsewhere, place foo.cpp in the `tests/` subdirectory, cd to the same directory, and invoke:

```sh
make foo.dac
```

Run the resulting executable via:

```sh
./foo.dac
```

(2) Use the [pkg-config](https://en.wikipedia.org/wiki/pkg-config) program which is installed by default on most Linux systems:

```sh
g++ `pkg-config --cflags meep` foo.cpp -o foo `pkg-config --libs meep`
```

Naturally, replace `g++` with the name of your C++ compiler if you are not using the GNU compilers.

(3) Compile with g++, this time invoking each library separately:

```sh
g++ -malign-double foo.cpp -o foo -lmeep -lhdf5 -lz -lgsl -lharminv -llapack -lcblas -latlas -lfftw3 -lm
```
