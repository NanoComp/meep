# Mode Decomposition

Meep contains a feature to decompose arbitrary fields into a superposition of the harmonic modes of a given structure via its integration with the eigenmode solver [MPB](https://mpb.readthedocs.io). Only dielectric structures with lossless, wavelength-independent, anisotropic $\varepsilon$ are supported by MPB. If dispersive materials are defined in Meep, only the real part of $\varepsilon$ is used. This section provides an overview of the theory and implementation of this feature. For examples, see [Tutorial/Mode Decomposition](Python_Tutorials/Mode_Decomposition.md).

[TOC]

## Theoretical Background

The theory underlying mode decomposition is described in Chapter 31 ("Modal methods for Maxwell's equations") of [Optical Waveguide Theory](http://www.springer.com/us/book/9780412099502) by Snyder and Love.

Consider a waveguide with propagation axis along the $x$ direction and constant cross section in the transverse directions $(y,z)$.   Let $\psi = (E_y,E_z,H_y,H_z)$ denote the transverse components of the electric and magnetic fields.  For a given angular frequency $\omega$ we can solve for the eigenmodes of the structure: solutions of the form $\psi^\pm_n(y,z) e^{\pm i \beta_n x - i\omega t}$, where $\beta_n$ are the propagation constants of the right ($+$) and left ($-$) traveling modes.   (There are also evanescent modes with complex $\beta_n$ values, but we will focus mainly here on the propagating modes with real $\beta_n$.)

Any *arbitrary* fields $\psi$, Fourier-transformed to a particular $\omega$, can be expressed in the basis of these eigenmodes:

$$
\psi(x,y,z) = \begin{pmatrix} E_y \\ E_z \\ H_y \\ H_z \end{pmatrix} =
 \sum_{n} \left\{   \alpha^+_n \mathbf \psi^+_n(\vec \rho)e^{+i\beta_n x}
                    + \alpha^-_n \mathbf \psi^-_n(\vec \rho)e^{-i\beta_n x}
            \right\}
$$

$\alpha^{\pm}_n$ are the expansion coefficients (the amplitudes of each mode present in the fields). Mode decomposition involves solving for these amplitudes. The following steps are involved in the computation:

1.  In Meep, compute the Fourier-transformed transverse fields $\psi$ on a surface that is transverse to the waveguide, stored in a `dft_flux` object.

2.  In MPB, compute the eigenmodes $\psi^\pm_n$ as well as the propagation constants $\beta_n$ for the same cross-sectional structure.

3.  Compute the coefficients $\alpha_n^\pm$ for any number of eigenmodes $n=1,2,...$.

This is all done automatically in Meep using the `get_eigenmode_coefficients` routine.
Meep normalizes the modes $\psi^\pm_n$ to unit power $\Re \int \mathbf{E}^*\times\mathbf{H} = 1$, so that $|\alpha^{\pm}_n|^2$ is equal to the **power**
(Poynting flux) carried by that mode in the Fourier-transformed field $\psi$.

## Function Description

The mode-decomposition feature is available via the `meep::fields::get_eigenmode_coefficients` function callable from Python or C++. This function makes use of several lower-level functions which are described in more detail below. The C++ header for this function is:

```c++
void fields::get_eigenmode_coefficients(dft_flux flux,
                                        const volume &eig_vol,
                                        int *bands,
                                        int num_bands,
                                        int parity,
                                        double eig_resolution,
                                        double eigensolver_tol,
                                        std::complex<double> *coeffs,
                                        double *vgrp,
                                        kpoint_func user_kpoint_func,
                                        void *user_kpoint_data,
                                        vec *kpoint_list,
                                        vec *kdom_list,
                                        double *cscale,
                                        direction d,
                                        diffractedplanewave *dp)
```
The following are the parameters:

+ `flux` is a `dft_flux` object containing the frequency-domain fields on a cross-sectional slice perpendicular to the waveguide axis.

+ `eig_vol` is the `volume` passed to [MPB](https://mpb.readthedocs.io) for the eigenmode calculation. In most cases this will simply be the volume over which the frequency-domain fields are tabulated (i.e. `flux.where`).

+ `bands` is an array of integers corresponding to the mode indices (equivalent to $n$ in the two formulas above).

+ `num_bands` is the length of the `bands` array.

+ `parity` is the parity of the mode to calculate, assuming the structure has $z$ and/or $y$ mirror symmetry in the source region. If the structure has both $y$ and $z$ mirror symmetry, you can combine more than one of these, e.g. `ODD_Z+EVEN_Y`. This is especially useful in 2d simulations to restrict yourself to a desired polarization.

+ `eig_resolution` is the spatial resolution to use in MPB for the eigenmode calculations.

+ `eigensolver_tol` is the tolerance to use in the MPB eigensolver. MPB terminates when the eigenvalues stop changing by less than this fractional tolerance.

+ `coeffs` is a user-allocated array of type `std::complex<double>` (shortened hereafter to `cdouble`) of length `2*num_freqs*num_bands` where `num_freqs` is the number of frequencies stored in the `flux` object (equivalent to `flux->Nfreq`) and `num_bands` is the length of the `bands` input array. The expansion coefficients for the mode with frequency `nf` and band index `nb`  are stored sequentially as $\alpha^+$, $\alpha^-$ starting at slot `2*nb*num_freqs+nf` of this array

+ `vgrp` is an optional user-allocated `double` array of length `num_freqs*num_bands`. On return, `vgrp[nb*num_freqs + nf]` is the group velocity of the mode with frequency `nf` and band index `nb.` If you do not need this information, simply pass `NULL` for this parameter.

+ `user_kpoint_func` is an optional function you supply to provide an initial guess of the wavevector of a mode with given frequency and band index having the following prototype:

```c++
vec (*kpoint_func)(double freq, int mode, void *user_data);
```

+ `user_kpoint_data` is the user data passed to the `user_kpoint_func`.

+ `kpoint_list` is a user allocated array of `meep::vec` objects of length (`num_bands*num_freqs`). If non-null, this array is filled in with the wavevectors of the eigenmode for each band from 1 to `num_bands*num_freqs`.

+ `kdom_list` is a user allocated array of `meep::vec` objects of length (`num_bands*num_freqs`). If non-null, this array is filled in with the wavevectors of the dominant planewave in the Fourier series expansion for each band from 1 to `num_bands*num_freqs`. `kdom_list[nb*num_freqs + nf]` is the dominant planewave of the mode with frequency `nf` and band index `nb` which defaults to `NULL`.  This is especially useful for interpreting the modes computed in a uniform medium, because those modes are exactly planewaves proportional to $\exp(2\pi i \vec{k}_{dom}\cdot \vec{x})$ where $\vec{k}_{dom}$ (`kdom`) is the wavevector.

+ `cscale` is a user allocated array of `double` objects of length (`num_bands*num_freqs`). If non-null, this array is filled in with scalar coefficients used for adjoint calculations.

+ `d` is the direction normal to the monitor plane.

+ `dp` is a user allocated `diffractedplanewave` used to specify a diffracted planewave in homogeneous media.

The following is a demonstration of typical usage:

```c++
 int num_bands = bands.size();
 int num_freqs = Flux->Nfreq;

 std::vector<cdouble> coeffs(2*num_bands*num_freqs);
 f.get_eigenmode_coefficients(...);

 for(int nb=0; nb<num_bands; nb++)
  for(int nf=0; nf<num_freqs++; nf++)
   {
     // get coefficients of forward- and backward-traveling
     // waves in eigenmode bands[nb] at frequency nf
     cdouble AlphaPlus = coeffs[2*nb*num_freqs+nf+0];
     cdouble AlphaMinus = coeffs[2*nb*num_freqs+nf+1];
     ...
```

## Normalization

The $\alpha$ coefficients computed by `get_eigenmode_coefficients` are normalized such that their squared magnitude equals the power carried by the corresponding eigenmode:

$$|\alpha_n^\pm|^2 = P_n^\pm$$

where $P_n^\pm$ is the power carried by the traveling eigenmode $n$ in the forward (+) or backward (-) direction. This is discussed in more detail below.

## Related Functions

Besides `get_eigenmode_coefficients`, there are a few computational routines in `libmeep` that you may find useful for problems like those considered above.

### Computing MPB Eigenmodes
````
  void *fields::get_eigenmode(double frequency,
                              direction d,
                              const volume where,
                              const volume eig_vol,
                              int band_num,
                              const vec kpoint,
                              bool match_frequency,
                              int parity,
                              double resolution,
                              double eigensolver_tol,
                              double *kdom,
                              void **user_mdata,
                              diffractedplanewave *dp);
````

Calls MPB to compute the `band_num`th eigenmode at frequency `frequency` for the portion of your geometry lying in `where` which is typically a cross-sectional slice of a waveguide. `kpoint` is an initial starting guess for what the propagation vector of the waveguide mode will be. `kdom`, if non-NULL and length 3, is filled in with the dominant planewave for the current band (see above). This is implemented in [src/mpb.cpp](https://github.com/NanoComp/meep/blob/master/src/mpb.cpp).

### Working with MPB Eigenmodes

The return value of `get_eigenmode` is an [opaque pointer](https://en.wikipedia.org/wiki/Opaque_pointer) to a data structure storing information about the computed eigenmode, which may be passed to the following routines:

````
// get a single component of the eigenmode field at a given point in space
std::complex<double> eigenmode_amplitude(const vec &p, void *vedata, component c);

// get the group velocity of the eigenmode
double get_group_velocity(void *vedata);

// free all memory associated with the eigenmode
void destroy_eigenmode_data(void *vedata);
````

These functions are implemented in [src/mpb.cpp](https://github.com/NanoComp/meep/blob/master/src/mpb.cpp).

### Exporting Frequency-Domain Fields

````
  void output_dft(dft_flux flux, const char *HDF5FileName);

````

`output_dft` exports the components of the frequency-domain fields stored in `flux` to an HDF5 file with the given filename In general, `flux` will store data for fields at multiple frequencies.

This function is implemented in [src/dft.cpp](https://github.com/NanoComp/meep/blob/master/src/dft.cpp#L1184-L1196).

### Computing Overlap Integrals
````
  std::complex<double> get_mode_flux_overlap(void *mode_data,
                                             dft_flux *flux,
                                             int num_freq,
                                             std::complex<double>overlap[2]);

  std::complex<double> get_mode_mode_overlap(void *mode1_data,
                                             void *mode2_data,
                                             dft_flux *flux,
                                             std::complex<double>overlap[2]);
````

`get_mode_flux_overlap` computes the overlap integral between the eigenmode described by `mode_data` and the fields stored in `flux` for the `num_freq`th stored frequency, where `num_freq` ranges from 0 to `flux->Nfreq-1`. `mode_data` should be the return value of a previous call to `get_eigenmode.`

`get_mode_mode_overlap` is similar, but computes the overlap integral between two eigenmodes. `mode1_data` and `mode2_data` may be identical, in which case you get the inner product of the mode with itself. This should equal the group velocity of the mode based on the MPB's normalization convention.

These functions are implemented in [src/dft.cpp](https://github.com/NanoComp/meep/blob/master/src/dft.cpp).

## How Mode Decomposition Works

The theoretical basis of the mode-decomposition algorithm is an orthogonality relation satisfied by the eigenmodes:

$$
\left\langle \psi_m^{\sigma}, \psi_n^{\tau} \right\rangle
=C_{m}\delta_{mn}\delta_{\sigma\tau}
   \qquad \{\sigma,\tau\}\in\{+,-\}
$$

where the (indefinite) inner product between two vectors of the transverse field components, e.g. $\psi=(E_y,E_z,H_y,H_z)$ and $\psi'=(E_y',E_z',H_y',H_z')$ for propagation in the $x$ direction, involves an integration over transverse coordinates:

$$ \left\langle \psi , \psi' \right\rangle
   \equiv
   \int_{S}
    \Big[ \mathbf{E}^*(\vec \rho) \times \mathbf{H}'(\vec \rho) + \mathbf{E}'(\vec \rho) \times \mathbf{H}^*(\vec \rho)\Big]
    \cdot \hat{\mathbf{n}} \, dA
  \tag{5}
$$

where $S$ is a cross-section transverse to the direction of propagation and $\hat{\mathbf{n}}$ is the unit normal vector to $S$ (i.e. $\hat{\mathbf{x}}$ in the case considered above).  The normalization constant $C_{m}$ is a matter of convention, but for MPB-computed modes (which normalizes proportional to unit energy) it is effectively the group velocity of the mode, $v_m$, times the area $A_S$ of the cross-sectional surface $S$: $$C_m = 2 v_m A_S$$   However, we can divide the MPB modes by $\sqrt{\int_S \Re[\mathbf{E}^* \times \mathbf{H}]\cdot \hat{\mathbf{n}}}$ to renormalize them to $C_m = 2$, which is equivalent to normalizing the MPB modes to unit power.

(There is some subtlety about the use of complex conjugation for orthogonality in the inner product above that arises for evanescent modes, which we don't consider here because Meep only computes coefficients of propagating modes.  The above definition has the nice property that $\langle \psi, \psi \rangle = 2 \int_S \Re[\mathbf{E}^* \times \mathbf{H}]\cdot \hat{\mathbf{n}} = 2P$ where $P$ is the Poynting flux.)

Now consider a Meep calculation in which we have accumulated transverse frequency-domain fields $\psi^{\text{meep}} = (\mathbf E^{\text{meep}}_\parallel, \mathbf H^{\text{meep}}_\parallel)$ on a `dft_flux` object located on a cross-sectional surface $S$. Invoking the eigenmode expansion and choosing the origin of the $x$ axis to be the position of the cross-sectional plane means that
Meep fields take the form:

$$
\psi^{\text{meep}} =
 \sum_{n} \left\{   \alpha^+_n \mathbf \psi^+_n(\vec \rho)e^{+i\beta_n x}
                    + \alpha^-_n \mathbf \psi^-_n(\vec \rho)e^{-i\beta_n x}
            \right\}
$$

Taking the inner product of this equation with the $\psi_m^\pm$ modes computed
from MPB (renormalized to unit power), we obtain the mode coefficients:

$$
\alpha^\pm_n = \left\langle \psi_m^\pm, \psi^{\text{meep}} \right\rangle \, .
$$

Some simplifications arise from the fact that, in a constant cross-section waveguide, the tangential components of the forward- and backward-traveling propagating modes are related by a mirror reflection in $x$:
$\mathbf{E}^+_{n\parallel} =+\mathbf{E}^-_{n\parallel}$
and $\mathbf{H}^+_{n\parallel} =-\mathbf{H}^-_{n\parallel}$.  This means that we only need to compute the MPB modes once for the forward (+) direction.  Also, the $\langle \psi_m^\pm, \psi^{\text{meep}} \rangle$ inner
product involves four integrals (of terms like $E_y H_z$ etc.) that are computed
once each and then combined with the correct signs to obtain $\alpha^\pm_n$.
These integrals are computed numerically by a trapezoidal-type rule on Meep's
Yee grid, interpolating the MPB-computed eigenmode fields as needed. This calculation is carried out by the routine `fields::get_mode_flux_overlap`. Although simple in principle, the implementation is complicated by the fact that, in multi-processor calculations, the Meep fields needed to evaluate the integrals are generally not all present on any one processor, but are instead distributed over multiple processors, requiring some interprocess communication to evaluate the full integral.

As mentioned above, the Poynting flux $P$ carried by the Meep fields may be expressed in the form:

$$
P = \frac{1}{2} \langle \psi, \psi \rangle = \sum_n \left( |a_n^+|^2 - |a_n^-|^2) \right)
$$

where the right-hand-side is obtained by substituting the modal expansion and
using the mode orthogonality and the normalization $C_n = 2$ chosen above.
Thus, the power carried by a given forward- or backward-traveling eigenmode is given by:

$$ |\textit{power}_n^\pm| = |\alpha_n^\pm|^2 $$

The $\alpha_n^\pm$ are the eigenmode coefficients returned by `get_eigenmode_coefficients`.
