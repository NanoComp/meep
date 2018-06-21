# Mode Decomposition

Meep contains a feature to decompose arbitrary fields into a superposition of the harmonic modes of a given structure using the eigenmode solver [MPB](https://mpb.readthedocs.io). This section provides a description of the analytical as well as implementation details of this feature. Tutorial examples are provided in [Tutorial/Mode Decomposition](Python_Tutorials/Mode_Decomposition/).

[TOC]

## Theoretical Background

The analytical theory for waveguide mode decomposition is described in Chapter 31 ("Modal methods for Maxwell's equations") of [Optical Waveguide Theory](http://www.springer.com/us/book/9780412099502) by Snyder and Love.

Consider a waveguide with propagation axis along the $x$ direction and constant cross section in the transverse direction $\vec\rho=(y,z)$. For a given angular frequency ω we can solve for the eigenmodes of the structure. Thus, arbitrary fields of the form $\mathbf{E}(\mathbf{r},t) = \mathbf{E}(\mathbf{r}) e^{-i\omega t}$ and $\mathbf{H}(\mathbf{r},t) = \mathbf{H}(\mathbf{r}) e^{-i\omega t}$ can be decomposed into a basis of these eigenmodes:

$$
   \mathbf{E}(\mathbf{r}) = 
   \mathbf{E}(x,\vec{\rho}) =
   \sum_{n} \left\{   \alpha^+_n \mathbf E^+_n(\vec \rho)e^{+i\beta_n x}
                    + \alpha^-_n \mathbf E^-_n(\vec \rho)e^{-i\beta_n x}
            \right\}
$$
$$
   \mathbf{H}(\mathbf{r}) = 
   \mathbf{H}(x,\vec{\rho}) =
   \sum_{n} \left\{   \alpha^+_n \mathbf H^+_n(\vec \rho)e^{+i\beta_n x}
                    + \alpha^-_n \mathbf H^-_n(\vec \rho)e^{-i\beta_n x}
            \right\}
$$

β$_n$ are the propagation wavevectors and α$^{\pm}_n$ are the basis coefficients. Mode decomposition involves solving for these unknown quantities. The following steps are involved in the computation:

1.  In Meep, compute the Fourier-transformed fields $\mathbf{E}(\mathbf{r})$ and $\mathbf{H}(\mathbf{r})$ on a surface that is transverse to the waveguide and stored in a `dft_flux` object.

2.  In MPB, compute the eigenmodes $\mathbf{E}^\pm_n$ and $\mathbf{H}^\pm_n$ as well as the propagation wavevectors β$_n$ for the same cross-sectional structure.

3.  Compute the coefficients α$_n^\pm$ for any number of eigenmodes n=1,2,....

This is all done automatically in Meep using the `get_eigenmode_coefficients` routine.

## Function Description

The mode-decomposition feature is available via the `meep::fields::get_eigenmode_coefficients` function callable from Python or C++. This function makes use of several lower-level functions which are described in more detail below. The C++ prototype for this routine is:

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
                                        void *user_kpoint_data)
```
The following are the parameters:

+ `flux` is a `dft_flux` object containing the frequency-domain fields on a cross-sectional slice perpendicular to the waveguide axis

+ `eig_vol` is the volume passed to [MPB](https://mpb.readthedocs.io) for the eigenmode calculation; in most cases this will simply be the volume over which the frequency-domain fields are tabulated (i.e. `flux.where`).

+ `bands` is an array of integers corresponding to the mode indices

+ `num_bands` is the length of the `bands` array

+ `parity` is the parity (= polarization in 2d) of the mode to calculate, assuming the structure has $z$ and/or $y$ mirror symmetry in the source region. If the structure has both $y$ and $z$ mirror symmetry, you can combine more than one of these, e.g. `EVEN_Z + ODD_Y`. This is especially useful in 2d simulations to restrict yourself to a desired polarization

+ `eig_resolution` is the spatial resolution to use in MPB for the eigenmode calculations

+ `eigensolver_tol` is the tolerance to use in the MPB eigensolver. MPB terminates when the eigenvalues stop changing to less than this fractional tolerance

+ `coeffs` is a user-allocated array of type `std::complex<double>` (shortened to `vector<cdouble>`) of length `2*num_freqs*num_bands` where `num_freqs` is the number of frequencies stored in the `flux` object (equivalent to `flux->Nfreq`) and `num_bands` is the length of the `bands` input array. The expansion coefficients for the mode with frequency `nf` and band index `nb`  are stored sequentially as α$^+$, α$^-$ starting at slot `2*nb*num_freqs+nf` of this array

+ `vgrp` is an optional user-allocated `double`-valued array of length `num_freqs*num_bands.` On return, `vgrp[nb*num_freqs + nf]` is the group velocity of the mode with frequency `nf` and band index `nb.` If you do not need this information, simply pass `NULL` for this parameter.

+ `user_kpoint_func` is an optional function you supply to provide an initial guess of the wavevector of a mode with given frequency and band index having the following prototype:

```c++
vec (*kpoint_func)(double freq, int mode, void *user_data);
```

+ `user_kpoint_data` is the user data passed to the `user_kpoint_func`


```c++
 int num_bands = bands.size();
 int num_freqs = Flux->Nfreq;

 std::vector<cdouble> coeffs(2 * num_bands * num_freqs);
 f.get_eigenmode_coefficients(...);

 for(int nb=0; nb<num_bands; nb++)
  for(int nf=0; nf<num_freqs++; nf++)
   { 
     // get coefficients of forward- and backward-traveling
     // waves in eigenmode bands[nb] at frequency #nf
     cdouble AlphaPlus  = coeffs[2*nb*num_freqs + nf + 0];
     cdouble AlphaMinus = coeffs[2*nb*num_freqs + nf + 1];
     ...
```

## Normalization

The α coefficients computed by `get_eigenmode_coefficients` are normalized to ensure that their squared magnitude equals the power carried by the corresponding eigenmode:

$$|\alpha_n^\pm|^2 = P_n^\pm$$

where P$_n^\pm$ is the power carried by $\pm$-traveling eigenmode n. This is discussed in more detail [below](#HowItWorks).

## Related Computational Routines

Besides `get_eigenmode_coefficients,` there are a few computational routines in `libmeep` that you may find useful for problems like those considered above.

### Computing MPB Eigenmodes
````
  void *fields::get_eigenmode(double &omega,
                              direction d, const volume &where,
                              const volume &eig_vol,
                              int band_num,
                              const vec &kpoint, bool match_frequency,
                              int parity,
                              double resolution,
                              double eigensolver_tol);
````

Calls MPB to compute the `band_num`th eigenmode at frequency `omega` for the portion of your geometry lying in `where` which is typically a cross-sectional slice of a waveguide. `kpoint` is an initial starting guess for what the propagation vector of the waveguide mode will be. This is implemented in [mpb.cpp](https://github.com/stevengj/meep/blob/master/src/mpb.cpp#L187-L485).

### Working with MPB Eigenmodes

The return value of `get_eigenmode` is an opaque pointer to a data structure storing information about the computed eigenmode, which may be passed to the following routines:

````
// get a single component of the eigenmode field at a given point in space
std::complex<double> eigenmode_amplitude(const vec &p, void *vedata, component c);

// get the group velocity of the eigenmode 
double get_group_velocity(void *vedata);

// free all memory associated with the eigenmode
void destroy_eigenmode_data(void *vedata);
````

These are implemented in [mpb.cpp](https://github.com/stevengj/meep/blob/master/src/mpb.cpp).

### Exporting Frequency-Domain Fields

````
  void output_dft(dft_flux flux, const char *HDF5FileName);

  void output_mode_fields(void *mode_data, dft_flux flux, const char *HDF5FileName);
````

`output_dft` exports the components of the frequency-domain fields stored in `flux` to an HDF5 file with the given filename In general, `flux` will store data for fields at multiple frequencies.

`output_mode_fields` is similar, but instead exports the components of the eigenmode described by `mode_data` which should be the return value of a call to `get_eigenmode`.

These are implemented in [dft.cpp](https://github.com/stevengj/meep/blob/master/src/dft.cpp#L1077-L1088).

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

These are implemented in [dft.cpp](https://github.com/stevengj/meep/blob/master/src/dft.cpp).

<a name="HowItWorks"></a>
## How Mode Decomposition Works

The theoretical basis of the mode-decomposition algorithm is the orthogonality relation satisfied by the normal modes:

$$ \left\langle \mathbf{E}_m^{\sigma} \right|
   \left.       \mathbf{H}^\tau_n     \right\rangle
   =C_{m}\delta_{mn}\delta_{\sigma\tau} 
   \qquad \{\sigma,\tau\}\in\{+,-\}
$$

where the inner product involves an integration over transverse coordinates:

$$ \left\langle \mathbf{f} \right| \left. \mathbf{g} \right\rangle 
   \equiv
   \int_{S}
    \Big[ \mathbf{f}^*(\vec \rho) \times \mathbf{g}(\vec \rho)\Big]
    \cdot \hat{\mathbf{n}} \, dA
  \tag{5}
$$

where $S$ is any surface transverse to the direction of propagation and $\hat{\mathbf{n}}$ is the unit normal vector to $S$ (i.e. $\hat{\mathbf{z}}$ in the case considered above). The normalization constant $C_{m}$ is a matter of convention, but in MPB it is taken to be the group velocity of the mode, $v_m$, times the area $A_S$ of the cross-sectional surface $S$: $$C_m = v_m A_S$$.

Now consider a Meep calculation in which we have accumulated frequency-domain fields $\mathbf E^{\text{meep}}$ and $\mathbf H^{\text{meep}}$ on a `dft_flux` object located on a cross-sectional surface $S$. Invoking the eigenmode expansion and choosing the origin of the $x$ axis to be the position of the cross-sectional plane, the tangential components of the frequency-domain Meep fields take the form:

$$ \mathbf E^{\text{meep}}_\parallel = \sum_{n} (a_n^+ + a_n^-)\mathbf{E}_{n\parallel}^+ $$

$$ \mathbf H^{\text{meep}}_\parallel = \sum_{n} (a_n^+ - a_n^-)\mathbf{H}_{n\parallel}^+ $$

We have used the well-known relations between the tangential components of the forward- and backward-traveling field modes: 

$$ \mathbf{E}^+_{n\parallel} =+\mathbf{E}^-_{n\parallel}
   \qquad
   \mathbf{H}^+_{n\parallel} =-\mathbf{H}^-_{n\parallel}
$$

Taking the inner product of both equations with the $\mathbf{H}$ and $\mathbf{E}$ fields of each eigenmode, we find

$$ \left\langle \mathbf{H}_m \right|\left. \mathbf{E}^{\text{meep}} \right\rangle =+(a_m^+ + a_m^-) v_m A_S $$

$$ \left\langle \mathbf{E}_m \right|\left. \mathbf{H}^{\text{meep}} \right\rangle =-(a_m^+ - a_m^+) v_m A_S $$

Thus, by evaluating the integrals on the left-hand side of these equations &mdash; numerically, using the MPB-computed eigenmode fields $\{\mathbf{E}_m, \mathbf{H}_m\}$ and the Meep-computed fields $\{\mathbf{E}^{\text{meep}}, \mathbf{H}^{\text{meep}}\}$ as tabulated on the computational grid &mdash; and combining the results appropriately, we can extract the coefficients $\alpha^\pm_m$. This calculation is carried out by the routine `meep::fields::get_mode_flux_overlap`. Although simple in principle, the implementation is complicated by the fact that, in multi-processor calculations, the Meep fields needed to evaluate the integrals are generally not all present on any one processor, but are instead distributed over multiple processors, requiring some interprocess communication to evaluate the full integral.

The Poynting flux carried by the Meep fields may be expressed in the form:

$$ S_x = \frac{1}{2}\text{Re }
         \left\langle \mathbf{E}^{\text{meep}}\right|
         \left.       \mathbf{H}^{\text{meep}}\right\rangle
       = \frac{1}{2}\sum_n \left\{ |a_n^+|^2 - |a_n^-|^2) \right\} v_n A_S
$$

Thus, the power carried by a given forward- or backward-traveling eigenmode is given by:

$$ \textit{power} = \frac{|a_n^\pm|^2 v_n A_S}{2S_x} $$

or alternatively,

$$ \textit{power} = |\alpha_n^\pm|^2 $$

where we have defined $\alpha_n^\pm \equiv \sqrt{v_n A_S}/\left({2S_x}\right)a_n^\pm$. These coefficients $\alpha_n^\pm$ are the quantities computed by `get_eigenmode_coefficients.`


