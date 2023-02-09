from typing import Union, Tuple
import functools
import math
import numbers
import operator
import warnings
from collections import namedtuple
from copy import deepcopy
from numbers import Number
import numpy as np
import meep as mp

FreqRange = namedtuple("FreqRange", ["min", "max"])


def check_nonnegative(prop, val):
    if val >= 0:
        return val
    else:
        raise ValueError(f"{prop} cannot be negative. Got {val}")


def init_do_averaging(mat_func):
    if not hasattr(mat_func, "do_averaging"):
        mat_func.do_averaging = False


class Vector3:
    """
    Properties:

    **`x`, `y`, `z` [`float` or `complex`]** — The `x`, `y`, and `z` components of the
    vector. Generally, functions that take a `Vector3` as an argument will accept an
    iterable (e.g., a tuple or list) and automatically convert to a `Vector3`.
    """

    def __init__(self, x: float = 0.0, y: float = 0.0, z: float = 0.0):
        """
        Create a new `Vector3` with the given components. All three components default to
        zero. This can also be represented simply as `(x,y,z)` or `[x,y,z]`.
        """
        self.x = float(x) if type(x) is int else x
        self.y = float(y) if type(y) is int else y
        self.z = float(z) if type(z) is int else z

    def __eq__(self, other):
        """
        Returns whether or not the two vectors are numerically equal. Beware of using this
        function after operations that may have some error due to the finite precision of
        floating-point numbers; use `close` instead.

        ```python
        v1 == v2
        ```
        """
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __ne__(self, other):
        """
        Returns whether or not the two vectors are numerically unequal. Beware of using
        this function after operations that may have some error due to the finite
        precision of floating-point numbers; use `close` instead.

        ```python
        v1 != v2
        ```
        """
        return not self == other

    def __add__(self, other):
        """
        Return the sum of the two vectors.

        ```python
        v3 = v1 + v2
        ```
        """
        if isinstance(other, GeometricObject):
            return NotImplemented

        x = self.x + other.x
        y = self.y + other.y
        z = self.z + other.z

        return Vector3(x, y, z)

    def __sub__(self, other):
        """
        Return the difference of the two vectors.

        ```python
        v3 = v1 - v2
        ```
        """
        x = self.x - other.x
        y = self.y - other.y
        z = self.z - other.z

        return Vector3(x, y, z)

    def __mul__(self, other):
        """
        If `other` is a `Vector3`, returns the dot product of `v1` and `other`. If `other`
        is a number, then `v1` is scaled by the number.

        ```python
        c = v1 * other
        ```
        """
        if type(other) is Vector3:
            return self.dot(other)
        elif isinstance(other, Number):
            return self.scale(other)
        else:
            raise TypeError(f"No operation known for 'Vector3 * {type(other)}'")

    def __truediv__(self, other):
        if type(other) is Vector3:
            return Vector3(self.x / other.x, self.y / other.y, self.z / other.z)
        elif isinstance(other, Number):
            return Vector3(self.x / other, self.y / other, self.z / other)
        else:
            raise TypeError(f"No operation known for 'Vector3 / {type(other)}'")

    def __rmul__(self, other):
        """
        If `other` is a `Vector3`, returns the dot product of `v1` and `other`. If `other`
        is a number, then `v1` is scaled by the number.

        ```python
        c = other * v1
        ```
        """
        if isinstance(other, Number):
            return self.scale(other)
        else:
            raise TypeError(f"No operation known for '{type(other)} * Vector3'")

    def __getitem__(self, i):
        if i == 0:
            return self.x
        elif i == 1:
            return self.y
        elif i == 2:
            return self.z
        else:
            raise IndexError(f"No value at index {i}")

    def __repr__(self):
        return f"Vector3<{self.x}, {self.y}, {self.z}>"

    def __array__(self):
        return np.array([self.x, self.y, self.z])

    def conj(self):
        return Vector3(self.x.conjugate(), self.y.conjugate(), self.z.conjugate())

    def scale(self, s):
        x = self.x * s
        y = self.y * s
        z = self.z * s

        return Vector3(x, y, z)

    def dot(self, v):
        """
        Returns the dot product of *`self`* and *`v`*.

        ```python
        v3 = v1.dot(v2)
        ```
        """
        return self.x * v.x + self.y * v.y + self.z * v.z

    def cdot(self, v):
        """Returns the conjugated dot product: `conj(self)` dot `v`."""
        return self.conj().dot(v)

    def cross(self, v):
        """
        Return the cross product of `self` and `v`.

        ```python
        v3 = v1.cross(v2)
        ```
        """
        x = self.y * v.z - self.z * v.y
        y = self.z * v.x - self.x * v.z
        z = self.x * v.y - self.y * v.x

        return Vector3(x, y, z)

    def norm(self):
        """
        Returns the length `math.sqrt(abs(self.dot(self)))` of the given vector.

        ```python
        v2 = v1.norm()
        ```
        """
        return math.sqrt(abs(self.cdot(self).real))

    def unit(self):
        """
        Returns a unit vector in the direction of the vector.

        ```python
        v2 = v1.unit()
        ```
        """
        return self.scale(1 / self.norm())

    def close(self, v, tol=1.0e-7):
        """
        Returns whether or not the corresponding components of the `self` and `v` vectors
        are within `tol` of each other. Defaults to 1e-7.

        ```python
        v1.close(v2, [tol])
        ```
        """
        return (
            abs(self.x - v.x) <= tol
            and abs(self.y - v.y) <= tol
            and abs(self.z - v.z) <= tol
        )

    def rotate(self, axis, theta):
        """
        Returns the vector rotated by an angle *`theta`* (in radians) in the right-hand
        direction around the *`axis`* vector (whose length is ignored). You may find the
        python functions `math.degrees` and `math.radians` useful to convert angles
        between degrees and radians.

        ```python
        v2 = v1.rotate(axis, theta)
        ```
        """
        u = axis.unit()
        vpar = u.scale(u.dot(self))
        vcross = u.cross(self)
        vperp = self - vpar
        return vpar + (vperp.scale(math.cos(theta)) + vcross.scale(math.sin(theta)))

    # rotate vectors in lattice/reciprocal coords (note that the axis
    # is also given in the corresponding basis):

    def rotate_lattice(self, axis, theta, lat):
        a = lattice_to_cartesian(axis, lat)
        v = lattice_to_cartesian(self, lat)
        return cartesian_to_lattice(v.rotate(a, theta), lat)

    def rotate_reciprocal(self, axis, theta, lat):
        a = reciprocal_to_cartesian(axis, lat)
        v = reciprocal_to_cartesian(self, lat)
        return cartesian_to_reciprocal(v.rotate(a, theta), lat)


class Medium:
    """
    This class is used to specify the materials that geometric objects are made of. It
    represents an electromagnetic medium which is possibly nonlinear and/or dispersive.
    See also [Materials](Materials.md). To model a perfectly-conducting metal, use the
    predefined `metal` object, above. To model imperfect conductors, use a dispersive
    dielectric material. See also the [Predefined Variables](#predefined-variables):
    `metal`, `perfect_electric_conductor`, and `perfect_magnetic_conductor`.

    **Material Function**

    Any function that accepts a `Medium` instance can also accept a user-defined Python
    function. This allows you to specify the material as an arbitrary function of
    position. The function must have one argument, the position `Vector3`, and return the
    material at that point, which should be a Python `Medium` instance. This is
    accomplished by passing a function to the `material_function` keyword argument in the
    `Simulation` constructor, or the `material` keyword argument in any `GeometricObject`
    constructor. For an example, see [Subpixel Smoothing/Enabling Averaging for Material
    Function](Subpixel_Smoothing.md#enabling-averaging-for-material-function).

    Instead of the `material` or `material_function` arguments, you can also use the
    `epsilon_func` keyword argument to `Simulation` and `GeometricObject`, which takes a
    function of position that returns the dielectric constant at that point.

    **Important:** If your material function returns nonlinear, dispersive (Lorentzian or
    conducting), or magnetic materials, you should also include a list of these materials
    in the `extra_materials` input variable (above) to let Meep know that it needs to
    support these material types in your simulation. For dispersive materials, you need to
    include a material with the *same* values of $\\gamma_n$ and $\\omega_n$, so
    you can only have a finite number of these, whereas $\\sigma_n$ can vary
    continuously and a matching $\\sigma_n$ need not be specified in
    `extra_materials`. For nonlinear or conductivity materials, your `extra_materials`
    list need not match the actual values of $\\sigma$ or $\\chi$ returned by your material function,
    which can vary continuously.

    **Complex $\\varepsilon$ and $\\mu$**: you cannot specify a
    frequency-independent complex $\\varepsilon$ or $\\mu$ in Meep where
    the imaginary part is a frequency-independent loss but there is an
    alternative.  That is because there are only two important
    physical situations. First, if you only care about the loss in a
    narrow bandwidth around some frequency, you can set the loss at
    that frequency via the
    [conductivity](Materials.md#conductivity-and-complex).  Second, if
    you care about a broad bandwidth, then all physical materials have
    a frequency-dependent complex $\\varepsilon$ and/or $\\mu$, and you
    need to specify that frequency dependence by fitting to Lorentzian
    and/or Drude resonances via the `LorentzianSusceptibility` or
    `DrudeSusceptibility` classes below.

    Dispersive dielectric and magnetic materials, above, are specified via a list of
    objects that are subclasses of type `Susceptibility`.
    """

    def __init__(
        self,
        epsilon_diag=Vector3(1, 1, 1),
        epsilon_offdiag=Vector3(),
        mu_diag=Vector3(1, 1, 1),
        mu_offdiag=Vector3(),
        E_susceptibilities=None,
        H_susceptibilities=None,
        E_chi2_diag=Vector3(),
        E_chi3_diag=Vector3(),
        H_chi2_diag=Vector3(),
        H_chi3_diag=Vector3(),
        D_conductivity_diag=Vector3(),
        D_conductivity_offdiag=Vector3(),
        B_conductivity_diag=Vector3(),
        B_conductivity_offdiag=Vector3(),
        epsilon=None,
        index=None,
        mu=None,
        chi2=None,
        chi3=None,
        D_conductivity=None,
        B_conductivity=None,
        E_chi2=None,
        E_chi3=None,
        H_chi2=None,
        H_chi3=None,
        valid_freq_range=FreqRange(min=-mp.inf, max=mp.inf),
    ):
        """
        Creates a `Medium` object.

        + **`epsilon` [`number`]** The frequency-independent isotropic relative
          permittivity or dielectric constant. Default is 1. You can also use `index=n` as
          a synonym for `epsilon=n*n`; note that this is not really the refractive index
          if you also specify μ, since the true index is $\\sqrt{\\mu\\varepsilon}$. Using
          `epsilon=ep` is actually a synonym for `epsilon_diag=mp.Vector3(ep, ep, ep)`.

        + **`epsilon_diag` and `epsilon_offdiag` [`Vector3`]** — These properties allow
          you to specify ε as an arbitrary real-symmetric tensor by giving the diagonal
          and offdiagonal parts. Specifying `epsilon_diag=Vector3(a, b, c)` and/or
          `epsilon_offdiag=Vector3(u, v, w)` corresponds to a relative permittivity ε
          tensor \\begin{pmatrix} a & u & v \\\\ u & b & w \\\\ v & w & c \\end{pmatrix}
          Default is the identity matrix ($a = b = c = 1$ and $u = v = w = 0$).

        + **`mu` [`number`]** — The frequency-independent isotropic relative permeability
          μ. Default is 1. Using `mu=pm` is actually a synonym for `mu_diag=mp.Vector3(pm,
          pm, pm)`.

        + **`mu_diag` and `mu_offdiag` [`Vector3`]** — These properties allow you to
          specify μ as an arbitrary real-symmetric tensor by giving the diagonal and
          offdiagonal parts exactly as for ε above. Default is the identity matrix.

        + **`D_conductivity` [`number`]** — The frequency-independent electric
          conductivity $\\sigma_D$. Default is 0. You can also specify a diagonal
          anisotropic conductivity tensor by using the property `D_conductivity_diag`
          which takes a `Vector3` to give the $\\sigma_D$ tensor diagonal. See also
          [Conductivity](Materials.md#conductivity-and-complex).

        + **`B_conductivity` [`number`]** — The frequency-independent magnetic
          conductivity $\\sigma_B$. Default is 0. You can also specify a diagonal
          anisotropic conductivity tensor by using the property `B_conductivity_diag`
          which takes a `Vector3` to give the $\\sigma_B$ tensor diagonal. See also
          [Conductivity](Materials.md#conductivity-and-complex).

        + **`chi2` [`number`]** — The nonlinear electric
          [Pockels](https://en.wikipedia.org/wiki/Pockels_effect) susceptibility
          $\\chi^{(2)}$ (quadratic nonlinearity). Default is 0. See also [Nonlinearity](Materials.md#nonlinearity).
          This is equivalent to setting `E_chi2`; alternatively, an analogous magnetic
          nonlinearity can be specified using `H_chi2`. These are isotropic nonlinearities,
          but *diagonal* anisotropic polarizations of the form $\\chi_i^{(2)} E_i^2$ can
          be specified with `E_chi2_diag` (which defaults to `[E_chi2,E_chi2,E_chi2]`).

        + **`chi3` [`number`]** — The nonlinear electric
          [Kerr](https://en.wikipedia.org/wiki/Kerr_effect) susceptibility $\\chi^{(3)}$
          (cubic nonlinearity). Default is 0. See also [Nonlinearity](Materials.md#nonlinearity).
          This is equivalent to setting `E_chi3`; alternatively, an analogous magnetic nonlinearity
          can be specified using `H_chi3`. These are isotropic nonlinearities, but *diagonal*
          anisotropic polarizations of the form $\\chi_i^{(3)} |E|^2 E_i$ can be specified with
          `E_chi3_diag` (which defaults to `[E_chi3,E_chi3,E_chi3]`).

        + **`E_susceptibilities` [ list of `Susceptibility` class ]** — List of dispersive
          susceptibilities (see below) added to the dielectric constant ε in order to
          model material dispersion. Defaults to none (empty list). See also [Material
          Dispersion](Materials.md#material-dispersion).

        + **`H_susceptibilities` [ list of `Susceptibility` class ]** — List of dispersive
          susceptibilities (see below) added to the permeability μ in order to model
          material dispersion. Defaults to none (empty list). See also [Material
          Dispersion](Materials.md#material-dispersion).
        """

        if epsilon:
            epsilon_diag = Vector3(epsilon, epsilon, epsilon)
        elif index:
            i2 = index * index
            epsilon_diag = Vector3(i2, i2, i2)

        if mu:
            mu_diag = Vector3(mu, mu, mu)

        if D_conductivity:
            D_conductivity_diag = Vector3(
                D_conductivity, D_conductivity, D_conductivity
            )
        if B_conductivity:
            B_conductivity_diag = Vector3(
                B_conductivity, B_conductivity, B_conductivity
            )

        if E_chi2:
            E_chi2_diag = Vector3(E_chi2, E_chi2, E_chi2)
        if E_chi3:
            E_chi3_diag = Vector3(E_chi3, E_chi3, E_chi3)
        if H_chi2:
            H_chi2_diag = Vector3(H_chi2, H_chi2, H_chi2)
        if H_chi3:
            H_chi3_diag = Vector3(H_chi3, H_chi3, H_chi3)

        self.epsilon_diag = Vector3(*epsilon_diag)
        self.epsilon_offdiag = Vector3(*epsilon_offdiag)
        self.mu_diag = Vector3(*mu_diag)
        self.mu_offdiag = Vector3(*mu_offdiag)
        self.E_susceptibilities = E_susceptibilities or []
        self.H_susceptibilities = H_susceptibilities or []
        self.E_chi2_diag = Vector3(chi2, chi2, chi2) if chi2 else Vector3(*E_chi2_diag)
        self.E_chi3_diag = Vector3(chi3, chi3, chi3) if chi3 else Vector3(*E_chi3_diag)
        self.H_chi2_diag = Vector3(*H_chi2_diag)
        self.H_chi3_diag = Vector3(*H_chi3_diag)
        self.D_conductivity_diag = Vector3(*D_conductivity_diag)
        self.D_conductivity_offdiag = Vector3(*D_conductivity_offdiag)
        self.B_conductivity_diag = Vector3(*B_conductivity_diag)
        self.B_conductivity_offdiag = Vector3(*D_conductivity_offdiag)
        self.valid_freq_range = valid_freq_range

    def __repr__(self):
        return "Medium()"

    def transform(self, m):
        """
        Transforms `epsilon`, `mu`, and `sigma` of any [susceptibilities](#susceptibility)
        by the 3×3 matrix `m`. If `m` is a [rotation
        matrix](https://en.wikipedia.org/wiki/Rotation_matrix), then the principal axes of
        the susceptibilities are rotated by `m`.  More generally, the susceptibilities χ
        are transformed to MχMᵀ/|det M|, which corresponds to [transformation
        optics](http://math.mit.edu/~stevenj/18.369/coordinate-transform.pdf) for an
        arbitrary curvilinear coordinate transformation with Jacobian matrix M. The
        absolute value of the determinant is to prevent inadvertent construction of
        left-handed materials, which are [problematic in nondispersive
        media](FAQ.md#why-does-my-simulation-diverge-if-0).
        """
        eps = Matrix(
            mp.Vector3(
                self.epsilon_diag.x, self.epsilon_offdiag.x, self.epsilon_offdiag.y
            ),
            mp.Vector3(
                self.epsilon_offdiag.x, self.epsilon_diag.y, self.epsilon_offdiag.z
            ),
            mp.Vector3(
                self.epsilon_offdiag.y, self.epsilon_offdiag.z, self.epsilon_diag.z
            ),
        )
        mu = Matrix(
            mp.Vector3(self.mu_diag.x, self.mu_offdiag.x, self.mu_offdiag.y),
            mp.Vector3(self.mu_offdiag.x, self.mu_diag.y, self.mu_offdiag.z),
            mp.Vector3(self.mu_offdiag.y, self.mu_offdiag.z, self.mu_diag.z),
        )

        new_eps = (m * eps * m.transpose()) / abs(m.determinant())
        new_mu = (m * mu * m.transpose()) / abs(m.determinant())
        self.epsilon_diag = mp.Vector3(new_eps.c1.x, new_eps.c2.y, new_eps.c3.z)
        self.epsilon_offdiag = mp.Vector3(new_eps.c2.x, new_eps.c3.x, new_eps.c3.y)
        self.mu_diag = mp.Vector3(new_mu.c1.x, new_mu.c2.y, new_mu.c3.z)
        self.mu_offdiag = mp.Vector3(new_mu.c2.x, new_mu.c3.x, new_mu.c3.y)

        for s in self.E_susceptibilities:
            s.transform(m)

        for s in self.H_susceptibilities:
            s.transform(m)

    def rotate(self, axis, theta):
        T = get_rotation_matrix(axis, theta)
        self.transform(T)

    def epsilon(self, freq):
        """
        Returns the medium's permittivity tensor as a 3x3 Numpy array at the specified
        frequency `freq` which can be either a scalar, list, or Numpy array. In the case
        of a list/array of N frequency points, a Numpy array of size Nx3x3 is returned.
        """
        return self._get_epsmu(
            self.epsilon_diag,
            self.epsilon_offdiag,
            self.E_susceptibilities,
            self.D_conductivity_diag,
            self.D_conductivity_offdiag,
            freq,
        )

    def mu(self, freq):
        """
        Returns the medium's permeability tensor as a 3x3 Numpy array at the specified
        frequency `freq` which can be either a scalar, list, or Numpy array. In the case
        of a list/array of N frequency points, a Numpy array of size Nx3x3 is returned.
        """
        return self._get_epsmu(
            self.mu_diag,
            self.mu_offdiag,
            self.H_susceptibilities,
            self.B_conductivity_diag,
            self.B_conductivity_offdiag,
            freq,
        )

    def _get_epsmu(
        self,
        diag,
        offdiag,
        susceptibilities,
        conductivity_diag,
        conductivity_offdiag,
        freq,
    ):
        # Clean the input
        if np.isscalar(freq):
            freqs = np.array(freq)[np.newaxis, np.newaxis, np.newaxis]
        else:
            freqs = np.squeeze(freq)
            freqs = freqs[:, np.newaxis, np.newaxis]

        # Check for values outside of allowed ranges
        if np.min(np.squeeze(freqs)) < self.valid_freq_range.min:
            raise ValueError(
                f"User specified frequency {np.min(np.squeeze(freqs))} is below the Medium's limit, {self.valid_freq_range.min}."
            )

        if np.max(np.squeeze(freqs)) > self.valid_freq_range.max:
            raise ValueError(
                f"User specified frequency {np.max(np.squeeze(freqs))} is above the Medium's limit, {self.valid_freq_range.max}."
            )

        # Initialize with instantaneous dielectric tensor
        epsmu = np.expand_dims(Matrix(diag=diag, offdiag=offdiag), axis=0)

        # Iterate through susceptibilities
        for i_sus in range(len(susceptibilities)):
            epsmu = epsmu + susceptibilities[i_sus].eval_susceptibility(freqs)

        # Account for conductivity term (only multiply if nonzero to avoid unnecessary complex numbers)
        conductivity = np.expand_dims(
            Matrix(diag=conductivity_diag, offdiag=conductivity_offdiag), axis=0
        )
        if np.count_nonzero(conductivity) > 0:
            epsmu = (1 + 1j / (2 * np.pi * freqs) * conductivity) * epsmu

        # Convert list matrix to 3D numpy array size [freqs,3,3]
        return np.squeeze(epsmu)


class MaterialGrid:
    """
    This class is used to specify materials on a rectilinear grid. A class object is passed
    as the `material` argument of a [`Block`](#block) geometric object or the `default_material`
    argument of the [`Simulation`](#Simulation) constructor (similar to a [material function](#medium)).
    """

    def check_weights(self, w):
        if np.amin(w) >= 0.0 and np.amax(w) <= 1.0:
            return w
        warnings.warn(
            "The weights parameter of MaterialGrid must be in the range [0,1]."
        )
        return np.clip(w, 0.0, 1.0)

    def __init__(
        self,
        grid_size: Union[Vector3, Tuple[float, ...]],
        medium1: Medium,
        medium2: Medium,
        weights: np.ndarray = None,
        grid_type: str = "U_DEFAULT",
        do_averaging: bool = True,
        beta: float = 0,
        eta: float = 0.5,
        damping: float = 0,
    ):
        """
        Creates a `MaterialGrid` object.

        The input are two materials `medium1` and `medium2` along with a weight function $u(x)$ which
        is defined on a rectilinear grid by the NumPy array `weights` of size `grid_size` (a 3-tuple or
        `Vector3` of integers $N_x$,$N_y$,$N_z$). The resolution of the grid may be nonuniform depending
        on the `size` property of the `Block` object as shown in the following example for a 2d `MaterialGrid`
        with $N_x=5$ and $N_y=4$. $N_z=0$ implies that the `MaterialGrid` is extruded in the $z$ direction.
        The grid points are defined at the corners of the voxels.

        ![](images/material_grid.png#center)

        Elements of the `weights` array must be in the range [0,1] where 0 is `medium1` and 1 is `medium2`.
        An array of boolean values `False` and `True` will be converted to 0 and 1, respectively.
        The `weights` array is used to define a linear interpolation from `medium1` to `medium2`.
        Two material types are supported: (1) frequency-independent isotropic $\\varepsilon$ (`epsilon_diag`
        and `epsilon_offdiag` are interpolated) and (2) `LorentzianSusceptibility` (`sigma` and `sigma_offdiag`
        are interpolated). `medium1` and `medium2` must both be the same type. The materials are
        [bilinearly interpolated](https://en.wikipedia.org/wiki/Bilinear_interpolation) from the rectilinear
        grid to Meep's [Yee grid](Yee_Lattice.md).

        For improving accuracy, [subpixel smoothing](Subpixel_Smoothing.md) can be enabled by specifying
        `do_averaging=True`. If you want to use a material grid to define a (nearly) discontinuous,
        piecewise-constant material that is *either* `medium1` or `medium2` almost everywhere, you can
        optionally enable a (smoothed) *projection* feature by setting the parameter `beta` to a
        positive value. The default is no projection (`beta=0`). When the projection feature is
        enabled, the weights $u(x)$ can be thought of as a
        [level-set function](https://en.wikipedia.org/wiki/Level-set_method) defining an interface at
        $u(x)=\\eta$ with a smoothing factor $\\beta$ where $\\beta=+\\infty$ gives an unsmoothed,
        discontinuous interface. The projection operator is $(\\tanh(\\beta\\times\\eta)
        +\\tanh(\\beta\\times(u-\\eta)))/(\\tanh(\\beta\\times\\eta)+\\tanh(\\beta\\times(1-\\eta)))$
        involving the parameters `beta` ($\\beta$: bias or "smoothness" of the turn on) and `eta`
        ($\\eta$: offset for erosion/dilation). The level set provides a general approach for defining
        a *discontinuous* function from otherwise continuously varying (via the bilinear interpolation)
        grid values. Subpixel smoothing is fast and accurate because it exploits an analytic formulation
        for level-set functions. Note that when subpixel smoothing is enabled via `do_averaging=True`,
        projecting the `weights` is done internally using the `beta` parameter. It is therefore not
        necessary to manually project the `weights` outside of `MaterialGrid`. However, visualizing
        the `weights` used to define the structure does require manually projecting the `weights` yourself.
        (Alternatively, you can output the actual structure using [`plot2D`](#data-visualization) or
        [`output_epsilon`](#output-functions_1).)

        A nonzero `damping` term creates an artificial conductivity $\\sigma = u(1-u)*$`damping`, which acts as
        dissipation loss that penalizes intermediate pixel values of non-binarized structures. The value of
        `damping` should be proportional to $2\\pi$ times the typical frequency of the problem.

        It is possible to overlap any number of different `MaterialGrid`s. This can be useful for defining
        grids which are symmetric (e.g., mirror, rotation). One way to set this up is by overlapping a
        given `MaterialGrid` object with a symmetrized copy of itself. In the case of spatially overlapping
        `MaterialGrid` objects (with no intervening objects), any overlapping points are computed using the
        method `grid_type` which is one of `"U_MIN"` (minimum of the overlapping grid values), `"U_PROD"`
        (product), `"U_MEAN"` (mean), `"U_DEFAULT"` (topmost material grid). In general, these `"U_*"` options
        allow you to combine any material grids that overlap in space with no intervening objects.
        """
        self.grid_size = mp.Vector3(*grid_size)
        self.medium1 = medium1
        self.medium2 = medium2

        def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
            return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

        if isclose(self.grid_size.x, 0):
            self.grid_size.x = 1
        if isclose(self.grid_size.y, 0):
            self.grid_size.y = 1
        if isclose(self.grid_size.z, 0):
            self.grid_size.z = 1
        self.num_params = int(self.grid_size.x * self.grid_size.y * self.grid_size.z)
        self.do_averaging = do_averaging
        self.beta = beta
        self.eta = eta
        self.damping = damping
        if weights is None:
            self.weights = np.zeros((self.num_params,))
        elif weights.size != self.num_params:
            raise ValueError(
                "weights of shape {} do not match user specified grid dimension: {}".format(
                    weights.size, self.grid_size
                )
            )
        else:
            self.weights = self.check_weights(weights.flatten().astype(np.float64))

        grid_type_dict = {"U_MIN": 0, "U_PROD": 1, "U_MEAN": 2, "U_DEFAULT": 3}
        if grid_type not in grid_type_dict:
            raise ValueError(
                "Invalid grid_type: {}. Must be either U_MIN, U_PROD, U_MEAN, or U_DEFAULT".format(
                    grid_type_dict
                )
            )
        self.grid_type = grid_type_dict[grid_type]

        self.swigobj = None

    def update_weights(self, x: np.ndarray):
        """
        Reset the `weights` to `x`.
        """
        if x.size != self.num_params:
            raise ValueError(
                f"weights of shape {self.weights.size} do not match user specified grid dimension: {self.grid_size}"
            )

        self.weights[:] = self.check_weights(x).flatten().astype(np.float64)


class Susceptibility:
    """
    Parent class for various dispersive susceptibility terms, parameterized by an
    anisotropic amplitude $\\sigma$. See [Material Dispersion](Materials.md#material-dispersion).
    """

    def __init__(self, sigma_diag=Vector3(), sigma_offdiag=Vector3(), sigma=None):
        """
        + **`sigma` [`number`]** — The scale factor $\\sigma$.

        You can also specify an anisotropic $\\sigma$ tensor by using the property `sigma_diag`
        which takes three numbers or a `Vector3` to give the $\\sigma_n$ tensor diagonal, and
        `sigma_offdiag` which specifies the offdiagonal elements (defaults to 0). That is,
        `sigma_diag=mp.Vector3(a, b, c)` and `sigma_offdiag=mp.Vector3(u, v, w)`
        corresponds to a $\\sigma$ tensor

        \\begin{pmatrix} a & u & v \\\\ u & b & w \\\\ v & w & c \\end{pmatrix}
        """
        self.sigma_diag = (
            Vector3(sigma, sigma, sigma) if sigma else Vector3(*sigma_diag)
        )
        self.sigma_offdiag = Vector3(*sigma_offdiag)

    def transform(self, m):
        sigma = Matrix(diag=self.sigma_diag, offdiag=self.sigma_offdiag)
        new_sigma = (m * sigma * m.transpose()) / abs(m.determinant())
        self.sigma_diag = mp.Vector3(new_sigma.c1.x, new_sigma.c2.y, new_sigma.c3.z)
        self.sigma_offdiag = mp.Vector3(new_sigma.c2.x, new_sigma.c3.x, new_sigma.c3.y)


class LorentzianSusceptibility(Susceptibility):
    """
    Specifies a single dispersive susceptibility of Lorentzian (damped harmonic
    oscillator) form. See [Material Dispersion](Materials.md#material-dispersion), with
    the parameters (in addition to $\\sigma$):
    """

    def __init__(self, frequency=0.0, gamma=0.0, **kwargs):
        """
        + **`frequency` [`number`]** — The resonance frequency $f_n = \\omega_n / 2\\pi$.

        + **`gamma` [`number`]** — The resonance loss rate $\\gamma_n / 2\\pi$.

        Note: multiple objects with identical values for the `frequency` and `gamma` but
        different `sigma` will appear as a *single* Lorentzian susceptibility term in the
        preliminary simulation info output.
        """
        super().__init__(**kwargs)
        self.frequency = frequency
        self.gamma = gamma

    def eval_susceptibility(self, freq):
        sigma = np.expand_dims(
            Matrix(diag=self.sigma_diag, offdiag=self.sigma_offdiag), axis=0
        )
        if self.gamma == 0:
            return (
                self.frequency
                * self.frequency
                / (self.frequency * self.frequency - freq * freq)
                * sigma
            )
        else:
            return (
                self.frequency
                * self.frequency
                / (
                    self.frequency * self.frequency
                    - freq * freq
                    - 1j * self.gamma * freq
                )
                * sigma
            )


class DrudeSusceptibility(Susceptibility):
    """
    Specifies a single dispersive susceptibility of Drude form. See [Material
    Dispersion](Materials.md#material-dispersion), with the parameters (in addition to $\\sigma$):
    """

    def __init__(self, frequency=0.0, gamma=0.0, **kwargs):
        """
        + **`frequency` [`number`]** — The frequency scale factor $f_n = \\omega_n / 2\\pi$
          which multiplies $\\sigma$ (not a resonance frequency).

        + **`gamma` [`number`]** — The loss rate $\\gamma_n / 2\\pi$.
        """
        super().__init__(**kwargs)
        self.frequency = frequency
        self.gamma = gamma

    def eval_susceptibility(self, freq):
        sigma = np.expand_dims(
            Matrix(diag=self.sigma_diag, offdiag=self.sigma_offdiag), axis=0
        )
        if self.gamma == 0:
            return -self.frequency * self.frequency / (freq * (freq)) * sigma
        else:
            return (
                -self.frequency
                * self.frequency
                / (freq * (freq + 1j * self.gamma))
                * sigma
            )


class NoisyLorentzianSusceptibility(LorentzianSusceptibility):
    """
    Specifies a single dispersive susceptibility of Lorentzian (damped harmonic
    oscillator) or Drude form. See [Material
    Dispersion](Materials.md#material-dispersion), with the same `sigma`, `frequency`, and
    `gamma` parameters, but with an additional Gaussian random noise term (uncorrelated in
    space and time, zero mean) added to the **P** damped-oscillator equation.
    """

    def __init__(self, noise_amp=0.0, **kwargs):
        """
        + **`noise_amp` [`number`]** — The noise has root-mean square amplitude σ $\\times$
          `noise_amp`.

        This is a somewhat unusual polarizable medium, a Lorentzian susceptibility with a
        random noise term added into the damped-oscillator equation at each point. This
        can be used to directly model thermal radiation in both the [far
        field](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.93.213905) and the
        [near field](http://math.mit.edu/~stevenj/papers/RodriguezIl11.pdf). Note, however
        that it is more efficient to [compute far-field thermal radiation using
        Kirchhoff's law](http://www.simpetus.com/projects.html#meep_thermal_radiation) of
        radiation, which states that emissivity equals absorptivity. Near-field thermal
        radiation can usually be computed more efficiently using frequency-domain methods,
        e.g. via [SCUFF-EM](https://github.com/HomerReid/scuff-em), as described e.g.
        [here](http://doi.org/10.1103/PhysRevB.92.134202) or
        [here](http://doi.org/10.1103/PhysRevB.88.054305).
        """
        super().__init__(**kwargs)
        self.noise_amp = noise_amp


class NoisyDrudeSusceptibility(DrudeSusceptibility):
    """
    Specifies a single dispersive susceptibility of Lorentzian (damped harmonic
    oscillator) or Drude form. See [Material
    Dispersion](Materials.md#material-dispersion), with the same `sigma`, `frequency`, and
    `gamma` parameters, but with an additional Gaussian random noise term (uncorrelated in
    space and time, zero mean) added to the **P** damped-oscillator equation.
    """

    def __init__(self, noise_amp=0.0, **kwargs):
        """
        + **`noise_amp` [`number`]** — The noise has root-mean square amplitude σ $\\times$
          `noise_amp`.

        This is a somewhat unusual polarizable medium, a Lorentzian susceptibility with a
        random noise term added into the damped-oscillator equation at each point. This
        can be used to directly model thermal radiation in both the [far
        field](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.93.213905) and the
        [near field](http://math.mit.edu/~stevenj/papers/RodriguezIl11.pdf). Note, however
        that it is more efficient to [compute far-field thermal radiation using
        Kirchhoff's law](http://www.simpetus.com/projects.html#meep_thermal_radiation) of
        radiation, which states that emissivity equals absorptivity. Near-field thermal
        radiation can usually be computed more efficiently using frequency-domain methods,
        e.g. via [SCUFF-EM](https://github.com/HomerReid/scuff-em), as described e.g.
        [here](http://doi.org/10.1103/PhysRevB.92.134202) or
        [here](http://doi.org/10.1103/PhysRevB.88.054305).
        """
        super().__init__(**kwargs)
        self.noise_amp = noise_amp


class GyrotropicLorentzianSusceptibility(LorentzianSusceptibility):
    """
    (**Experimental feature**) Specifies a single dispersive [gyrotropic
    susceptibility](Materials.md#gyrotropic-media) of [Lorentzian (damped harmonic
    oscillator) or Drude form](Materials.md#gyrotropic-drude-lorentz-model). Its
    parameters are `sigma`, `frequency`, and `gamma`, which have the [usual
    meanings](#susceptibility), and an additional 3-vector `bias`:
    """

    def __init__(self, bias=Vector3(), **kwargs):
        """
        + **`bias` [`Vector3`]** — The gyrotropy vector.  Its direction determines the
          orientation of the gyrotropic response, and the magnitude is the precession
          frequency $|\\mathbf{b}_n|/2\\pi$.
        """
        super().__init__(**kwargs)
        self.bias = bias


class GyrotropicDrudeSusceptibility(DrudeSusceptibility):
    """
    (**Experimental feature**) Specifies a single dispersive [gyrotropic
    susceptibility](Materials.md#gyrotropic-media) of [Lorentzian (damped harmonic
    oscillator) or Drude form](Materials.md#gyrotropic-drude-lorentz-model). Its
    parameters are `sigma`, `frequency`, and `gamma`, which have the [usual
    meanings](#susceptibility), and an additional 3-vector `bias`:
    """

    def __init__(self, bias=Vector3(), **kwargs):
        """
        + **`bias` [`Vector3`]** — The gyrotropy vector.  Its direction determines the
          orientation of the gyrotropic response, and the magnitude is the precession
          frequency $|\\mathbf{b}_n|/2\\pi$.
        """
        super().__init__(**kwargs)
        self.bias = bias


class GyrotropicSaturatedSusceptibility(Susceptibility):
    """
    (**Experimental feature**) Specifies a single dispersive [gyrotropic
    susceptibility](Materials.md#gyrotropic-media) governed by a [linearized
    Landau-Lifshitz-Gilbert
    equation](Materials.md#gyrotropic-saturated-dipole-linearized-landau-lifshitz-gilbert-model).
    This class takes parameters `sigma`, `frequency`, and `gamma`, whose meanings are
    different from the Lorentzian and Drude case. It also takes a 3-vector `bias`
    parameter and an `alpha` parameter:
    """

    def __init__(self, bias=Vector3(), frequency=0.0, gamma=0.0, alpha=0.0, **kwargs):
        """
        + **`sigma` [`number`]** — The coupling factor $\\sigma_n / 2\\pi$ between the
          polarization and the driving field. In [magnetic
          ferrites](https://en.wikipedia.org/wiki/Ferrite_(magnet)), this is the Larmor
          precession frequency at the saturation field.

        + **`frequency` [`number`]** — The [Larmor
          precession](https://en.wikipedia.org/wiki/Larmor_precession) frequency,
          $f_n = \\omega_n / 2\\pi$.

        + **`gamma` [`number`]** — The loss rate $\\gamma_n / 2\\pi$ in the off-diagonal
          response.

        + **`alpha` [`number`]** — The loss factor $\\alpha_n$ in the diagonal response.
          Note that this parameter is dimensionless and contains no 2π factor.

        + **`bias` [`Vector3`]** — Vector specifying the orientation of the gyrotropic
          response. Unlike the similarly-named `bias` parameter for the [gyrotropic
          Lorentzian/Drude
          susceptibilities](#gyrotropiclorentziansusceptibility-or-gyrotropicdrudesusceptibility),
          the magnitude is ignored; instead, the relevant precession frequencies are
          determined by the `sigma` and `frequency` parameters.
        """
        super().__init__(**kwargs)
        self.frequency = frequency
        self.gamma = gamma
        self.bias = bias
        self.alpha = alpha


class MultilevelAtom(Susceptibility):
    """
    Specifies a multievel atomic susceptibility for modeling saturable gain and
    absorption. This is a subclass of `E_susceptibilities` which contains two objects: (1)
    `transitions`: a list of atomic `Transition`s (defined below), and (2)
    `initial_populations`: a list of numbers defining the initial population of each
    atomic level. See [Materials/Saturable Gain and
    Absorption](Materials.md#saturable-gain-and-absorption).
    """

    def __init__(self, initial_populations=None, transitions=None, **kwargs):
        super().__init__(**kwargs)
        self.initial_populations = initial_populations or []
        self.transitions = transitions or []


class Transition:
    """ """

    def __init__(
        self,
        from_level,
        to_level,
        transition_rate=0,
        frequency=0,
        sigma_diag=Vector3(1, 1, 1),
        gamma=0,
        pumping_rate=0,
    ):
        """
        Construct a `Transition`.

        + **`frequency` [`number`]** — The radiative transition frequency $f = \\omega / 2\\pi$.

        + **`gamma` [`number`]** — The loss rate $\\gamma = \\gamma / 2\\pi$.

        + **`sigma_diag` [`Vector3`]** — The per-polarization coupling strength $\\sigma$.

        + **`from_level` [`number`]** — The atomic level from which the transition occurs.

        + **`to_level` [`number`]** — The atomic level to which the transition occurs.

        + **`transition_rate` [`number`]** — The non-radiative transition rate
          $f = \\omega / 2\\pi$. Default is 0.

        + **`pumping_rate` [`number`]** — The pumping rate $f = \\omega / 2\\pi$. Default is 0.
        """
        self.from_level = check_nonnegative("from_level", from_level)
        self.to_level = check_nonnegative("to_level", to_level)
        self.transition_rate = transition_rate
        self.frequency = frequency
        self.sigma_diag = sigma_diag
        self.gamma = gamma
        self.pumping_rate = pumping_rate


class GeometricObject:
    """
    This class, and its descendants, are used to specify the solid geometric objects that
    form the dielectric structure being simulated.

    In a 2d calculation, only the intersections of the objects with the $xy$ plane are
    considered.

    **Geometry Utilities**

    See the [MPB documentation](https://mpb.readthedocs.io/en/latest/Python_User_Interface/#geometry-utilities)
    for utility functions to help manipulate geometric objects.

    **Examples**

    These are some examples of geometric objects created using some `GeometricObject`
    subclasses:

    ```python
    # A cylinder of infinite radius and height 0.25 pointing along the x axis,
    # centered at the origin:
    cyl = mp.Cylinder(center=mp.Vector3(0,0,0), height=0.25, radius=mp.inf,
                    axis=mp.Vector3(1,0,0), material=mp.Medium(index=3.5))
    ```

    ```python
    # An ellipsoid with its long axis pointing along (1,1,1), centered on
    # the origin (the other two axes are orthogonal and have equal semi-axis lengths):
    ell = mp.Ellipsoid(center=mp.Vector3(0,0,0), size=mp.Vector3(0.8,0.2,0.2),
                    e1=Vector3(1,1,1), e2=Vector3(0,1,-1), e3=Vector3(-2,1,1),
                    material=mp.Medium(epsilon=13))
    ```

    ```python
    # A unit cube of material metal with a spherical air hole of radius 0.2 at
    # its center, the whole thing centered at (1,2,3):
    geometry=[mp.Block(center=Vector3(1,2,3), size=Vector3(1,1,1), material=mp.metal),
            mp.Sphere(center=Vector3(1,2,3), radius=0.2, material=mp.air)]
    ```

    ```python
    # A hexagonal prism defined by six vertices centered on the origin
    # of material crystalline silicon (from the materials library)
    vertices = [mp.Vector3(-1,0),
                mp.Vector3(-0.5,math.sqrt(3)/2),
                mp.Vector3(0.5,math.sqrt(3)/2),
                mp.Vector3(1,0),
                mp.Vector3(0.5,-math.sqrt(3)/2),
                mp.Vector3(-0.5,-math.sqrt(3)/2)]

    geometry = [mp.Prism(vertices, height=1.5, center=mp.Vector3(), material=cSi)]
    ```
    """

    def __init__(self, material=Medium(), center=Vector3(), epsilon_func=None):
        """
        Construct a `GeometricObject`.

        + **`material` [`Medium` class or function ]** — The material that the object is
          made of (usually some sort of dielectric). Uses default `Medium`. If a function
          is supplied, it must take one argument and return a Python `Medium`.

        + **`epsilon_func` [ function ]** — A function that takes one argument (a
          `Vector3`) and returns the dielectric constant at that point. Can be used
          instead of `material`. Default is `None`.

        + **`center` [`Vector3`]** — Center point of the object. Defaults to `(0,0,0)`.

        One normally does not create objects of type `GeometricObject` directly, however;
        instead, you use one of the following subclasses. Recall that subclasses inherit
        the properties of their superclass, so these subclasses automatically have the
        `material` and `center` properties and can be specified in a subclass's
        constructor via keyword arguments.
        """
        if type(material) is not Medium and callable(material):
            init_do_averaging(material)
            material.eps = False
        elif epsilon_func:
            init_do_averaging(epsilon_func)
            epsilon_func.eps = True
            material = epsilon_func

        self.material = material
        self.center = Vector3(*center)

    def __contains__(self, point):
        return mp.is_point_in_object(Vector3(*point), self)

    def __add__(self, vec):
        return self.shift(Vector3(*vec))

    def __radd__(self, vec):
        return self.shift(Vector3(*vec))

    def __iadd__(self, vec):
        self.center += Vector3(*vec)
        return self

    def shift(self, vec):
        """
        Shifts the object's `center` by `vec` (`Vector3`), returning a new object.
        This can also be accomplished via the `+` operator:

        ```python
        geometric_obj + Vector3(10,10,10)
        ```

        Using `+=` will shift the object in place.
        """
        c = deepcopy(self)
        c.center += Vector3(*vec)
        return c

    def info(self, indent_by=0):
        """
        Displays all properties and current values of a `GeometricObject`, indented by
        `indent_by` spaces (default is 0).
        """
        mp.display_geometric_object_info(indent_by, self)


class Sphere(GeometricObject):
    """
    Represents a sphere.

    **Properties:**

    + **`radius` [`number`]** — Radius of the sphere. No default value.
    """

    def __init__(self, radius, **kwargs):
        """Constructs a `Sphere`"""
        self.radius = float(radius)
        super().__init__(**kwargs)

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, val):
        self._radius = check_nonnegative("Sphere.radius", val)


class Cylinder(GeometricObject):
    """
    A cylinder, with circular cross-section and finite height.

    **Properties:**

    + **`radius` [`number`]** — Radius of the cylinder's cross-section. No default value.

    + **`height` [`number`]** — Length of the cylinder along its axis. No default value.

    + **`axis` [`Vector3`]** — Direction of the cylinder's axis; the length of this vector
      is ignored. Defaults to `Vector3(x=0, y=0, z=1)`.
    """

    def __init__(self, radius, axis=Vector3(0, 0, 1), height=1e20, **kwargs):
        """
        Constructs a `Cylinder`.
        """
        self.axis = Vector3(*axis)
        self.radius = float(radius)
        self.height = float(height)
        super().__init__(**kwargs)

    @property
    def radius(self):
        return self._radius

    @property
    def height(self):
        return self._height

    @radius.setter
    def radius(self, val):
        self._radius = check_nonnegative("Cylinder.radius", val)

    @height.setter
    def height(self, val):
        self._height = check_nonnegative("Cylinder.height", val)


class Wedge(Cylinder):
    """
    Represents a cylindrical wedge.
    """

    def __init__(
        self, radius, wedge_angle=2 * math.pi, wedge_start=Vector3(1, 0, 0), **kwargs
    ):
        """
        Constructs a `Wedge`.
        """
        self.wedge_angle = wedge_angle
        self.wedge_start = Vector3(*wedge_start)
        super().__init__(radius, **kwargs)


class Cone(Cylinder):
    """
    A cone, or possibly a truncated cone. This is actually a subclass of `Cylinder`, and
    inherits all of the same properties, with one additional property. The radius of the
    base of the cone is given by the `radius` property inherited from `Cylinder`, while
    the radius of the tip is given by the new property, `radius2`. The `center` of a cone
    is halfway between the two circular ends.
    """

    def __init__(self, radius, radius2=0, **kwargs):
        """
        Construct a `Cone`.

        **`radius2` [`number`]**
        —
        Radius of the tip of the cone (i.e. the end of the cone pointed to by the `axis` vector). Defaults to zero (a "sharp" cone).
        """
        self.radius2 = radius2
        super().__init__(radius, **kwargs)


class Block(GeometricObject):
    """
    A parallelepiped (i.e., a brick, possibly with non-orthogonal axes).
    """

    def __init__(
        self,
        size,
        e1=Vector3(1, 0, 0),
        e2=Vector3(0, 1, 0),
        e3=Vector3(0, 0, 1),
        **kwargs,
    ):
        """
        Construct a `Block`.

        + **`size` [`Vector3`]** — The lengths of the block edges along each of its three
          axes. Not really a 3-vector, but it has three components, each of which should
          be nonzero. No default value.

        + **`e1`, `e2`, `e3` [`Vector3`]** — The directions of the axes of the block; the
          lengths of these vectors are ignored. Must be linearly independent. They default
          to the three lattice directions.
        """

        self.size = Vector3(*size)
        self.e1 = Vector3(*e1)
        self.e2 = Vector3(*e2)
        self.e3 = Vector3(*e3)
        super().__init__(**kwargs)


class Ellipsoid(Block):
    """
    An ellipsoid. This is actually a subclass of `Block`, and inherits all the same
    properties, but defines an ellipsoid inscribed inside the block.
    """

    def __init__(self, **kwargs):
        """
        Construct an `Ellipsoid`.
        """
        super().__init__(**kwargs)


class Prism(GeometricObject):
    """
    Polygonal prism type.
    """

    def __init__(
        self,
        vertices,
        height,
        axis=Vector3(z=1),
        center=None,
        sidewall_angle=0,
        **kwargs,
    ):
        """
        Construct a `Prism`.

        + **`vertices` [list of `Vector3`]** — The vertices that make up the prism. They
          must lie in a plane that's perpendicular to the `axis`. Note that infinite
          lengths are not supported. To simulate infinite geometry, just extend the edge
          of the prism beyond the cell.

        + **`height` [`number`]** — The prism thickness, extruded in the direction of
          `axis`. `mp.inf` can be used for infinite height. No default value.

        + **`axis` [`Vector3`]** — The axis perpendicular to the prism. Defaults to
          `Vector3(0,0,1)`.

        + **`center` [`Vector3`]** — If `center` is not specified, then the coordinates of
          the `vertices` define the *bottom* of the prism with the top of the prism being
          at the same coordinates shifted by `height*axis`. If `center` is specified, then
          `center` is the coordinates of the
          [centroid](https://en.wikipedia.org/wiki/Centroid) of all the vertices (top and
          bottom) of the resulting 3d prism so that the coordinates of the `vertices` are
          shifted accordingly.

        + **`sidewall_angle` [`number`]** — The sidewall angle of the prism in units of
          radians. Default is 0.
        """

        centroid = sum(vertices, Vector3(0)) * (
            1.0 / len(vertices)
        )  # centroid of floor polygon
        original_center = (
            centroid + (0.5 * height) * axis
        )  # center as computed from vertices, height, axis
        if center is not None and len(vertices):
            center = Vector3(*center)
            # translate vertices to center prism at requested center
            shift = center - original_center
            vertices = list(map(lambda v: v + shift, vertices))
        else:
            center = original_center
        self.vertices = vertices
        self.height = height
        self.axis = axis
        self.sidewall_angle = sidewall_angle

        super().__init__(center=center, **kwargs)


class Matrix:
    """
    The `Matrix` class represents a 3x3 matrix with c1, c2, and c3 as its columns.

    ```python
    m.transpose()
    m.getH() or m.H
    m.determinant()
    m.inverse()
    ```

    Return the transpose, adjoint (conjugate transpose), determinant, or inverse of the
    given matrix.

    ```python
    m1 + m2
    m1 - m2
    m1 * m2
    ```

    Return the sum, difference, or product of the given matrices.

    ```python
    v * m
    m * v
    ```

    Returns the `Vector3` product of the matrix `m` by the vector `v`, with the vector
    multiplied on the left or the right respectively.

    ```python
    s * m
    m * s
    ```

    Scales the matrix `m` by the number `s`.
    """

    def __init__(
        self,
        c1=Vector3(),
        c2=Vector3(),
        c3=Vector3(),
        diag=Vector3(),
        offdiag=Vector3(),
    ):
        """
        Constructs a `Matrix`.
        """
        self.c1 = Vector3(*c1)
        self.c2 = Vector3(*c2)
        self.c3 = Vector3(*c3)
        if np.all(c1 == c2) and np.all(c2 == c3) and np.all(c3 == Vector3()):
            self.c1 = Vector3(diag[0], offdiag[0], offdiag[1])
            self.c2 = Vector3(np.conj(offdiag[0]), diag[1], offdiag[2])
            self.c3 = Vector3(np.conj(offdiag[1]), np.conj(offdiag[2]), diag[2])

    def __getitem__(self, i):
        return self.row(i)

    def __mul__(self, m):
        if type(m) is Matrix:
            return self.mm_mult(m)
        elif type(m) is Vector3:
            return self.mv_mult(m)
        elif isinstance(m, Number):
            return self.scale(m)
        else:
            raise TypeError(f"No operation known for 'Matrix * {type(m)}'")

    def __rmul__(self, left_arg):
        if isinstance(left_arg, Number):
            return self.scale(left_arg)
        else:
            raise TypeError(f"No operation known for 'Matrix * {type(left_arg)}'")

    def __truediv__(self, scalar):
        return Matrix(self.c1 / scalar, self.c2 / scalar, self.c3 / scalar)

    def __add__(self, m):
        return Matrix(self.c1 + m.c1, self.c2 + m.c2, self.c3 + m.c3)

    def __sub__(self, m):
        return Matrix(self.c1 - m.c1, self.c2 - m.c2, self.c3 - m.c3)

    def __repr__(self):
        r0 = self.row(0)
        r1 = self.row(1)
        r2 = self.row(2)
        return f"<<{r0[0]} {r0[1]} {r0[2]}>\n <{r1[0]} {r1[1]} {r1[2]}>\n <{r2[0]} {r2[1]} {r2[2]}>>"

    def __array__(self):
        return np.array(
            [self.row(0).__array__(), self.row(1).__array__(), self.row(2).__array__()]
        )

    def row(self, i):
        return Vector3(self.c1[i], self.c2[i], self.c3[i])

    def mm_mult(self, m):
        c1 = Vector3(
            self.row(0).dot(m.c1), self.row(1).dot(m.c1), self.row(2).dot(m.c1)
        )
        c2 = Vector3(
            self.row(0).dot(m.c2), self.row(1).dot(m.c2), self.row(2).dot(m.c2)
        )
        c3 = Vector3(
            self.row(0).dot(m.c3), self.row(1).dot(m.c3), self.row(2).dot(m.c3)
        )

        return Matrix(c1, c2, c3)

    def mv_mult(self, v):
        return Vector3(*[self.row(i).dot(Vector3(*v)) for i in range(3)])

    def scale(self, s):
        return Matrix(self.c1.scale(s), self.c2.scale(s), self.c3.scale(s))

    def determinant(self):
        sum1 = sum(
            [
                functools.reduce(operator.mul, [self[x][x] for x in range(3)]),
                functools.reduce(operator.mul, [self[0][1], self[1][2], self[2][0]]),
                functools.reduce(operator.mul, [self[1][0], self[2][1], self[0][2]]),
            ]
        )
        sum2 = sum(
            [
                functools.reduce(operator.mul, [self[0][2], self[1][1], self[2][0]]),
                functools.reduce(operator.mul, [self[0][1], self[1][0], self[2][2]]),
                functools.reduce(operator.mul, [self[1][2], self[2][1], self[0][0]]),
            ]
        )
        return sum1 - sum2

    def conj(self):
        return Matrix(self.c1.conj(), self.c2.conj(), self.c3.conj())

    def transpose(self):
        return Matrix(self.row(0), self.row(1), self.row(2))

    def getH(self):
        return self.transpose().conj()

    def inverse(self):
        v1x = self[1][1] * self[2][2] - self[1][2] * self[2][1]
        v1y = self[1][2] * self[2][0] - self[1][0] * self[2][2]
        v1z = self[1][0] * self[2][1] - self[1][1] * self[2][0]
        v1 = mp.Vector3(v1x, v1y, v1z)

        v2x = self[2][1] * self[0][2] - self[0][1] * self[2][2]
        v2y = self[0][0] * self[2][2] - self[0][2] * self[2][0]
        v2z = self[0][1] * self[2][0] - self[0][0] * self[2][1]
        v2 = mp.Vector3(v2x, v2y, v2z)

        v3x = self[0][1] * self[1][2] - self[1][1] * self[0][2]
        v3y = self[1][0] * self[0][2] - self[0][0] * self[1][2]
        v3z = self[1][1] * self[0][0] - self[1][0] * self[0][1]
        v3 = mp.Vector3(v3x, v3y, v3z)

        m = Matrix(v1, v2, v3)

        return m.scale(1 / self.determinant())

    H = property(getH, None)


class Lattice:
    def __init__(
        self,
        size=Vector3(1, 1, 1),
        basis_size=Vector3(1, 1, 1),
        basis1=Vector3(1, 0, 0),
        basis2=Vector3(0, 1, 0),
        basis3=Vector3(0, 0, 1),
    ):

        self.size = Vector3(*size)
        self.basis_size = Vector3(*basis_size)
        self.basis1 = Vector3(*basis1)
        self.basis2 = Vector3(*basis2)
        self.basis3 = Vector3(*basis3)

    @property
    def basis1(self):
        return self._basis1

    @basis1.setter
    def basis1(self, val):
        self._basis1 = val.unit()

    @property
    def basis2(self):
        return self._basis2

    @basis2.setter
    def basis2(self, val):
        self._basis2 = val.unit()

    @property
    def basis3(self):
        return self._basis3

    @basis3.setter
    def basis3(self, val):
        self._basis3 = val.unit()

    @property
    def b1(self):
        return self.basis1.scale(self.basis_size.x)

    @property
    def b2(self):
        return self.basis2.scale(self.basis_size.y)

    @property
    def b3(self):
        return self.basis3.scale(self.basis_size.z)

    @property
    def basis(self):
        B = Matrix(self.b1, self.b2, self.b3)

        if B.determinant() == 0:
            raise ValueError("Lattice basis vectors must be linearly independent.")

        return B

    @property
    def metric(self):
        B = self.basis
        return B.transpose() * B


def lattice_to_cartesian(x, lat):
    if isinstance(x, Vector3):
        return lat.basis * x

    return (lat.basis * x) * lat.basis.inverse()


def cartesian_to_lattice(x, lat):
    if isinstance(x, Vector3):
        return lat.basis.inverse() * x

    return (lat.basis.inverse() * x) * lat.basis


def reciprocal_to_cartesian(x, lat):
    s = Vector3(*[1 if v == 0 else v for v in lat.size])

    m = Matrix(Vector3(s.x), Vector3(y=s.y), Vector3(z=s.z))
    Rst = (lat.basis * m).transpose()

    return Rst.inverse() * x if isinstance(x, Vector3) else (Rst.inverse() * x) * Rst


def cartesian_to_reciprocal(x, lat):
    s = Vector3(*[1 if v == 0 else v for v in lat.size])

    m = Matrix(Vector3(s.x), Vector3(y=s.y), Vector3(z=s.z))
    Rst = (lat.basis * m).transpose()

    return Rst * x if isinstance(x, Vector3) else (Rst * x) * Rst.inverse()


def lattice_to_reciprocal(x, lat):
    return cartesian_to_reciprocal(lattice_to_cartesian(x, lat), lat)


def reciprocal_to_lattice(x, lat):
    return cartesian_to_lattice(reciprocal_to_cartesian(x, lat), lat)


def geometric_object_duplicates(shift_vector, min_multiple, max_multiple, go):

    shift_vector = Vector3(*shift_vector)

    def _dup(min_multiple, lst):
        if min_multiple > max_multiple:
            return lst
        shifted = go.shift(shift_vector.scale(min_multiple))
        return _dup(min_multiple + 1, [shifted] + lst)

    return _dup(min_multiple, [])


def geometric_objects_duplicates(shift_vector, min_multiple, max_multiple, go_list):
    dups = []
    shift_vector = Vector3(*shift_vector)
    for go in go_list:
        dups += geometric_object_duplicates(
            shift_vector, min_multiple, max_multiple, go
        )
    return dups


def geometric_objects_lattice_duplicates(lat, go_list, *usize):
    def lat_to_lattice(v):
        return cartesian_to_lattice(lat.basis * v, lat)

    u1 = usize[0] if usize else 1
    u2 = usize[1] if len(usize) >= 2 else 1
    u3 = usize[2] if len(usize) >= 3 else 1
    s = lat.size

    b1 = lat_to_lattice(mp.Vector3(u1))
    b2 = lat_to_lattice(mp.Vector3(0, u2, 0))
    b3 = lat_to_lattice(mp.Vector3(0, 0, u3))

    n1 = math.ceil((s.x if s.x else 1e-20) / u1)
    n2 = math.ceil((s.y if s.y else 1e-20) / u2)
    n3 = math.ceil((s.z if s.z else 1e-20) / u3)

    min3 = -math.floor((n3 - 1) / 2)
    max3 = math.ceil((n3 - 1) / 2)
    d3 = geometric_objects_duplicates(b3, int(min3), int(max3), go_list)

    min2 = -math.floor((n2 - 1) / 2)
    max2 = math.ceil((n2 - 1) / 2)
    d2 = geometric_objects_duplicates(b2, int(min2), int(max2), d3)

    min1 = -math.floor((n1 - 1) / 2)
    max1 = math.ceil((n1 - 1) / 2)

    return geometric_objects_duplicates(b1, int(min1), int(max1), d2)


# Return a 'memoized' version of the function f, which caches its
# arguments and return values so as never to compute the same thing twice.
def memoize(f):
    f_memo_tab = {}

    def _mem(y=None):
        tab_val = f_memo_tab.get(y, None)
        if tab_val:
            return tab_val

        fy = f(y)
        f_memo_tab[y] = fy
        return fy

    return _mem


# Find a root by Newton's method with bounds and bisection,
# given a function f that returns a pair of (value . derivative)
def find_root_deriv(f, tol, x_min, x_max, x_guess=None):
    # Some trickiness: we only need to evaluate the function at x_min and
    # x_max if a Newton step fails, and even then only if we haven't already
    # bracketed the root, so do this via lazy evaluation.
    f_memo = memoize(f)

    def lazy(x):
        return x if isinstance(x, numbers.Number) else x()

    def pick_bound(which):
        def _pb():
            fmin_tup = f_memo(x_min)
            fmax_tup = f_memo(x_max)
            fmin = fmin_tup[0]
            fmax = fmax_tup[0]

            if which(fmin):
                return x_min
            elif which(fmax):
                return x_max
            else:
                raise ValueError("failed to bracket the root in find_root_deriv")

        return _pb

    def in_bounds(x, f, df, a, b):
        return (f - (df * (x - a))) * (f - (df * (x - b))) < 0

    def newton(x, a, b, dx):
        if abs(dx) < abs(tol * x):
            return x

        fx_tup = f_memo(x)
        f = fx_tup[0]
        df = fx_tup[1]

        if f == 0:
            return x

        a_prime = x if f < 0 else a
        b_prime = x if f > 0 else b

        if (
            dx != x_max - x_min
            and dx * (f / df) < 0
            and f_memo(lazy(a_prime))[0] * f_memo(lazy(b_prime))[0] > 0
        ):
            raise ValueError("failed to bracket the root in find_root_deriv")

        if isinstance(a, numbers.Number) and isinstance(b, numbers.Number):
            is_in_bounds = in_bounds(x, f, df, a, b)
        else:
            is_in_bounds = in_bounds(x, f, df, x_min, x_max)

        if is_in_bounds:
            return newton(x - (f / df), a_prime, b_prime, f / df)

        av = lazy(a)
        bv = lazy(b)
        dx_prime = 0.5 * (bv - av)
        a_pp = av if a == a_prime else a_prime
        b_pp = bv if b == b_prime else b_prime

        return newton((av + bv) * 0.5, a_pp, b_pp, dx_prime)

    if x_guess is None:
        x_guess = (x_min + x_max) * 0.5

    return newton(
        x_guess,
        pick_bound(lambda aa: aa < 0),
        pick_bound(lambda aa: aa > 0),
        x_max - x_min,
    )


def get_rotation_matrix(axis, theta):
    """
    Like `Vector3.rotate`, except returns the (unitary) rotation matrix that performs the
    given rotation. i.e., `get_rotation_matrix(axis, theta) * v` produces the same result
    as `v.rotate(axis, theta)`.

    + `axis` [`Vector3`] — The vector around which the rotation is applied in the right-hand direction.

    + `theta` [`number`] — The rotation angle (in radians).
    """
    return Matrix(
        Vector3(x=1).rotate(axis, theta),
        Vector3(y=1).rotate(axis, theta),
        Vector3(z=1).rotate(axis, theta),
    )
