"""
General filter functions to be used in other projection and morphological transform routines.
"""
import meep as mp
import numpy as np
from autograd import numpy as npa
from scipy import signal, special


def _centered(arr, newshape):
    """Helper function that reformats the padded array of the fft filter operation.
    Borrowed from scipy:
    https://github.com/scipy/scipy/blob/v1.4.1/scipy/signal/signaltools.py#L263-L270
    """
    # Return the center newshape portion of the array.
    newshape = np.asarray(newshape)
    currshape = np.array(arr.shape)
    startind = (currshape - newshape) // 2
    endind = startind + newshape
    myslice = [slice(startind[k], endind[k]) for k in range(len(endind))]
    return arr[tuple(myslice)]


def _proper_pad(arr, pad_to):
    """Return a zero-padded/expanded version of filter arr (≈1/4 of kernel) to size pad_to, for convolution.

    Parameters
    ----------
    arr : 2d input array representing the nonnegative-coordinate ≈1/4 of a filter kernel.
    pad_to : 1d array composed of two integers indicating the total size to be padded to.
    """

    pad_size = pad_to - 2 * np.array(arr.shape) + 1

    top = np.zeros((pad_size[0], arr.shape[1]))
    bottom = np.zeros((pad_size[0], arr.shape[1] - 1))
    middle = np.zeros((pad_to[0], pad_size[1]))

    top_left = arr[:, :]
    top_right = npa.flipud(arr[1:, :])
    bottom_left = npa.fliplr(arr[:, 1:])
    bottom_right = npa.flipud(
        npa.fliplr(arr[1:, 1:])
    )  # equivalent to flip, but flip is incompatible with autograd

    return npa.concatenate(
        (
            npa.concatenate((top_left, top, top_right)),
            middle,
            npa.concatenate((bottom_left, bottom, bottom_right)),
        ),
        axis=1,
    )


def _edge_pad(arr, pad):

    # fill sides
    left = npa.tile(arr[0, :], (pad[0][0], 1))  # left side
    right = npa.tile(arr[-1, :], (pad[0][1], 1))  # right side
    top = npa.tile(arr[:, 0], (pad[1][0], 1)).transpose()  # top side
    bottom = npa.tile(arr[:, -1], (pad[1][1], 1)).transpose()  # bottom side)

    # fill corners
    top_left = npa.tile(arr[0, 0], (pad[0][0], pad[1][0]))  # top left
    top_right = npa.tile(arr[-1, 0], (pad[0][1], pad[1][0]))  # top right
    bottom_left = npa.tile(arr[0, -1], (pad[0][0], pad[1][1]))  # bottom left
    bottom_right = npa.tile(arr[-1, -1], (pad[0][1], pad[1][1]))  # bottom right

    return npa.concatenate(
        (
            npa.concatenate((top_left, top, top_right)),
            npa.concatenate((left, arr, right)),
            npa.concatenate((bottom_left, bottom, bottom_right)),
        ),
        axis=1,
    )


def simple_2d_filter(x, h):
    """A simple 2d filter algorithm that is differentiable with autograd.
    Uses a 2D fft approach since it is typically faster and preserves the shape
    of the input and output arrays.

    The ffts pad the operation to prevent any circular convolution garbage.

    Parameters
    ----------
    x : array_like (2D)
        Input array to be filtered. Must be 2D.
    h : array_like (2D)
        Filter kernel (before the DFT). Must be same size as `x`

    Returns
    -------
    array_like (2D)
        The output of the 2d convolution.
    """
    (kx, ky) = x.shape
    h = _proper_pad(h, 3 * np.array([kx, ky]))
    h = h / npa.sum(h)  # Normalize the kernel
    x = _edge_pad(x, ((kx, kx), (ky, ky)))

    return _centered(
        npa.real(npa.fft.ifft2(npa.fft.fft2(x) * npa.fft.fft2(h))), (kx, ky)
    )


def cylindrical_filter(x, radius, Lx, Ly, resolution):
    """A uniform cylindrical filter [1]. Typically allows for sharper transitions.

    Parameters
    ----------
    x : array_like (2D)
        Design parameters
    radius : float
        Filter radius (in "meep units")
    Lx : float
        Length of design region in X direction (in "meep units")
    Ly : float
        Length of design region in Y direction (in "meep units")
    resolution : int
        Resolution of the design grid (not the meep simulation resolution)

    Returns
    -------
    array_like (2D)
        Filtered design parameters.

    References
    ----------
    [1] Lazarov, B. S., Wang, F., & Sigmund, O. (2016). Length scale and manufacturability in
    density-based topology optimization. Archive of Applied Mechanics, 86(1-2), 189-218.
    """
    Nx = int(round(Lx * resolution))
    Ny = int(round(Ly * resolution))
    x = x.reshape(Nx, Ny)  # Ensure the input is 2D

    xv = np.arange(0, Lx / 2, 1 / resolution)
    yv = np.arange(0, Ly / 2, 1 / resolution)

    X, Y = np.meshgrid(xv, yv, sparse=True, indexing="ij")
    h = np.where(X**2 + Y**2 < radius**2, 1, 0)

    # Filter the response
    return simple_2d_filter(x, h)


def conic_filter(x, radius, Lx, Ly, resolution):
    """A linear conic filter, also known as a "Hat" filter in the literature [1].

    Parameters
    ----------
    x : array_like (2D)
        Design parameters
    radius : float
        Filter radius (in "meep units")
    Lx : float
        Length of design region in X direction (in "meep units")
    Ly : float
        Length of design region in Y direction (in "meep units")
    resolution : int
        Resolution of the design grid (not the meep simulation resolution)

    Returns
    -------
    array_like (2D)
        Filtered design parameters.

    References
    ----------
    [1] Lazarov, B. S., Wang, F., & Sigmund, O. (2016). Length scale and manufacturability in
    density-based topology optimization. Archive of Applied Mechanics, 86(1-2), 189-218.
    """
    Nx = int(round(Lx * resolution))
    Ny = int(round(Ly * resolution))
    x = x.reshape(Nx, Ny)  # Ensure the input is 2D

    xv = np.arange(0, Lx / 2, 1 / resolution)
    yv = np.arange(0, Ly / 2, 1 / resolution)

    X, Y = np.meshgrid(xv, yv, sparse=True, indexing="ij")
    h = np.where(
        X**2 + Y**2 < radius**2, (1 - np.sqrt(abs(X**2 + Y**2)) / radius), 0
    )

    # Filter the response
    return simple_2d_filter(x, h)


def gaussian_filter(x, sigma, Lx, Ly, resolution):
    """A simple gaussian filter of the form exp(-x **2 / sigma ** 2) [1].

    Parameters
    ----------
    x : array_like (2D)
        Design parameters
    sigma : float
        Filter radius (in "meep units")
    Lx : float
        Length of design region in X direction (in "meep units")
    Ly : float
        Length of design region in Y direction (in "meep units")
    resolution : int
        Resolution of the design grid (not the meep simulation resolution)

    Returns
    -------
    array_like (2D)
        Filtered design parameters.

    References
    ----------
    [1] Wang, E. W., Sell, D., Phan, T., & Fan, J. A. (2019). Robust design of
    topology-optimized metasurfaces. Optical Materials Express, 9(2), 469-482.
    """
    Nx = int(round(Lx * resolution))
    Ny = int(round(Ly * resolution))
    x = x.reshape(Nx, Ny)  # Ensure the input is 2D

    xv = np.arange(0, Lx / 2, 1 / resolution)
    yv = np.arange(0, Ly / 2, 1 / resolution)

    X, Y = np.meshgrid(xv, yv, sparse=True, indexing="ij")
    h = np.exp(-(X**2 + Y**2) / sigma**2)

    # Filter the response
    return simple_2d_filter(x, h)


"""
# ------------------------------------------------------------------------------------ #
Erosion and dilation operators
"""


def exponential_erosion(x, radius, beta, Lx, Ly, resolution):
    """Performs and exponential erosion operation.

    Parameters
    ----------
    x : array_like
        Design parameters
    radius : float
        Filter radius (in "meep units")
    beta : float
        Thresholding parameter
    Lx : float
        Length of design region in X direction (in "meep units")
    Ly : float
        Length of design region in Y direction (in "meep units")
    resolution : int
        Resolution of the design grid (not the meep simulation resolution)

    Returns
    -------
    array_like
        Eroded design parameters.

    References
    ----------
    [1] Sigmund, O. (2007). Morphology-based black and white filters for topology optimization.
    Structural and Multidisciplinary Optimization, 33(4-5), 401-424.
    [2] Schevenels, M., & Sigmund, O. (2016). On the implementation and effectiveness of
    morphological close-open and open-close filters for topology optimization. Structural
    and Multidisciplinary Optimization, 54(1), 15-21.
    """

    x_hat = npa.exp(beta * (1 - x))
    return (
        1
        - npa.log(cylindrical_filter(x_hat, radius, Lx, Ly, resolution).flatten())
        / beta
    )


def exponential_dilation(x, radius, beta, Lx, Ly, resolution):
    """Performs a exponential dilation operation.

    Parameters
    ----------
    x : array_like
        Design parameters
    radius : float
        Filter radius (in "meep units")
    beta : float
        Thresholding parameter
    Lx : float
        Length of design region in X direction (in "meep units")
    Ly : float
        Length of design region in Y direction (in "meep units")
    resolution : int
        Resolution of the design grid (not the meep simulation resolution)

    Returns
    -------
    array_like
        Dilated design parameters.

    References
    ----------
    [1] Sigmund, O. (2007). Morphology-based black and white filters for topology optimization.
    Structural and Multidisciplinary Optimization, 33(4-5), 401-424.
    [2] Schevenels, M., & Sigmund, O. (2016). On the implementation and effectiveness of
    morphological close-open and open-close filters for topology optimization. Structural
    and Multidisciplinary Optimization, 54(1), 15-21.
    """

    x_hat = npa.exp(beta * x)
    return (
        npa.log(cylindrical_filter(x_hat, radius, Lx, Ly, resolution).flatten()) / beta
    )


def heaviside_erosion(x, radius, beta, Lx, Ly, resolution):
    """Performs a heaviside erosion operation.

    Parameters
    ----------
    x : array_like
        Design parameters
    radius : float
        Filter radius (in "meep units")
    beta : float
        Thresholding parameter
    Lx : float
        Length of design region in X direction (in "meep units")
    Ly : float
        Length of design region in Y direction (in "meep units")
    resolution : int
        Resolution of the design grid (not the meep simulation resolution)

    Returns
    -------
    array_like
        Eroded design parameters.

    References
    ----------
    [1] Guest, J. K., Prévost, J. H., & Belytschko, T. (2004). Achieving minimum length scale in topology
    optimization using nodal design variables and projection functions. International journal for
    numerical methods in engineering, 61(2), 238-254.
    """

    x_hat = cylindrical_filter(x, radius, Lx, Ly, resolution).flatten()
    return npa.exp(-beta * (1 - x_hat)) + npa.exp(-beta) * (1 - x_hat)


def heaviside_dilation(x, radius, beta, Lx, Ly, resolution):
    """Performs a heaviside dilation operation.

    Parameters
    ----------
    x : array_like
        Design parameters
    radius : float
        Filter radius (in "meep units")
    beta : float
        Thresholding parameter
    Lx : float
        Length of design region in X direction (in "meep units")
    Ly : float
        Length of design region in Y direction (in "meep units")
    resolution : int
        Resolution of the design grid (not the meep simulation resolution)

    Returns
    -------
    array_like
        Dilated design parameters.

    References
    ----------
    [1] Guest, J. K., Prévost, J. H., & Belytschko, T. (2004). Achieving minimum length scale in topology
    optimization using nodal design variables and projection functions. International journal for
    numerical methods in engineering, 61(2), 238-254.
    """

    x_hat = cylindrical_filter(x, radius, Lx, Ly, resolution).flatten()
    return 1 - npa.exp(-beta * x_hat) + npa.exp(-beta) * x_hat


def geometric_erosion(x, radius, alpha, Lx, Ly, resolution):
    """Performs a geometric erosion operation.

    Parameters
    ----------
    x : array_like
        Design parameters
    radius : float
        Filter radius (in "meep units")
    beta : float
        Thresholding parameter
    Lx : float
        Length of design region in X direction (in "meep units")
    Ly : float
        Length of design region in Y direction (in "meep units")
    resolution : int
        Resolution of the design grid (not the meep simulation resolution)

    Returns
    -------
    array_like
        Eroded design parameters.

    References
    ----------
    [1] Svanberg, K., & Svärd, H. (2013). Density filters for topology optimization based on the
    Pythagorean means. Structural and Multidisciplinary Optimization, 48(5), 859-875.
    """
    x_hat = npa.log(x + alpha)
    return (
        npa.exp(cylindrical_filter(x_hat, radius, Lx, Ly, resolution)).flatten() - alpha
    )


def geometric_dilation(x, radius, alpha, Lx, Ly, resolution):
    """Performs a geometric dilation operation.

    Parameters
    ----------
    x : array_like
        Design parameters
    radius : float
        Filter radius (in "meep units")
    beta : float
        Thresholding parameter
    Lx : float
        Length of design region in X direction (in "meep units")
    Ly : float
        Length of design region in Y direction (in "meep units")
    resolution : int
        Resolution of the design grid (not the meep simulation resolution)

    Returns
    -------
    array_like
        Dilated design parameters.

    References
    ----------
    [1] Svanberg, K., & Svärd, H. (2013). Density filters for topology optimization based on the
    Pythagorean means. Structural and Multidisciplinary Optimization, 48(5), 859-875.
    """

    x_hat = npa.log(1 - x + alpha)
    return (
        -npa.exp(cylindrical_filter(x_hat, radius, Lx, Ly, resolution)).flatten()
        + alpha
        + 1
    )


def harmonic_erosion(x, radius, alpha, Lx, Ly, resolution):
    """Performs a harmonic erosion operation.

    Parameters
    ----------
    x : array_like
        Design parameters
    radius : float
        Filter radius (in "meep units")
    beta : float
        Thresholding parameter
    Lx : float
        Length of design region in X direction (in "meep units")
    Ly : float
        Length of design region in Y direction (in "meep units")
    resolution : int
        Resolution of the design grid (not the meep simulation resolution)

    Returns
    -------
    array_like
        Eroded design parameters.

    References
    ----------
    [1] Svanberg, K., & Svärd, H. (2013). Density filters for topology optimization based on the
    Pythagorean means. Structural and Multidisciplinary Optimization, 48(5), 859-875.
    """

    x_hat = 1 / (x + alpha)
    return 1 / cylindrical_filter(x_hat, radius, Lx, Ly, resolution).flatten() - alpha


def harmonic_dilation(x, radius, alpha, Lx, Ly, resolution):
    """Performs a harmonic dilation operation.

    Parameters
    ----------
    x : array_like
        Design parameters
    radius : float
        Filter radius (in "meep units")
    beta : float
        Thresholding parameter
    Lx : float
        Length of design region in X direction (in "meep units")
    Ly : float
        Length of design region in Y direction (in "meep units")
    resolution : int
        Resolution of the design grid (not the meep simulation resolution)

    Returns
    -------
    array_like
        Dilated design parameters.

    References
    ----------
    [1] Svanberg, K., & Svärd, H. (2013). Density filters for topology optimization based on the
    Pythagorean means. Structural and Multidisciplinary Optimization, 48(5), 859-875.
    """

    x_hat = 1 / (1 - x + alpha)
    return (
        1 - 1 / cylindrical_filter(x_hat, radius, Lx, Ly, resolution).flatten() + alpha
    )


"""
# ------------------------------------------------------------------------------------ #
Projection filters
"""


def tanh_projection(x, beta, eta):
    """Projection filter that thresholds the input parameters between 0 and 1. Typically
    the "strongest" projection.

    Parameters
    ----------
    x : array_like
        Design parameters
    beta : float
        Thresholding parameter (0 to infinity). Dictates how "binary" the output will be.
    eta: float
        Threshold point (0 to 1)

    Returns
    -------
    array_like
        Projected and flattened design parameters.
    References
    ----------
    [1] Wang, F., Lazarov, B. S., & Sigmund, O. (2011). On projection methods, convergence and robust
    formulations in topology optimization. Structural and Multidisciplinary Optimization, 43(6), 767-784.
    """

    return (npa.tanh(beta * eta) + npa.tanh(beta * (x - eta))) / (
        npa.tanh(beta * eta) + npa.tanh(beta * (1 - eta))
    )


def heaviside_projection(x, beta, eta):
    """Projection filter that thresholds the input parameters between 0 and 1.

    Parameters
    ----------
    x : array_like
        Design parameters
    beta : float
        Thresholding parameter (0 to infinity). Dictates how "binary" the output will be.
    eta: float
        Threshold point (0 to 1)

    Returns
    -------
    array_like
        Projected and flattened design parameters.

    References
    ----------
    [1] Lazarov, B. S., Wang, F., & Sigmund, O. (2016). Length scale and manufacturability in
    density-based topology optimization. Archive of Applied Mechanics, 86(1-2), 189-218.
    """

    case1 = eta * npa.exp(-beta * (eta - x) / eta) - (eta - x) * npa.exp(-beta)
    case2 = (
        1
        - (1 - eta) * npa.exp(-beta * (x - eta) / (1 - eta))
        - (eta - x) * npa.exp(-beta)
    )
    return npa.where(x < eta, case1, case2)


"""
# ------------------------------------------------------------------------------------ #
Length scale operations
"""


def get_threshold_wang(delta, sigma):
    """Calculates the threshold point according to the gaussian filter radius (`sigma`) and
    the perturbation parameter (`sigma`) needed to ensure the proper length
    scale and morphological transformation according to Wang et. al. [2].

    Parameters
    ----------
    sigma : float
        Smoothing radius (in meep units)
    delta : float
        Perturbation parameter (in meep units)

    Returns
    -------
    float
        Threshold point (`eta`)

    References
    ----------
    [1] Wang, F., Jensen, J. S., & Sigmund, O. (2011). Robust topology optimization of
    photonic crystal waveguides with tailored dispersion properties. JOSA B, 28(3), 387-397.
    [2] Wang, E. W., Sell, D., Phan, T., & Fan, J. A. (2019). Robust design of
    topology-optimized metasurfaces. Optical Materials Express, 9(2), 469-482.
    """

    return 0.5 - special.erf(delta / sigma)


def get_eta_from_conic(b, R):
    """Extracts the eroded threshold point (`eta_e`) for a conic filter given the desired
    minimum length (`b`) and the filter radius (`R`). This only works for conic filters.

    Note that the units for `b` and `R` can be arbitrary so long as they are consistent.

    Results in paper were thresholded using a "tanh" Heaviside projection.

    Parameters
    ----------
    b : float
        Desired minimum length scale.
    R : float
        Conic filter radius

    Returns
    -------
    float
        The eroded threshold point (1-eta)

    References
    ----------
    [1] Qian, X., & Sigmund, O. (2013). Topological design of electromechanical actuators with
    robustness toward over-and under-etching. Computer Methods in Applied
    Mechanics and Engineering, 253, 237-251.
    [2] Wang, F., Lazarov, B. S., & Sigmund, O. (2011). On projection methods, convergence and
    robust formulations in topology optimization. Structural and Multidisciplinary
    Optimization, 43(6), 767-784.
    [3] Lazarov, B. S., Wang, F., & Sigmund, O. (2016). Length scale and manufacturability in
    density-based topology optimization. Archive of Applied Mechanics, 86(1-2), 189-218.
    """

    norm_length = b / R
    if norm_length < 0:
        return 0
    elif norm_length < 1:
        return 0.25 * norm_length**2 + 0.5
    elif norm_length < 2:
        return -0.25 * norm_length**2 + norm_length
    else:
        return 1


def get_conic_radius_from_eta_e(b, eta_e):
    """Calculates the corresponding filter radius given the minimum length scale (b)
    and the desired eroded threshold point (eta_e).

    Parameters
    ----------
    b : float
        Desired minimum length scale.
    eta_e : float
        Eroded threshold point (1-eta)

    Returns
    -------
    float
        Conic filter radius.

    References
    ----------
    [1] Qian, X., & Sigmund, O. (2013). Topological design of electromechanical actuators with
    robustness toward over-and under-etching. Computer Methods in Applied
    Mechanics and Engineering, 253, 237-251.
    [2] Wang, F., Lazarov, B. S., & Sigmund, O. (2011). On projection methods, convergence and
    robust formulations in topology optimization. Structural and Multidisciplinary
    Optimization, 43(6), 767-784.
    [3] Lazarov, B. S., Wang, F., & Sigmund, O. (2016). Length scale and manufacturability in
    density-based topology optimization. Archive of Applied Mechanics, 86(1-2), 189-218.
    """
    if (eta_e >= 0.5) and (eta_e < 0.75):
        return b / (2 * np.sqrt(eta_e - 0.5))
    elif (eta_e >= 0.75) and (eta_e <= 1):
        return b / (2 - 2 * np.sqrt(1 - eta_e))
    else:
        raise ValueError(
            "The erosion threshold point (eta_e) must be between 0.5 and 1."
        )


def indicator_solid(x, c, filter_f, threshold_f, resolution):
    """Calculates the indicator function for the void phase needed for minimum length optimization [1].

    Parameters
    ----------
    x : array_like
        Design parameters
    c : float
        Decay rate parameter (1e0 - 1e8)
    eta_e : float
        Erosion threshold limit (0-1)
    filter_f : function_handle
        Filter function. Must be differntiable by autograd.
    threshold_f : function_handle
        Threshold function. Must be differntiable by autograd.

    Returns
    -------
    array_like
        Indicator value

    References
    ----------
    [1] Zhou, M., Lazarov, B. S., Wang, F., & Sigmund, O. (2015). Minimum length scale in topology optimization by
    geometric constraints. Computer Methods in Applied Mechanics and Engineering, 293, 266-282.
    """

    filtered_field = filter_f(x)
    design_field = threshold_f(filtered_field)
    gradient_filtered_field = npa.gradient(filtered_field)
    grad_mag = (gradient_filtered_field[0] * resolution) ** 2 + (
        gradient_filtered_field[1] * resolution
    ) ** 2
    if grad_mag.ndim != 2:
        raise ValueError(
            "The gradient fields must be 2 dimensional. Check input array and filter functions."
        )
    return design_field * npa.exp(-c * grad_mag)


def constraint_solid(x, c, eta_e, filter_f, threshold_f, resolution):
    """Calculates the constraint function of the solid phase needed for minimum length optimization [1].

    Parameters
    ----------
    x : array_like
        Design parameters
    c : float
        Decay rate parameter (1e0 - 1e8)
    eta_e : float
        Erosion threshold limit (0-1)
    filter_f : function_handle
        Filter function. Must be differntiable by autograd.
    threshold_f : function_handle
        Threshold function. Must be differntiable by autograd.

    Returns
    -------
    float
        Constraint value

    Example
    -------
    >> g_s = constraint_solid(x,c,eta_e,filter_f,threshold_f) # constraint
    >> g_s_grad = grad(constraint_solid,0)(x,c,eta_e,filter_f,threshold_f) # gradient

    References
    ----------
    [1] Zhou, M., Lazarov, B. S., Wang, F., & Sigmund, O. (2015). Minimum length scale in topology optimization by
    geometric constraints. Computer Methods in Applied Mechanics and Engineering, 293, 266-282.
    """

    filtered_field = filter_f(x)
    I_s = indicator_solid(
        x.reshape(filtered_field.shape), c, filter_f, threshold_f, resolution
    ).flatten()
    return npa.mean(I_s * npa.minimum(filtered_field.flatten() - eta_e, 0) ** 2)


def indicator_void(x, c, filter_f, threshold_f, resolution):
    """Calculates the indicator function for the void phase needed for minimum length optimization [1].

    Parameters
    ----------
    x : array_like
        Design parameters
    c : float
        Decay rate parameter (1e0 - 1e8)
    eta_d : float
        Dilation threshold limit (0-1)
    filter_f : function_handle
        Filter function. Must be differntiable by autograd.
    threshold_f : function_handle
        Threshold function. Must be differntiable by autograd.

    Returns
    -------
    array_like
        Indicator value

    References
    ----------
    [1] Zhou, M., Lazarov, B. S., Wang, F., & Sigmund, O. (2015). Minimum length scale in topology optimization by
    geometric constraints. Computer Methods in Applied Mechanics and Engineering, 293, 266-282.
    """

    filtered_field = filter_f(x).reshape(x.shape)
    design_field = threshold_f(filtered_field)
    gradient_filtered_field = npa.gradient(filtered_field)
    grad_mag = (gradient_filtered_field[0] * resolution) ** 2 + (
        gradient_filtered_field[1] * resolution
    ) ** 2
    if grad_mag.ndim != 2:
        raise ValueError(
            "The gradient fields must be 2 dimensional. Check input array and filter functions."
        )
    return (1 - design_field) * npa.exp(-c * grad_mag)


def constraint_void(x, c, eta_d, filter_f, threshold_f, resolution):
    """Calculates the constraint function of the void phase needed for minimum length optimization [1].

    Parameters
    ----------
    x : array_like
        Design parameters
    c : float
        Decay rate parameter (1e0 - 1e8)
    eta_d : float
        Dilation threshold limit (0-1)
    filter_f : function_handle
        Filter function. Must be differntiable by autograd.
    threshold_f : function_handle
        Threshold function. Must be differntiable by autograd.

    Returns
    -------
    float
        Constraint value

    Example
    -------
    >> g_v = constraint_void(p,c,eta_d,filter_f,threshold_f) # constraint
    >> g_v_grad = tensor_jacobian_product(constraint_void,0)(p,c,eta_d,filter_f,threshold_f,g_s) # gradient

    References
    ----------
    [1] Zhou, M., Lazarov, B. S., Wang, F., & Sigmund, O. (2015). Minimum length scale in topology optimization by
    geometric constraints. Computer Methods in Applied Mechanics and Engineering, 293, 266-282.
    """

    filtered_field = filter_f(x)
    I_v = indicator_void(
        x.reshape(filtered_field.shape), c, filter_f, threshold_f, resolution
    ).flatten()
    return npa.mean(I_v * npa.minimum(eta_d - filtered_field.flatten(), 0) ** 2)


def gray_indicator(x):
    """Calculates a measure of "grayness" according to [1].

    Lower numbers ( < 2%) indicate a good amount of binarization [1].

    Parameters
    ----------
    x : array_like
        Filtered and thresholded design parameters (between 0 and 1)

    Returns
    -------
    float
        Measure of "grayness" (in percent)

    References
    ----------
    [1] Lazarov, B. S., Wang, F., & Sigmund, O. (2016). Length scale and manufacturability in
    density-based topology optimization. Archive of Applied Mechanics, 86(1-2), 189-218.
    """
    return npa.mean(4 * x.flatten() * (1 - x.flatten())) * 100
