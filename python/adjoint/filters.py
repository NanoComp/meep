"""
A collection of routines for use in topology optimization comprising
convolution filters, projection operators, and morphological transforms.
"""
from typing import List, Tuple, Union

from autograd import numpy as npa
import meep as mp
import numpy as np
from scipy import signal, special

ArrayLikeType = Union[List, Tuple, np.ndarray]


def _centered(arr: np.ndarray, newshape: ArrayLikeType) -> np.ndarray:
    """Reformats the zero-padded array of the FFT filter.

    A helper function borrowed from SciPy:
    https://github.com/scipy/scipy/blob/v1.4.1/scipy/signal/signaltools.py#L263-L270

    Args:
        arr: the zero-padded array to be used as the FFT filter.
        newshape: a 1d array containing the dimensions of the new array.

    Returns:
        The reformatted array with center based on the new array dimensions.
    """
    newshape = np.asarray(newshape)
    currshape = np.array(arr.shape)
    startind = (currshape - newshape) // 2
    endind = startind + newshape
    myslice = [slice(startind[k], endind[k]) for k in range(len(endind))]

    return arr[tuple(myslice)]


def _proper_pad(arr: np.ndarray, pad_to: np.ndarray) -> np.ndarray:
    """Zero-pads the array as preprocessing for convolution filter.

    Args:
        arr: 2d input array representing the nonnegative coordinates of a
            filter kernel with C4v symmetry.
        pad_to: 1d array of two elements (integers) indicating the size
            of the zero-padded array.

    Returns:
        A 2d array with zero padding.
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


def _edge_pad(arr: np.ndarray, pad: np.ndarray) -> np.ndarray:
    """Zero-pads the edges of the array as preprocessing for convolution filter.

    Args:
        arr: 2d array representing the nonnegative coordinates of a
            filter kernel with C4v symmetry.
        pad_to: 2x2 array of integers indicating the size
            of the zero-padded array.

    Returns:
        A 2d array with zero padding.
    """
    # fill sides
    left = npa.tile(arr[0, :], (pad[0][0], 1))
    right = npa.tile(arr[-1, :], (pad[0][1], 1))
    top = npa.tile(arr[:, 0], (pad[1][0], 1)).transpose()
    bottom = npa.tile(arr[:, -1], (pad[1][1], 1)).transpose()

    # fill corners
    top_left = npa.tile(arr[0, 0], (pad[0][0], pad[1][0]))
    top_right = npa.tile(arr[-1, 0], (pad[0][1], pad[1][0]))
    bottom_left = npa.tile(arr[0, -1], (pad[0][0], pad[1][1]))
    bottom_right = npa.tile(arr[-1, -1], (pad[0][1], pad[1][1]))

    if pad[0][0] > 0 and pad[0][1] > 0 and pad[1][0] > 0 and pad[1][1] > 0:
        return npa.concatenate(
            (
                npa.concatenate((top_left, top, top_right)),
                npa.concatenate((left, arr, right)),
                npa.concatenate((bottom_left, bottom, bottom_right)),
            ),
            axis=1,
        )
    elif pad[0][0] == 0 and pad[0][1] == 0 and pad[1][0] > 0 and pad[1][1] > 0:
        return npa.concatenate((top, arr, bottom), axis=1)
    elif pad[0][0] > 0 and pad[0][1] > 0 and pad[1][0] == 0 and pad[1][1] == 0:
        return npa.concatenate((left, arr, right), axis=0)
    elif pad[0][0] == 0 and pad[0][1] == 0 and pad[1][0] == 0 and pad[1][1] == 0:
        return arr
    else:
        raise ValueError("At least one of the padding numbers is invalid.")


def simple_2d_filter(
    x: np.ndarray, h: np.ndarray, periodic_axes: ArrayLikeType = None
) -> np.ndarray:
    """A simple 2d filter function.

    Uses a 2d FFT which is typically faster and preserves the shape
    of the input and output arrays. The FFTs pad the operation to prevent
    unwanted effects from the edges of the array.

    Args:
        x: 2d design weights (unfiltered).
        h: filter kernel prior to applying the FFT. Must be same size as `x`
        periodic_axes: list of axes (x, y = 0, 1) that are to be treated as
            periodic. Default is None (all axes are non-periodic).

    Returns:
        The output of the 2d convolution.
    """
    (sx, sy) = x.shape

    if periodic_axes is None:
        h = _proper_pad(h, 3 * np.array([sx, sy]))
        x = _edge_pad(x, ((sx, sx), (sy, sy)))
    else:
        (kx, ky) = h.shape

        npx = int(
            np.ceil((2 * kx - 1) / sx)
        )  # 2*kx-1 is the size of a complete kernel in the x direction
        npy = int(
            np.ceil((2 * ky - 1) / sy)
        )  # 2*ky-1 is the size of a complete kernel in the y direction
        if npx % 2 == 0:
            npx += 1  # Ensure npx is an odd number
        if npy % 2 == 0:
            npy += 1  # Ensure npy is an odd number

        # Repeat the design pattern in periodic directions according to
        # the kernel size
        x = npa.tile(
            x, (npx if 0 in periodic_axes else 1, npy if 1 in periodic_axes else 1)
        )

        npadx = 0 if 0 in periodic_axes else sx
        npady = 0 if 1 in periodic_axes else sy
        x = _edge_pad(
            x, ((npadx, npadx), (npady, npady))
        )  # pad only in nonperiodic directions
        h = _proper_pad(
            h,
            np.array(
                [
                    npx * sx if 0 in periodic_axes else 3 * sx,
                    npy * sy if 1 in periodic_axes else 3 * sy,
                ]
            ),
        )

    h = h / npa.sum(h)  # Normalize the kernel

    return _centered(
        npa.real(npa.fft.ifft2(npa.fft.fft2(x) * npa.fft.fft2(h))), (sx, sy)
    )


def cylindrical_filter(
    x: np.ndarray,
    radius: float,
    Lx: float,
    Ly: float,
    resolution: int,
    periodic_axes: ArrayLikeType = None,
) -> np.ndarray:
    """A uniform cylindrical filter.

    Typically allows for sharper features.

    Ref: B.S. Lazarov, F. Wang, & O. Sigmund, Length scale and
    manufacturability in density-based topology optimization,
    Archive of Applied Mechanics, 86(1-2), pp. 189-218 (2016).

    Args:
        x: 2d design weights to be filtered.
        radius: filter radius (in Meep units).
        Lx: length of design region in X direction (in Meep units).
        Ly: length of design region in Y direction (in Meep units).
        resolution: resolution of the design grid (not the Meep grid
            resolution)
        periodic_axes: list of axes (x, y = 0, 1) that are to be treated as
            periodic. Default is None (all axes are non-periodic).

    Returns:
        The filtered design weights.
    """
    Nx = int(round(Lx * resolution)) + 1
    Ny = int(round(Ly * resolution)) + 1
    x = x.reshape(Nx, Ny)  # Ensure the input is 2D

    xv = np.arange(0, Lx / 2, 1 / resolution)
    yv = np.arange(0, Ly / 2, 1 / resolution)

    # If the design weights are periodic in a direction,
    # the size of the kernel in that direction needs to be adjusted
    # according to the filter radius.
    if periodic_axes is not None:
        periodic_axes = np.array(periodic_axes)
        if 0 in periodic_axes:
            xv = np.arange(0, np.ceil(2 * radius / Lx) * Lx / 2, 1 / resolution)
        if 1 in periodic_axes:
            yv = np.arange(0, np.ceil(2 * radius / Ly) * Ly / 2, 1 / resolution)

    X, Y = np.meshgrid(xv, yv, sparse=True, indexing="ij")
    h = np.where(X**2 + Y**2 < radius**2, 1, 0)

    return simple_2d_filter(x, h, periodic_axes)


def conic_filter(
    x: np.ndarray,
    radius: float,
    Lx: float,
    Ly: float,
    resolution: int,
    periodic_axes: ArrayLikeType = None,
) -> np.ndarray:
    """A linear conic (or "hat") filter.

    Ref: B.S. Lazarov, F. Wang, & O. Sigmund, Length scale and
    manufacturability in density-based topology optimization.
    Archive of Applied Mechanics, 86(1-2), pp. 189-218 (2016).

    Args:
        x: 2d design weights (unfiltered).
        radius: filter radius (in Meep units).
        Lx: length of design region in X direction (in Meep units).
        Ly: length of design region in Y direction (in Meep units).
        resolution: resolution of the design grid (not the Meep grid
            resolution).
        periodic_axes: list of axes (x, y = 0, 1) that are to be treated as
            periodic. Default is None (all axes are non-periodic).

    Returns:
        The filtered design weights.
    """
    Nx = int(round(Lx * resolution)) + 1
    Ny = int(round(Ly * resolution)) + 1
    x = x.reshape(Nx, Ny)  # Ensure the input is 2D

    xv = np.arange(0, Lx / 2, 1 / resolution)
    yv = np.arange(0, Ly / 2, 1 / resolution)

    # If the design pattern is periodic in a direction,
    # the size of the kernel in that direction needs to be adjusted
    # according to the filter radius.
    if periodic_axes is not None:
        periodic_axes = np.array(periodic_axes)
        if 0 in periodic_axes:
            xv = np.arange(0, np.ceil(2 * radius / Lx) * Lx / 2, 1 / resolution)
        if 1 in periodic_axes:
            yv = np.arange(0, np.ceil(2 * radius / Ly) * Ly / 2, 1 / resolution)

    X, Y = np.meshgrid(xv, yv, sparse=True, indexing="ij")
    h = np.where(
        X**2 + Y**2 < radius**2, (1 - np.sqrt(abs(X**2 + Y**2)) / radius), 0
    )

    return simple_2d_filter(x, h, periodic_axes)


def gaussian_filter(
    x: np.ndarray,
    sigma: float,
    Lx: float,
    Ly: float,
    resolution: int,
    periodic_axes: ArrayLikeType = None,
):
    """A simple Gaussian filter of the form exp(-x**2 / sigma**2).

    Ref: E. W. Wang, D. Sell, T. Phan, & J. A. Fan, Robust design of
    topology-optimized metasurfaces, Optical Materials Express, 9(2),
    pp. 469-482 (2019).

    Args:
        x: 2d design weights (unfiltered).
        sigma: filter radius (in Meep units).
        Lx: length of design region in X direction (in Meep units).
        Ly: length of design region in Y direction (in Meep units).
        resolution: resolution of the design grid (not the Meep grid
            resolution).
        periodic_axes: list of axes (x, y = 0, 1) that are to be treated as
            periodic. Default is None (all axes are non-periodic).

    Returns:
        The filtered design weights.
    """
    Nx = int(round(Lx * resolution)) + 1
    Ny = int(round(Ly * resolution)) + 1
    x = x.reshape(Nx, Ny)  # Ensure the input is 2D

    xv = np.arange(0, Lx / 2, 1 / resolution)
    yv = np.arange(0, Ly / 2, 1 / resolution)

    # If the design pattern is periodic in a direction,
    # the size of the kernel in that direction needs to be
    # adjusted according to 3σ.
    if periodic_axes is not None:
        periodic_axes = np.array(periodic_axes)
        if 0 in periodic_axes:
            xv = np.arange(0, np.ceil(6 * sigma / Lx) * Lx / 2, 1 / resolution)
        if 1 in periodic_axes:
            yv = np.arange(0, np.ceil(6 * sigma / Ly) * Ly / 2, 1 / resolution)

    X, Y = np.meshgrid(xv, yv, sparse=True, indexing="ij")
    h = np.exp(-(X**2 + Y**2) / sigma**2)

    return simple_2d_filter(x, h, periodic_axes)


def exponential_erosion(
    x: np.ndarray,
    radius: float,
    beta: float,
    Lx: float,
    Ly: float,
    resolution: int,
    periodic_axes: ArrayLikeType = None,
):
    """Morphological erosion using exponential projection.

    Refs:
    O. Sigmund, Morphology-based black and white filters for topology
    optimization. Structural and Multidisciplinary Optimization,
    33(4-5), pp. 401-424 (2007).
    M. Schevenels & O. Sigmund, On the implementation and effectiveness of
    morphological close-open and open-close filters for topology optimization.
    Structural and Multidisciplinary Optimization, 54(1), pp. 15-21 (2016).

    Args:
        x: 2d design weights (unfiltered).
        radius: filter radius (in Meep units).
        beta: threshold value for projection. Range of [0, inf].
        Lx: length of design region in X direction (in Meep units).
        Ly: length of design region in Y direction (in Meep units).
        resolution: resolution of the design grid (not the Meep grid
            resolution).
        periodic_axes: list of axes (x, y = 0, 1) that are to be treated as
            periodic. Default is None (all axes are non-periodic).

    Returns:
        The eroded design weights.
    """
    x_hat = npa.exp(beta * (1 - x))
    return (
        1
        - npa.log(
            cylindrical_filter(
                x_hat, radius, Lx, Ly, resolution, periodic_axes
            ).flatten()
        )
        / beta
    )


def exponential_dilation(x, radius, beta, Lx, Ly, resolution, periodic_axes=None):
    """Morphological dilation using exponential projection.

    Refs:
    O. Sigmund, Morphology-based black and white filters for topology
    optimization. Structural and Multidisciplinary Optimization,
    33(4-5), pp. 401-424 (2007).
    M. Schevenels & O. Sigmund, On the implementation and effectiveness of
    morphological close-open and open-close filters for topology optimization.
    Structural and Multidisciplinary Optimization, 54(1), pp. 15-21 (2016).

    Args:
        x: 2d design weights (unfiltered).
        radius: filter radius (in Meep units).
        beta: threshold value for projection. Range of [0, inf].
        Lx: length of design region in X direction (in Meep units).
        Ly: length of design region in Y direction (in Meep units).
        resolution: resolution of the design grid (not the Meep grid
            resolution).
        periodic_axes: list of axes (x, y = 0, 1) that are to be treated as
            periodic. Default is None (all axes are non-periodic).

    Returns:
        The dilated design weights.
    """
    x_hat = npa.exp(beta * x)
    return (
        npa.log(
            cylindrical_filter(
                x_hat, radius, Lx, Ly, resolution, periodic_axes
            ).flatten()
        )
        / beta
    )


def heaviside_erosion(x, radius, beta, Lx, Ly, resolution, periodic_axes=None):
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
    periodic_axes: array_like (1D)
        List of axes (x, y = 0, 1) that are to be treated as periodic (default is none: all axes are non-periodic)

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

    x_hat = cylindrical_filter(x, radius, Lx, Ly, resolution, periodic_axes).flatten()
    return npa.exp(-beta * (1 - x_hat)) + npa.exp(-beta) * (1 - x_hat)


def heaviside_dilation(x, radius, beta, Lx, Ly, resolution, periodic_axes=None):
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
    periodic_axes: array_like (1D)
        List of axes (x, y = 0, 1) that are to be treated as periodic (default is none: all axes are non-periodic)

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

    x_hat = cylindrical_filter(x, radius, Lx, Ly, resolution, periodic_axes).flatten()
    return 1 - npa.exp(-beta * x_hat) + npa.exp(-beta) * x_hat


def geometric_erosion(x, radius, alpha, Lx, Ly, resolution, periodic_axes=None):
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
    periodic_axes: array_like (1D)
        List of axes (x, y = 0, 1) that are to be treated as periodic (default is none: all axes are non-periodic)

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
        npa.exp(
            cylindrical_filter(x_hat, radius, Lx, Ly, resolution, periodic_axes)
        ).flatten()
        - alpha
    )


def geometric_dilation(x, radius, alpha, Lx, Ly, resolution, periodic_axes=None):
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
    periodic_axes: array_like (1D)
        List of axes (x, y = 0, 1) that are to be treated as periodic (default is none: all axes are non-periodic)

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
        -npa.exp(
            cylindrical_filter(x_hat, radius, Lx, Ly, resolution, periodic_axes)
        ).flatten()
        + alpha
        + 1
    )


def harmonic_erosion(x, radius, alpha, Lx, Ly, resolution, periodic_axes=None):
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
    periodic_axes: array_like (1D)
        List of axes (x, y = 0, 1) that are to be treated as periodic (default is none: all axes are non-periodic)

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
    return (
        1
        / cylindrical_filter(x_hat, radius, Lx, Ly, resolution, periodic_axes).flatten()
        - alpha
    )


def harmonic_dilation(x, radius, alpha, Lx, Ly, resolution, periodic_axes=None):
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
    periodic_axes: array_like (1D)
        List of axes (x, y = 0, 1) that are to be treated as periodic (default is none: all axes are non-periodic)

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
        1
        - 1
        / cylindrical_filter(x_hat, radius, Lx, Ly, resolution, periodic_axes).flatten()
        + alpha
    )


def tanh_projection(x: np.ndarray, beta: float, eta: float) -> np.ndarray:
    """Sigmoid projection filter.

    Ref: F. Wang, B. S. Lazarov, & O. Sigmund, On projection methods,
    convergence and robust formulations in topology optimization.
    Structural and Multidisciplinary Optimization, 43(6), pp. 767-784 (2011).

    Args:
        x: 2d design weights to be filtered.
        beta: thresholding parameter in the range [0, inf]. Determines the
            degree of binarization of the output.
        eta: threshold point in the range [0, 1].

    Returns:
        The filtered design weights.
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


def indicator_solid(x, c, filter_f, threshold_f, resolution, periodic_axes=None):
    """Calculates the indicator function for the void phase needed for minimum length optimization [1].

    Parameters
    ----------
    x : array_like
        Design parameters
    c : float
        Decay rate parameter (1e0 - 1e8)
    filter_f : function_handle
        Filter function. Must be differntiable by autograd.
    threshold_f : function_handle
        Threshold function. Must be differntiable by autograd.
    periodic_axes: array_like (1D)
        List of axes (x, y = 0, 1) that are to be treated as periodic (default is none: all axes are non-periodic)

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

    if periodic_axes is None:
        gradient_filtered_field = npa.gradient(filtered_field)
    else:
        periodic_axes = np.array(periodic_axes)
        if 0 in periodic_axes:
            filtered_field = npa.tile(filtered_field, (3, 1))
        if 1 in periodic_axes:
            filtered_field = npa.tile(filtered_field, (1, 3))
        gradient_filtered_field = _centered(
            npa.array(npa.gradient(filtered_field)), (2,) + x.shape
        )

    grad_mag = (gradient_filtered_field[0] * resolution) ** 2 + (
        gradient_filtered_field[1] * resolution
    ) ** 2
    if grad_mag.ndim != 2:
        raise ValueError(
            "The gradient fields must be 2 dimensional. Check input array and filter functions."
        )
    return design_field * npa.exp(-c * grad_mag)


def constraint_solid(
    x, c, eta_e, filter_f, threshold_f, resolution, periodic_axes=None
):
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
    periodic_axes: array_like (1D)
        List of axes (x, y = 0, 1) that are to be treated as periodic (default is none: all axes are non-periodic)

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
        x.reshape(filtered_field.shape),
        c,
        filter_f,
        threshold_f,
        resolution,
        periodic_axes,
    ).flatten()
    return npa.mean(I_s * npa.minimum(filtered_field.flatten() - eta_e, 0) ** 2)


def indicator_void(x, c, filter_f, threshold_f, resolution, periodic_axes=None):
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
    periodic_axes: array_like (1D)
        List of axes (x, y = 0, 1) that are to be treated as periodic (default is none: all axes are non-periodic)

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

    if periodic_axes is None:
        gradient_filtered_field = npa.gradient(filtered_field)
    else:
        periodic_axes = np.array(periodic_axes)
        if 0 in periodic_axes:
            filtered_field = npa.tile(filtered_field, (3, 1))
        if 1 in periodic_axes:
            filtered_field = npa.tile(filtered_field, (1, 3))
        gradient_filtered_field = _centered(
            npa.array(npa.gradient(filtered_field)), (2,) + x.shape
        )

    grad_mag = (gradient_filtered_field[0] * resolution) ** 2 + (
        gradient_filtered_field[1] * resolution
    ) ** 2
    if grad_mag.ndim != 2:
        raise ValueError(
            "The gradient fields must be 2 dimensional. Check input array and filter functions."
        )
    return (1 - design_field) * npa.exp(-c * grad_mag)


def constraint_void(x, c, eta_d, filter_f, threshold_f, resolution, periodic_axes=None):
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
    periodic_axes: array_like (1D)
        List of axes (x, y = 0, 1) that are to be treated as periodic (default is none: all axes are non-periodic).

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
        x.reshape(filtered_field.shape),
        c,
        filter_f,
        threshold_f,
        resolution,
        periodic_axes,
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
