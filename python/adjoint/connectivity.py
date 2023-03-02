import numpy as np
from scipy.sparse.linalg import cg, spsolve
from scipy.sparse import kron, diags, csr_matrix, eye, csc_matrix
from autograd import numpy as npa
from autograd import grad
from typing import List

solvers = [spsolve, cg]


def constraint_connectivity(
    rho: List[float] = None,
    nx: int = None,
    ny: int = None,
    nz: int = None,
    cond_v: float = 1.0,
    cond_s: float = 1e4,
    src_v: float = 0.0,
    src_s: float = 1.0,
    solver_option: int = 0,
    thresh: float = 50.0,
    p: float = 3.0,
    need_grad: bool = True,
):
    """Computes its connectivity constraint value and the gradients with
    respect to each pixel value.

    Given a structure, assembles the finite difference matrices of the heat equation
    and solves for an auxiliary temperature field. Computes the p-norm of the field
    and compares it with the threshold to determine whether the structure is connected
    or not. Reference: "Li, Q., Chen, W., Liu, S. et al. Structural topology optimization considering
    connectivity constraint. Struct Multidisc Optim 54, 971–984 (2016).
    https://doi.org/10.1007/s00158-016-1459-5"

    Args:
        rho: the structure as (filtered and projected) design variables
        nx: size of the design variable in the x direction
        ny: size of the design variable in the y direction; should set to 1 for 2D structures
        nz: size of the design variable in the z direction; the supporting layer is outside the
          last slice of pixels in the z direction
        cond_v: heat conductivity for the void pixels; should NOT be changed
        cond_s: heat conductivity for solid pixels; IMPORTANT hyperparameter: changes the overall magnitude of
          heat constraint; good values generally between 1e3 to 1e5
        src_v: heat source value at void pixels; should NOT be changed
        src_s: heat source value at solid pixels; should NOT be changed
        solver_option: sparse solver option for solving the linear system. 0 for spsolve and 1 for cg
        thresh: threshold value against which the p-norm of the temperature field is compared. VERY IMPORTANT
          hyperparameter. Good values depend on cond_s and p values. Should be tuned given the problem setup and
          after choosing cond_s and p values.
        p: which p-norm of the temperature field to compute; good values generally between 3 to 5
        need_grad: True if gradients are needed; False if only want the forward constraint value

    Returns:
        If need_grad = True, returns T, heat, grad; if need_grad = False, only returns heat
        T is the computed auxiliary temperature field, heat is the constraint value (negative
        for connected and positive for disconnected structures), grad is the gradient of the
        constraint with respect to the design variables, i.e. d(heat)/d(rho).
    """

    rho = np.reshape(rho, (nz, ny, nx))
    n = nx * ny * nz
    # gradient matrices
    gx = diags([-1, 1], [0, 1], shape=(nx - 1, nx), format="csr")
    gy = diags([-1, 1], [0, 1], shape=(ny - 1, ny), format="csr")
    gz = diags([-1, 1], [0, 1], shape=(nz, nz), format="csr")
    # -div matrices
    dx, dy, dz = gx.copy().transpose(), gy.copy().transpose(), gz.copy().transpose()
    # 1D -> 3D
    Ix, Iy, Iz = eye(nx), eye(ny), eye(nz)
    gx, gy, gz = kron(Iz, kron(Iy, gx)), kron(Iz, kron(gy, Ix)), kron(gz, kron(Iy, Ix))
    dx, dy, dz = kron(Iz, kron(Iy, dx)), kron(Iz, kron(dy, Ix)), kron(dz, kron(Iy, Ix))

    # harmonic mean as mid-point conductivity and derivatives
    f = lambda x, y: 2 * x * y / (x + y)
    fx = lambda x, y: 2 * (y / (x + y)) ** 2
    fy = lambda x, y: 2 * (x / (x + y)) ** 2

    cond = cond_v + (cond_s - cond_v) * rho
    condx = [
        f(cond[k, j, i], cond[k, j, i + 1])
        for k in range(nz)
        for j in range(ny)
        for i in range(nx - 1)
    ]
    condy = [
        f(cond[k, j, i], cond[k, j + 1, i])
        for k in range(nz)
        for j in range(ny - 1)
        for i in range(nx)
    ]
    condz = [
        f(cond[k, j, i], cond[k + 1, j, i])
        for k in range(nz - 1)
        for j in range(ny)
        for i in range(nx)
    ]
    condz = np.append(condz, [cond_s] * nx * ny)  # conductivity at Dirichlet boundary
    condx, condy, condz = diags(condx, 0), diags(condy, 0), diags(condz, 0)

    src = src_v + (src_s - src_v) * rho.flatten()
    eq = dx @ condx @ gx + dy @ condy @ gy + dz @ condz @ gz
    solver = solvers[solver_option]
    if solver == spsolve:
        T = solver(eq, src)
    else:
        T, info = (solver(eq, src)).reshape(1, -1)

    heat_func = lambda x: npa.sum(x**p) ** (1 / p) / thresh
    # Instead of constraining heat - threshold <= 0, we are constraining heat/threshold - 1 <= 0
    heat_constraint = heat_func(T) - 1
    if not need_grad:
        return heat_constraint

    dgdx = grad(heat_func)(T)
    if solver == spsolve:
        aT = (solver(eq, dgdx)).reshape(1, -1)
    else:
        aT, sinfo = (solver(eq, dgdx)).reshape(1, -1)

    # conductivity matrix derivative w.r.t 1st and 2nd argument of each entries (e.g. fx and fy)
    dcondx_r = [
        k * (nx - 1) * ny + j * (nx - 1) + i
        for k in range(nz)
        for j in range(ny)
        for i in range(nx - 1)
    ]
    dcondx_c1 = [
        k * nx * ny + j * nx + i
        for k in range(nz)
        for j in range(ny)
        for i in range(nx - 1)
    ]
    dcondx_d1 = [
        fx(cond[k, j, i], cond[k, j, i + 1])
        for k in range(nz)
        for j in range(ny)
        for i in range(nx - 1)
    ]
    dcondx1 = csc_matrix(
        (dcondx_d1, (dcondx_r, dcondx_c1)), shape=(nz * ny * (nx - 1), nx * ny * nz)
    )

    dcondx_c2 = [
        k * nx * ny + j * nx + i + 1
        for k in range(nz)
        for j in range(ny)
        for i in range(nx - 1)
    ]
    dcondx_d2 = [
        fy(cond[k, j, i], cond[k, j, i + 1])
        for k in range(nz)
        for j in range(ny)
        for i in range(nx - 1)
    ]
    dcondx2 = csc_matrix(
        (dcondx_d2, (dcondx_r, dcondx_c2)), shape=(nz * ny * (nx - 1), nx * ny * nz)
    )

    dcondy_r = [
        k * nx * (ny - 1) + j * nx + i
        for k in range(nz)
        for j in range(ny - 1)
        for i in range(nx)
    ]
    dcondy_c1 = [
        k * nx * ny + j * nx + i
        for k in range(nz)
        for j in range(ny - 1)
        for i in range(nx)
    ]
    dcondy_d1 = [
        fx(cond[k, j, i], cond[k, j + 1, i])
        for k in range(nz)
        for j in range(ny - 1)
        for i in range(nx)
    ]
    dcondy1 = csc_matrix(
        (dcondy_d1, (dcondy_r, dcondy_c1)), shape=(nz * (ny - 1) * nx, nx * ny * nz)
    )

    dcondy_c2 = [
        k * nx * ny + (j + 1) * nx + i
        for k in range(nz)
        for j in range(ny - 1)
        for i in range(nx)
    ]
    dcondy_d2 = [
        fy(cond[k, j, i], cond[k, j + 1, i])
        for k in range(nz)
        for j in range(ny - 1)
        for i in range(nx)
    ]
    dcondy2 = csc_matrix(
        (dcondy_d2, (dcondy_r, dcondy_c2)), shape=(nz * (ny - 1) * nx, nx * ny * nz)
    )

    dcondz_r = [
        k * nx * ny + j * nx + i
        for k in range(nz - 1)
        for j in range(ny)
        for i in range(nx)
    ]
    dcondz_c1 = [
        k * nx * ny + j * nx + i
        for k in range(nz - 1)
        for j in range(ny)
        for i in range(nx)
    ]
    dcondz_d1 = [
        fx(cond[k, j, i], cond[k + 1, j, i])
        for k in range(nz - 1)
        for j in range(ny)
        for i in range(nx)
    ]
    dcondz1 = csc_matrix(
        (dcondz_d1, (dcondz_r, dcondz_c1)), shape=(nx * ny * nz, nx * ny * nz)
    )

    dcondz_c2 = [
        (k + 1) * nx * ny + j * nx + i
        for k in range(nz - 1)
        for j in range(ny)
        for i in range(nx)
    ]
    dcondz_d2 = [
        fy(cond[k, j, i], cond[k + 1, j, i])
        for k in range(nz - 1)
        for j in range(ny)
        for i in range(nx)
    ]
    dcondz2 = csc_matrix(
        (dcondz_d2, (dcondz_r, dcondz_c2)), shape=(nx * ny * nz, nx * ny * nz)
    )

    dcondx, dcondy, dcondz = dcondx1 + dcondx2, dcondy1 + dcondy2, dcondz1 + dcondz2
    dAdp_x = (
        dz @ dcondz.multiply(gz @ T.reshape(-1, 1))
        + dy @ dcondy.multiply(gy @ T.reshape(-1, 1))
        + dx @ dcondx.multiply(gx @ T.reshape(-1, 1))
    )

    gradient = aT * (src_s - src_v) - (cond_s - cond_v) * (aT @ dAdp_x)
    return T, heat_constraint, gradient


# Finite difference gradient for debugging.
def cc_fd(
    rho,
    nx,
    ny,
    nz,
    cond_v=1,
    cond_s=1e6,
    src_v=0,
    src_s=1,
    solver_option=0,
    thresh=None,
    p=4,
    num_grad=6,
    db=1e-6,
):
    n = nx * ny * nz
    fdidx = np.random.choice(n, num_grad)
    fdgrad = []
    for k in fdidx:
        rho[k] += db
        fp = constraint_connectivity(
            rho.reshape((nz, ny, nx)),
            nx,
            ny,
            nz,
            cond_v,
            cond_s,
            src_v,
            src_s,
            solver_option,
            thresh,
            p,
            need_grad=False,
        )
        rho[k] -= 2 * db
        fm = constraint_connectivity(
            rho.reshape((nz, ny, nx)),
            nx,
            ny,
            nz,
            cond_v,
            cond_s,
            src_v,
            src_s,
            solver_option,
            thresh,
            p,
            need_grad=False,
        )
        fdgrad.append((fp - fm) / (2 * db))
        rho[k] += db
    return fdidx, fdgrad
