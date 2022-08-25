import numpy as np
from scipy.sparse.linalg import cg, spsolve
from scipy.sparse import kron, diags, csr_matrix, eye, csc_matrix
from autograd import numpy as npa
from autograd import grad


def constraint_connectivity(
    rho,
    nx,
    ny,
    nz,
    cond_v=1,
    cond_s=1e6,
    src_v=0,
    src_s=1,
    solver=spsolve,
    thresh=None,
    p=4,
    need_grad=True,
):
    rho = np.reshape(rho, (nz, ny, nx))
    n = nx * ny * nz

    if ny == 1:
        path = nx * nz / 2
        phi = 0.5 * path * path / cond_s  # estimate of warmest connected structure
    else:
        path = nx * ny * nz / 4
        phi = 0.5 * path * path / cond_s
    if not thresh:
        thresh = phi

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
    if solver == spsolve:
        T = solver(eq, src)
    else:
        T, info = (solver(eq, src)).reshape(1, -1)

    heat_func = lambda x: npa.sum(x**p) ** (1 / p) / thresh
    heat = heat_func(T)
    if not need_grad:
        return heat / thresh - 1

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
    return T, heat / thresh - 1, gradient


def cc_fd(
    rho,
    nx,
    ny,
    nz,
    cond_v=1,
    cond_s=1e6,
    src_v=0,
    src_s=1,
    solver=spsolve,
    thresh=None,
    p=4,
    num_grad=6,
    db=1e-4,
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
            solver,
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
            solver,
            thresh,
            p,
            need_grad=False,
        )
        fdgrad.append((fp - fm) / (2 * db))
        rho[k] += db
    return fdidx, fdgrad
