"""
Connectivity constraint inspired by https://link.springer.com/article/10.1007/s00158-016-1459-5
Solve and find adjoint gradients for [-div (k grad) + alpha^2 (1-rho)k + alpha0^2]T=0.
BC: Dirichlet T0 on last row of rho, 0 outside the first row, and Neumann on sides.
Mo Chen <mochen@mit.edu>
"""

import numpy as np
from scipy.sparse.linalg import *
from scipy.sparse import kron, diags, csr_matrix, eye, csc_matrix, lil_matrix
from scipy.linalg import norm


class ConnectivityConstraint(object):
    def __init__(self, nx, ny, k0=1000, T0=1000, zeta=0, sp_solver=spsolve, alpha=None, alpha0=None, thresh=0.9):
        #zeta is to prevent singularity when damping is zero; with damping, zeta should be zero
        self.nx, self.ny = nx, ny
        self.n = nx*ny
        self.m = nx*(ny-1)
        self.solver = sp_solver
        self.k0, self.T0, self.zeta = k0, T0, zeta
        self.thresh = thresh

        #default alpha and alpha0
        if alpha != None:
            self.alpha=alpha
        else:
            self.alpha = 0.1*min(1/nx, 1/ny)
        if alpha0 != None:
            self.alpha0 = alpha0
        else:
            self.alpha0 = -np.log(thresh)/min(nx, ny)


    def forward(self, rho_vector):
        self.rho_vector = rho_vector
        # gradient and -div operator for x and y directions
        gx = diags([-1,1], [0,1], shape=(self.nx-1, self.nx), format='csr')
        dx = gx.copy().transpose()
        gy = diags([1,-1], [0, -1], shape=(self.ny, self.ny), format='csr')
        dy = diags([1,-1], [0, 1], shape=(self.ny-1, self.ny), format='csr')

        # kron product for 2D
        Ix, Iy = eye(self.nx), eye(self.ny-1)
        self.gx, self.gy = kron(Iy, gx), kron(gy, Ix)
        self.dx, self.dy = kron(Iy, dx), kron(dy, Ix)

        #conductivity based on rho
        rho_pad = np.reshape(rho_vector, (self.ny, self.nx))
        rhox = np.array([0.5*(rho_pad[j, i]+rho_pad[j, i+1]) for j in range(self.ny-1) for i in range(self.nx-1)])
        rhoy = np.array([0.5*(rho_pad[j, i]+rho_pad[j+1, i]) for j in range(self.ny-1) for i in range(self.nx)])
        rhoy = np.insert(rhoy, [0]*self.nx, 0)#0 outside first row
        kx, ky = diags((self.zeta+(1-self.zeta)*rhox)*self.k0), diags((self.zeta+(1-self.zeta)*rhoy)*self.k0)
        #matrices in x and y
        Lx, Ly = self.dx * kx * self.gx, self.dy * ky * self.gy

        # Dirichlet condition on the last row becomes term on the RHS
        By = csc_matrix(Ly)[:,-self.nx:]
        rhs = -self.T0*By.sum(axis=1)

        #LHS operator after moving the boundary term to the RHS
        eq = Ly[:, :-self.nx]+Lx
        #add damping
        damping = self.k0*self.alpha**2*diags(1-rho_vector[:-self.nx]) + diags([self.alpha0**2], shape=(self.m, self.m))
        self.A = eq + damping
        self.T = self.solver(csr_matrix(self.A), rhs)
        #exclude last row of rho and calculate weighted average of temperature
        self.rho_vec = rho_vector[:-self.nx]
        return (sum(self.T*self.rho_vec)+self.T0*sum(rho_vector[-self.nx:]))/sum(rho_vector)

    def adjoint(self):
        dg_dT = self.rho_vec/sum(self.rho_vector)
        return self.solver(csr_matrix(self.A.transpose()), dg_dT)
    def calculate_grad(self):
        dg_dp = np.zeros(self.n)
        dg_dp[:-self.nx] = (self.T*sum(self.rho_vector))/sum(self.rho_vector)**2
        dg_dp[-self.nx:] = (self.T0*np.ones(self.nx)*sum(self.rho_vector))/sum(self.rho_vector)**2
        dg_dp = dg_dp - (sum(self.T*self.rho_vec)+self.T0*sum(self.rho_vector[-self.nx:]))/sum(self.rho_vector)**2

        dAx = lil_matrix((self.m, self.n))
        gxT = np.reshape(self.gx * self.T, (-1,1))
        drhox = kron(eye(self.ny-1), diags([0.5,0.5], [0, 1], shape=[self.nx-1,self.nx]))
        dAx[:, :-self.nx] =  (1-self.zeta)*self.k0*lil_matrix(self.dx * drhox.multiply(gxT)) #element-wise product

        Ty = np.pad(self.T, (0, self.nx), 'constant', constant_values=self.T0)
        gyTy = np.reshape(self.gy * Ty, (-1,1))
        drhoy = diags([0.5,0.5], [0, -1], shape=[self.ny,self.ny], format="lil")
        drhoy[0,0]=0
        drhoy = kron(drhoy,eye(self.nx))
        dAy =  (1-self.zeta)*self.k0*self.dy * drhoy.multiply(gyTy)

        d_damping = self.k0*self.alpha**2*diags(-self.T, shape=(self.m, self.n))

        self.grad = dg_dp + self.adjoint().reshape(1, -1) * csr_matrix( - dAy - dAx - d_damping)
        return -self.grad[0]

    def __call__(self, rho_vector):
        Tmean = self.forward(rho_vector)
        grad = self.calculate_grad()
        return -Tmean+self.thresh*self.T0, grad
    def calculate_fd_grad(self, rho_vector, num, db=1e-4):
        fdidx = np.random.choice(self.n, num)
        fdgrad = []
        for k in fdidx:
            rho_vector[k]+=db
            fp = self.forward(rho_vector)
            rho_vector[k]-=2*db
            fm = self.forward(rho_vector)
            fdgrad.append((fm-fp)/(2*db))
            rho_vector[k]+=db
        return fdidx, fdgrad
    def calculate_all_fd_grad(self, rho_vector, db=1e-4):
        fdgrad = []
        for k in range(self.n):
            rho_vector[k]+=db
            fp = self.forward(rho_vector)
            rho_vector[k]-=2*db
            fm = self.forward(rho_vector)
            fdgrad.append((fp-fm)/(2*db))
            rho_vector[k]+=db
        return range(self.n), np.array(fdgrad)
