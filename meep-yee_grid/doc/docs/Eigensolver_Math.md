---
# Frequency-Domain Eigensolver
---

In these notes, we describe the mathematical background for the frequency-domain eigensolver algorithms in Meep (used by the `solve_eigfreq` function), which build on the [frequency-domain solver](Python_Tutorials/Frequency_Domain_Solver.md) algorithm described in the [Meep paper](http://doi.org/10.1016/j.cpc.2009.11.008).  See also [these notes](https://github.com/mitmath/18369/blob/9d61d1731af4ad32ea924e1af57b89e7e6a6c488/notes/time-evolution.pdf) for more background on the notation for Maxwell's equations used below.

Maxwell Equations in the Frequency Domain
-----------------------------------------

The Meep frequency-domain solver (at a frequency $\omega$) is essentially solving the linear system of equations $\hat{M}\psi=\xi$, where:

$$\hat{M}(\omega)\psi=\overbrace{\left[\underbrace{\left(\begin{array}{cc}
 & \nabla\times\\
-\nabla\times
\end{array}\right)}_{\hat{C}}+i\omega(1+\hat{\chi})\right]}^{\hat{M}(\omega)}\underbrace{\left(\begin{array}{c}
\mathbf{E}\\
\mathbf{H}
\end{array}\right)}_{\psi}=\underbrace{\left(\begin{array}{c}
\mathbf{J}\\
\mathbf{K}
\end{array}\right)}_{\xi}.$$

That is, $\psi$ is the six-component field state, $\hat{\chi}(\omega)$ is the $6\times6$ susceptibility at $\omega$, $\xi$ is the six-component current, and $\hat{C}$ is the $6\times6$ "curl" operator. (For these notes "natural" units are used in which $\varepsilon_{0}=\mu_{0}=1$.) Furthermore, in most physical circumstances the matrix $\hat{\chi}$ block-diagonalizes, and $1+\hat{\chi}$ is related to the electric permittivity $\varepsilon$ and the magnetic permeability $\mu$:

$$1+\hat{\chi}(\omega)=\left(\begin{array}{cc}
\varepsilon(\omega)\\
 & \mu(\omega)
\end{array}\right).$$

Equivalently, we can write things in terms of the $\mathbf{D}$ and $\mathbf{B}$ fields
$$\Psi=\left(\begin{array}{c}
\mathbf{D}\\
\mathbf{B}
\end{array}\right)=(1+\hat{\chi})\psi,$$

yielding

$$\overbrace{\left[\hat{C}(1+\hat{\chi})^{-1}+i\omega\right]}^{\hat{A}(\omega)}\Psi=\underbrace{(1+\hat{\chi})^{-1}\left(\begin{array}{c}
\mathbf{J}\\
\mathbf{K}
\end{array}\right)}_{\Xi}.$$

In fact, this equation $\hat{A}\Psi=\Xi$ is precisely the linear system of equations that Meep is actually solving internally.

Maxwell Eigenproblems
---------------------

The eigenvalue problem is simply $$\hat{M}(\omega)\psi=0\Longleftrightarrow\hat{A}(\omega)\Psi=0$$ i.e. to find a (complex resonant) frequency $\omega$ for which $\hat{M}$ is singular, and a corresponding eigenvector $\psi$ in the nullspace of $\hat{M}(\omega)$. If we have *non-dispersive* materials, those where $1+\hat{\chi}$ is *independent* of $\omega$, then this is a linear generalized eigenproblem $$\hat{C}\psi=-i\omega(1+\hat{\chi})\psi$$ or equivalently a linear eigenproblem $$i(1+\hat{\chi})^{-1}\hat{C}\psi=\omega\psi\Longleftrightarrow i\hat{C}(1+\hat{\chi})^{-1}\Psi=\omega\Psi.$$ For lossless transparent materials, i.e. real $\hat{\chi}>-1$, these problems are Hermitian under an inner product weighted by $1+\hat{\chi}$ (for $\psi$) or $(1+\hat{\chi})^{-1}$ (for $\Psi$), which leads to real $\omega$. More generally, for *dispersive* ($\omega$-dependent) materials, $\hat{M}(\omega)\psi=0$ or $\hat{A}(\omega)\Psi=0$ is a "*nonlinear* eigenvalue problem."

Iterative Eigenvalue Algorithms
-------------------------------

There are many algorithms for linear and nonlinear eigenvalue problems, but let us focus on the case where we have a **good initial guess** $\omega_{0}$ for the desired eigenvalue $\omega$. That is suppose we want the *closest eigenvalue* to $\omega_{0}$, and that there is a single eigenvalue *much closer* to $\omega_{0}$ than any other eigenvalue. (In a time-domain solver like Meep, we get good estimates for many eigenvalues simultaneously by signal processing analyses of the response to a short pulse input.)

In fact, suppose that $|\omega-\omega_{0}|\ll\omega_{0}$ is so small that we can approximate $\hat{\chi}(\omega_{0})\approx\hat{\chi}(\omega)$, allowing us to *neglect material dispersion* when computing $\omega$. In this case, we can use the standard [shift-and-invert power method](https://en.wikipedia.org/wiki/Inverse_iteration) to repeatedly solve $$\left[i(1-\hat{\chi})^{-1}\hat{C}-\omega_{0}\right]\psi_{n+1}=\psi_{n}\Longleftrightarrow\hat{M}(\omega_{0})\psi_{n+1}=-i(1+\hat{\chi})\psi_{n}\Longleftrightarrow\hat{A}(\omega_{0})\Psi_{n+1}=-i\Psi_{n}$$ with some arbitrary $\psi_{0}$ (e.g. random or a point source). This is, in fact, just a Maxwell solve, where the "current source" is $$\xi=-i(1+\hat{\chi})\psi_{n}=-i\Psi_{n},$$ i.e. the $\mathbf{D}$ and $\mathbf{B}$ fields of the previous solve. The $-i$ factor is essentially irrelevant, since we can scale eigenfunctions arbitrarily, and in fact one ordinarily wants to renormalize $\Psi_{n}\to\Psi_{n}/\Vert\Psi_{n}\Vert$ on each power iteration to prevent the iterations $\Psi_{n}$ from blowing up (or decaying to zero).

Equivalently, we are solving the shifted eigenproblem $$\hat{A}(\omega_{0})\Psi=-i(\omega-\omega_{0})\Psi$$ whose eigenvalue is $-i(\omega-\omega_{0})$ instead of $\omega$.

Estimating the Eigenvalue
-------------------------

Given an estimated eigenvector $\Psi_{n}$, the typical way to estimate the corresponding eigenvalue $-i(\omega_{n}-\omega_{0})$ is to compute a [Rayleigh quotient](https://en.wikipedia.org/wiki/Rayleigh_quotient)

$$-i(\omega_{n}-\omega_{0})=\frac{\langle\Psi_{n},\hat{A}(\omega_{0})\Psi_{n}\rangle}{\langle\Psi_{n},\Psi_{n}\rangle}$$

using some inner product $\langle\cdot,\cdot\rangle$. For a general non-normal $\hat{A}$ where we have arbitrary complex eigenvalues, it doesn't matter too much which inner product we choose, e.g. the obvious inner product $\langle\Psi,\Phi\rangle=\int\Psi^{*}\Phi$ is fine.

In the case of lossless media (Hermitian positive-definite $1+\hat{\chi}$) with real $\omega$, the accuracy of $\omega_{n}$ can be improved by using an inner product where $i\hat{A}$ is Hermitian, i.e. $\langle\Psi,\Phi\rangle_{\hat{\chi}}=\int\Psi^{*}(1+\hat{\chi})^{-1}\Phi$, which implies that $\langle\Psi,\Psi\rangle_{\hat{\chi}}=\int\Psi^{*}\psi=\int(\mathbf{D}^{*}\mathbf{E}+\mathbf{B}^{*}\mathbf{H})$ which corresponds physically to electromagnetic energy. Doing this essentially squares the error (i.e. it doubles the number of digits in the eigenvalue estimate) because eigenvalues of Hermitian operators are extrema of the Rayleigh quotient. Unfortunately, if the medium is not lossless, you can run into problems because this $\langle\Psi,\Phi\rangle_{\hat{\chi}}$ is not an inner product, and one can even have $\langle\Psi,\Psi\rangle_{\hat{\chi}}=0$ at exceptional points. Since the main utility of computing eigenvalues in Meep is arguably for computing resonance modes of non-Hermitian problems (since most Hermitian cases can be handled more efficiently in MPB), we should probably just stick with the $\langle\Psi,\Phi\rangle=\int\Psi^{*}\Phi$ inner product.

Correcting for Time Discretization
----------------------------------

Meep does not compute $\frac{\partial}{\partial t}$ exactly, of course—it uses a finite-difference approximation $\hat{D}$: $$\left.\frac{\partial\Phi}{\partial t}\right|_{\Delta t/2}\approx\hat{D}\Phi=\frac{\Phi(\Delta t)-\Phi(0)}{\Delta t}.$$ So, whereas a time-harmonic field $\Phi(t)=e^{-i\omega t}\Phi(0)$ would have $\frac{\partial\Phi}{\partial t}=-i\omega\Phi$, we instead have $$\hat{D}\Phi=\underbrace{\frac{e^{-i\omega\Delta t}-1}{\Delta t}}_{-i\hat{\omega}}\Phi.$$ Note that $-i\hat{\omega}=-i\omega+O(\Delta t)$, so that the two agree for $\Delta t\to0$.

In all of the analyses above, we simply replace $\omega$ with $\hat{\omega}$ (for $\omega$, $\omega_{0}$, and $\omega_{n}$), and everything carries through in the same way. At the end of an eigenvalue calculation, we compute $\omega$ from $\hat{\omega}$ using the formula $$\omega=\frac{\log(1-i\hat{\omega}\Delta t)}{-i\Delta t}.$$

Further Improvements
--------------------

There are many ways to potentially improve a numerical eigensolver beyond the simple shift-and-invert power method describe above. For example, the most common technique would be to plug the same $\hat{A}(\omega_{0})^{-1}$ solves into an [Arnoldi iteration](https://en.wikipedia.org/wiki/Arnoldi_iteration), e.g. as implemented by a library like [ARPACK](https://www.caam.rice.edu/software/ARPACK/). An advantage of Arnoldi iterations, beyond accelerated convergence (especially if our shift estimate $\omega_{0}$ is not so accurate) is that it can compute multiple eigenvalues simulaneously (albeit with increased computational expense).
