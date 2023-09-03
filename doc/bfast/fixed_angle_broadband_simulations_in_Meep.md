---
author:
- 'Daniel Lloyd-Jones'
date: 28th July 2023
title: Fixed angle broadband simulations in MEEP
---

Introduction
============

Currently in MEEP, Bloch Periodic boundary conditions are implemented,
which fix the wave vector of an incident wave [1]. As a result, the
angle of an oblique incident wave becomes frequency dependent. Following
the procedure detailed by B. Liang et al [2], all fields can be
redefined so that the boundary conditions become periodic and the angle
of the incident wave can be fixed over a broad frequency spectrum. This
requires the addition of a new field. It is assumed that the reader is
already familiar with the UPML formulation in MEEP [3], from which
the equations will be modified.

Boundary conditions
===================

The fields from section 3 of *Notes on the UPML implementation in MEEP*
[3] are first redefined as:
$ \begin{equation}
\text{field}'(x,y,z) = \text{field}(x,y,z)e^{-i(k_{x}x+k_{y}y)},
\end{equation} $
where $k_{x}$ and $k_{y}$ are the wave vector components in the x and y
directions. This is for a structure which is periodic in these
directions. Taking the electric field $E$ as an example, the new
boundary condition can be expressed as
$ \begin{equation}
E'(x+a,y+b,z) = E(x+a,y+b,z)e^{-i(k_{x}(x+a)+k_{y}(y+b))}
\end{equation} $
where a is the length of the unit cell in the x direction and b in the y direction.
Substituting in the original Bloch periodic boundary conditions gives
$ \begin{equation}
E'(x+a,y+b,z) = E(x,y,z)e^{i(k_{x}a+k_{y}b)}e^{-i(k_{x}(x+a)+k_{y}(y+b))}.
\end{equation} $
Cancelling the $a$ and $b$ terms gives
$ \begin{equation}
E'(x+a,y+b,z) =E(x,y,z)e^{-i(k_{x}x+k_{y}y)}=E'(x,y,z),
\end{equation} $
and so the boundary conditions are now periodic.

Formulation
===========

Equation (5) from section 3 of *Notes on the UPML implementation in
MEEP* [3] is
$ \begin{equation}
\vec{K} = \nabla \times \vec{H}=-i\omega(1+\frac{i\sigma_{D}}{\omega})\vec{C},
\end{equation} $
where $\vec{H}$ is the magnetic field, $\sigma_{D}$ the conductivity and
$\vec{C}$ an auxiliary field. When the magnetic field is redefined, the
curl of a product must be carried out:
$ \begin{equation}
\nabla\times \vec{H'} = \nabla\times (\vec{H} e^{-i(k_{x}x+k_{y}y)})
\end{equation} $
so,
$ \begin{equation}
\nabla\times \vec{H'} = e^{-i(k_{x}x+k_{y}y)} \nabla\times \vec{H} + \begin{pmatrix} -ik_{x} \\-ik_{y}\\0 \end{pmatrix} \ \times \vec{H'}
\end{equation} $
where the complex exponential in the second term has been absorbed by
$\vec{H'}$. Substituting in equation (5) gives
$ \begin{equation}
\nabla\times \vec{H'} = \vec{K'} = -i\omega (1+\frac{i\sigma_{D}}{\omega}) \vec{C'} + \begin{pmatrix} -ik_{x} \\-ik_{y}\\0 \end{pmatrix} \ \times \vec{H'}.
\end{equation} $
From here on in, the prime notation can be dropped since this applies to
all fields. By introducing a new field $\vec{F}$, equation
(8) can
be written as
$ \begin{equation}
\vec{K} = -i\omega(1+\frac{i\sigma_{D}}{\omega})\vec{C} - i\omega\vec{F}.
\end{equation} $
This new field satisfies the equation:
$ \begin{equation}
\vec{F} = \vec{\bar{k}}\times\vec{H},
\end{equation} $
where
$ \begin{equation}
\vec{\bar{k}} =\frac{1}{\omega}\begin{pmatrix} k_{x} \\ k_{y}\\0 \end{pmatrix} \ =\begin{pmatrix} \sin{\theta}\cos{\phi} \\ \sin{\theta}\sin{\phi}\\0 \end{pmatrix} \
\end{equation} $
and so $\vec{\bar{k}}$ is the wave vector with its frequency dependence
removed. $\theta$ and $\phi$ are the propagating direction angles and c,
the speed of light is taken to be 1. Therefore by defining $\vec{F}$,
the angle of the incident wave is fixed. Equation
(10) can be discretized as:
$ \begin{equation}
\vec{F}^{n+1}=2\bar{\vec{k}}\times\vec{H}^{n+0.5} -\vec{F}^{n} .
\end{equation} $
Transforming equation (9)
to the time domain gives:
$ \begin{equation}
\vec{K} = \frac{\partial \vec{C}}{\partial t}+\sigma_{D}\vec{C}+\frac{\partial \vec{F}}{\partial t} .
\end{equation} $
This can be discretized as:
$ \begin{equation}
\vec{K}^{n+0.5}=\frac{\vec{C}^{n+1}-\vec{C}^n}{\Delta t}+\sigma_{D}\frac{\vec{C}^{n+1}+\vec{C}^n}{2} + \frac{\vec{F}^{n+1}-\vec{F}^{n}}{\Delta t}
\end{equation} $
and then solved to update the value of $\vec{C}$ using:
$ \begin{equation}
\vec{C}^{n+1}=(1+\frac{\sigma_{D} \Delta t}{2})^{-1} [(1-\frac{\sigma_D \Delta t}{2}) \vec{C}^n+\Delta t\vec{K}^{n+0.5}+\vec{F}^{n}-\vec{F}^{n+1}] .
\end{equation} $
All other equations are unaffected by these changes.

A new field must be introduced because $\vec{H}$ is defined at
$n+\frac{1}{2}$ timesteps whereas $\vec{C}$ is defined at $n$ timesteps,
where $n$ is an integer. As a result, if the derivative in $\vec{F}$ in
equation (14) was replaced with
$ \begin{equation}
\vec{\bar{k}}\times(\frac{\vec{H}^{n+0.5}-\vec{H}^{n-0.5}}{\Delta t}),
\end{equation} $
only first order accuracy would be achieved, since this is a backward
difference scheme. To achieve second order accuracy would require
$\vec{H}^{n+1.5}$ to be known.

Stability
=========

As the incident angle increases, the maximum possible $\Delta t$ value
decreases, following the formula:
$ \begin{equation}
\frac{c\Delta t}{\Delta x} \leq \frac{(1-sin(\theta))}{\sqrt{D}}
\end{equation} $
where D is the number of dimensions [2].

References
=========

[1] Taflove A., Oskooi A., Johnson S.. *Advances in FDTD Computational
Electrodynamics: Photonics and Nanotechnology*. Artech House, Inc.; 2013

[2] Liang B., Bai M., Ma H., Ou N., Miao J.. Wideband Analysis of Periodic
Structures at Oblique Incidence by Material Independent FDTD Algorithm.
*IEEE Transactions on Antennas and Propagation*, vol. 62, no. 1, pp.
354-360, Jan. 2014, doi: 10.1109/TAP.2013.2287896.

[3] Johnson S. *Notes on the UPML implementation in Meep*. Massachusetts
Institute of Technology. Posted August 17, 2009; updated March 10, 2010.
http://ab-initio.mit.edu/meep/pml-meep.pdf
