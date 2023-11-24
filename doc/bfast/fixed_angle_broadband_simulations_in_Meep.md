---
# Fixed Angle Broadband Simulations in Meep
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
the equations will be modified.  B. Liang et al applied this method by
splitting the electric and magnetic fields into multiple terms. This is
different from the UPML formulation in MEEP, which allows simulation of
more complex media, and so the method of B. Liang et al has to be adapted.

Boundary conditions
===================

The fields from section 3 of *Notes on the UPML implementation in MEEP*
[3] are first redefined as:
```math
\text{field}'(x,y,z) = \text{field}(x,y,z)e^{-i(k_{x}x+k_{y}y)}, \tag{1}
```
where $k_{x}$ and $k_{y}$ are the wave vector components in the x and y
directions. This is for a structure which is periodic in these
directions. Taking the electric field $E$ as an example, the new
boundary condition can be expressed as
```math
E'(x+a,y+b,z) = E(x+a,y+b,z)e^{-i(k_{x}(x+a)+k_{y}(y+b))} \tag{2}
```
where a is the length of the unit cell in the x direction and b in the y direction.
Substituting in the original Bloch periodic boundary conditions gives
```math
E'(x+a,y+b,z) = E(x,y,z)e^{i(k_{x}a+k_{y}b)}e^{-i(k_{x}(x+a)+k_{y}(y+b))}. \tag{3}
```
Cancelling the $a$ and $b$ terms gives
```math
E'(x+a,y+b,z) =E(x,y,z)e^{-i(k_{x}x+k_{y}y)}=E'(x,y,z), \tag{4}
```
and so the boundary conditions are now periodic.

Formulation
===========

Equation (5) from section 3 of *Notes on the UPML implementation in
MEEP* [3] is
```math
\vec{K} = \nabla \times \vec{H}=-i\omega(1+\frac{i\sigma_{D}}{\omega})\vec{C}, \tag{5}
```
where $\vec{H}$ is the magnetic field, $\sigma_{D}$ the conductivity and
$\vec{C}$ an auxiliary field. When the magnetic field is redefined, the
curl of a product must be carried out:
```math
\nabla\times \vec{H'} = \nabla\times (\vec{H} e^{-i(k_{x}x+k_{y}y)}) \tag{6}
```
so,
```math
\nabla\times \vec{H'} = e^{-i(k_{x}x+k_{y}y)} \nabla\times \vec{H} + \begin{pmatrix} -ik_{x} \\-ik_{y}\\0 \end{pmatrix} \ \times \vec{H'} \tag{7}
```
where the complex exponential in the second term has been absorbed by
$\vec{H'}$. Substituting in equation (5) gives
```math
\nabla\times \vec{H'} = \vec{K'} = -i\omega (1+\frac{i\sigma_{D}}{\omega}) \vec{C'} + \begin{pmatrix} -ik_{x} \\-ik_{y}\\0 \end{pmatrix} \ \times \vec{H'}. \tag{8}
```
From here on in, the prime notation can be dropped since this applies to
all fields. By introducing a new field $\vec{F}$, equation
(8) can be written as
```math
\vec{K} = -i\omega(1+\frac{i\sigma_{D}}{\omega})\vec{C} - i\omega\vec{F}. \tag{9}
```
This new field satisfies the equation:
```math
\vec{F} = \vec{\bar{k}}\times\vec{H}, \tag{10}
```
where
```math
\vec{\bar{k}} =\frac{1}{\omega}\begin{pmatrix} k_{x} \\ k_{y}\\0 \end{pmatrix} \ = n \begin{pmatrix} \sin{\theta}\cos{\phi} \\ \sin{\theta}\sin{\phi}\\0 \end{pmatrix} \ \tag{11}
```
and so $\vec{\bar{k}}$ is the wave vector with its frequency dependence
removed. $\theta$ and $\phi$ are the propagating direction angles, $n$ is the refractive index of the source medium and c,
the speed of light is taken to be 1. Therefore by defining $\vec{F}$,
the angle of the incident wave is fixed. Equation
(10) can be discretized as:
```math
\vec{F}^{n+1}=2\bar{\vec{k}}\times\vec{H}^{n+0.5} -\vec{F}^{n} . \tag{12}
```
Transforming equation (9)
to the time domain gives:
```math
\vec{K} = \frac{\partial \vec{C}}{\partial t}+\sigma_{D}\vec{C}+\frac{\partial \vec{F}}{\partial t} . \tag{13}
```
This can be discretized as:
```math
\vec{K}^{n+0.5}=\frac{\vec{C}^{n+1}-\vec{C}^n}{\Delta t}+\sigma_{D}\frac{\vec{C}^{n+1}+\vec{C}^n}{2} + \frac{\vec{F}^{n+1}-\vec{F}^{n}}{\Delta t} \tag{14}
```
and then solved to update the value of $\vec{C}$ using:
```math
\vec{C}^{n+1}=(1+\frac{\sigma_{D} \Delta t}{2})^{-1} [(1-\frac{\sigma_D \Delta t}{2}) \vec{C}^n+\Delta t\vec{K}^{n+0.5}+\vec{F}^{n}-\vec{F}^{n+1}] . \tag{15}
```
All other equations are unaffected by these changes.

A new field must be introduced because $\vec{H}$ is defined at
$n+\frac{1}{2}$ timesteps whereas $\vec{C}$ is defined at $n$ timesteps,
where $n$ is an integer. As a result, if the derivative in $\vec{F}$ in
equation (14) was replaced with
```math
\vec{\bar{k}}\times(\frac{\vec{H}^{n+0.5}-\vec{H}^{n-0.5}}{\Delta t}), \tag{16}
```
only first order accuracy would be achieved, since this is a backward
difference scheme. To achieve second order accuracy would require
$\vec{H}^{n+1.5}$ to be known.

Application to the source code
=========

$\vec{C}$ is first time-stepped using the original function in the MEEP source code, i.e. using
```math
\vec{C}^{n+1}=(1+\frac{\sigma_{D} \Delta t}{2})^{-1} [(1-\frac{\sigma_D \Delta t}{2}) \vec{C}^n+\Delta t\vec{K}^{n+0.5}]. \tag{17}
```
When fixed angle broadband simulations are required, a new function is called which recovers equation (15) using
```math
\vec{C'}^{n+1} = (1+\frac{\sigma_{D} \Delta t}{2})^{-1} (\vec{F}^{n}-\vec{F}^{n+1}) + \vec{C}^{n+1}. \tag{18}
```
Similarly, all other fields which depend on the value of $\vec{C}$ have new terms added to them, multiplied by the relevant conductivity. Therefore the original functions in MEEP remain unchanged.

Stability
=========

As the incident angle increases, the maximum possible $\Delta t$ value
decreases, following the formula:
```math
\frac{c\Delta t}{\Delta x} \leq \frac{(1-\sin(\theta))}{\sqrt{D}} \tag{17}
```
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
