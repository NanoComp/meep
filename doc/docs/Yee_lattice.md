---
# Yee lattice
---


![400px|center|thumb|Illustration of Yee lattice in 3d for a single grid voxel.](images/Yee-cube.png)

 In order to discretize Maxwell's equations with second-order accuracy (for homogeneous regions where there no discontinuous material boundaries), FDTD methods *store different field components for different grid locations*. This discretization is known as a **Yee lattice**.

The form of the Yee lattice in 3d is shown in the illustration above for a single cubic grid voxel ($\Delta x \times \Delta x \times \Delta x$). The basic idea is that the three components of **E** are stored for the *edges* of the cube in the corresponding directions, while the components of **H** are stored for the *faces* of the cube.

More precisely, let a coordinate $(i,j,k)$ in the grid correspond to:

$$\textbf{x} = (i \hat\textbf{e}_1 + j \hat\textbf{e}_2 + k \hat\textbf{e}_3) \Delta x$$,

where $\hat\textbf{e}_k$ denotes the unit vector in the *k*-th coordinate direction. Then, the $\ell$<sup>th</sup> component of $\textbf{E}$ or $\textbf{D}$ (or $\textbf{P}$) is stored for the locations

$$(i,j,k)+ \frac{1}{2} \hat\textbf{e}_\ell  \Delta x$$.

The $\ell$<sup>th</sup> component of $\textbf{H}$, on the other hand, is stored for the locations

$$(i+\frac{1}{2},j+\frac{1}{2},k+\frac{1}{2})-\frac{1}{2} \hat\textbf{e}_\ell  \Delta x$$.

In two dimensions, the idea is similar except that we set $\hat\textbf{e}_3=0$. The 2d Yee lattice for the <i>P</i>-polarization (**E** in the *xy* plane and **H** in the *z* direction) is shown in the figure below. 
![thumb|center|250px|Yee lattice in 2d for the TE polarization.](images/Yee-te.png)



The consequence of the Yee lattice is that, whenever you need to compare or combine different field components, e.g. to find the energy density $(\textbf{E}^* \cdot \textbf{D} + |\textbf{H}|^2)/2$ or the flux $\textrm{Re}\, \textbf{E}^* \times \textbf{H}$, then the components need to be **interpolated** to some common point (in order to remain second-order accurate). Meep automatically does this interpolation for you wherever necessary—in particular, whenever you compute energy density or flux, or whenever you output a field to a file, it is stored for the locations $(i+0.5,j+0.5,k+0.5)$: the centers of each grid voxel.

[Category:Meep](Meep.md)
