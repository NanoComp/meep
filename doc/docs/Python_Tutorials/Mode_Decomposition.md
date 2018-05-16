---
# Mode Decomposition
---

This tutorial demonstrates Meep's mode-decomposition feature which is used to separate a given field profile into a superposition of harmonic basis modes. The example involves computing the reflection coefficient of the fundamental mode of a linear waveguide taper as shown in the schematic below. We will verify that the scaling of the reflection coefficient with the taper length is quadratic, consistent with analytical results from [Optics Express, Vol. 16, pp. 11376-92, 2008](http://www.opticsinfobase.org/abstract.cfm?URI=oe-16-15-11376).

<center>
![](../images/waveguide-taper.png)
</center>

The 2d structure consists of a single-mode waveguide of width 1 μm (w1) with ε=12 in vacuum coupled to a second waveguide of width 2 μm (w2) via a linearly-sloped taper of length Lt. PML absorbing boundaries surround the computational cell. An eigenmode source is used to launch the fundamental mode. There is an eigenmode-expansion monitor placed at the midpoint of the first waveguide. This is a line monitor which extends beyond the waveguide in order to capture the entire mode profile including its evanescent tails. The Fourier-transformed fields along the line monitor are used to compute the basis coefficients of the harmonic modes which are obtained separately via the eigenmode solver MPB. The details of this calculation are provided in [Mode Decomposition](../Mode_Decomposition). An alternative method for computing the reflection coefficient involves [calculating the Poynting flux](Basics/#transmission-spectrum-of-a-waveguide-bend). This approach requires a separate normalization run to compute the incident fields due to the source and does not provide phase information. The mode-decomposition feature uses only a single run to compute the complex reflection coefficient. This enables the computation of [S parameters](https://en.wikipedia.org/wiki/Scattering_parameters). The simulation script is shown below.

```py
import meep as mp
import math
import argparse
import numpy as np

def main(args):

    resolution = args.res

    w1 = 1            # width of waveguide 1
    w2 = args.w2      # width of waveguide 2
    Lw = 10           # length of waveguide 1 and 2
    Lt = args.Lt      # taper length

    Si = mp.Medium(epsilon=12.0)
    
    dair = 3.0
    dpml = 3.0
    
    sx = dpml + Lw + Lt + Lw + dpml
    sy = dpml + dair + w2 + dair + dpml
    cell_size = mp.Vector3(sx, sy, 0)

    geometry = [ mp.Block(material=Si, center=mp.Vector3(0,0,0), size=mp.Vector3(mp.inf,w1,mp.inf)), 
                 mp.Block(material=Si, center=mp.Vector3(0.5*sx-0.5*(Lt+Lw+dpml),0,0), size=mp.Vector3(Lt+Lw+dpml,w2,mp.inf)) ]

    # form linear taper

    hh = w2
    ww = 2*Lt
    
    # taper angle (CCW, relative to +X axis)
    rot_theta = math.atan(0.5*(w2-w1)/Lt)
    
    pvec = mp.Vector3(-0.5*sx+dpml+Lw,0.5*w1,0)
    cvec = mp.Vector3(-0.5*sx+dpml+Lw+0.5*ww,0.5*hh+0.5*w1,0)
    rvec = cvec-pvec
    rrvec = rvec.rotate(mp.Vector3(0,0,1), rot_theta)

    geometry.append(mp.Block(material=mp.air, center=pvec+rrvec, size=mp.Vector3(ww,hh,mp.inf),
                             e1=mp.Vector3(1,0,0).rotate(mp.Vector3(0,0,1),rot_theta),
                             e2=mp.Vector3(0,1,0).rotate(mp.Vector3(0,0,1),rot_theta),
                             e3=mp.Vector3(0,0,1)))

    pvec = mp.Vector3(-0.5*sx+dpml+Lw,-0.5*w1,0)
    cvec = mp.Vector3(-0.5*sx+dpml+Lw+0.5*ww,-(0.5*hh+0.5*w1),0)
    rvec = cvec-pvec
    rrvec = rvec.rotate(mp.Vector3(0,0,1),-rot_theta)

    geometry.append(mp.Block(material=mp.air, center=pvec+rrvec, size=mp.Vector3(ww,hh,mp.inf),
                             e1=mp.Vector3(1,0,0).rotate(mp.Vector3(0,0,1),-rot_theta),
                             e2=mp.Vector3(0,1,0).rotate(mp.Vector3(0,0,1),-rot_theta),
                             e3=mp.Vector3(0,0,1)))

    boundary_layers = [ mp.PML(dpml) ]

    # mode frequency
    fcen = 0.15
    
    sources = [ mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=0.5*fcen),
                                   component=mp.Ez,
                                   size=mp.Vector3(0,sy-2*dpml,0),
                                   center=mp.Vector3(-0.5*sx+dpml,0,0),
                                   eig_match_freq=True,
                                   eig_parity=mp.ODD_Z,
                                   eig_kpoint=mp.Vector3(0.4,0,0),
                                   eig_resolution=32) ]
    
    symmetries = [ mp.Mirror(mp.Y,+1) ]
    
    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=boundary_layers,
                        geometry=geometry,
                        sources=sources,
                        symmetries=symmetries)

    xm = -0.5*sx + dpml + 0.5*Lw; # x-coordinate of monitor
    mflux = sim.add_eigenmode(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(xm,0), size=mp.Vector3(0,sy-2*dpml)));
        
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(xm,0), 1e-10))

    bands = np.array([1],dtype=np.int32) # indices of modes for which to compute expansion coefficients
    num_bands = 1;
    alpha = np.zeros(2*num_bands,dtype=np.complex128) # preallocate array to store coefficients
    vgrp = np.zeros(num_bands,dtype=np.float64)       # also store mode group velocities
    mvol = mp.volume(mp.vec(xm,-0.5*sy+dpml),mp.vec(xm,+0.5*sy-dpml))
    sim.fields.get_eigenmode_coefficients(mflux, mp.X, mvol, bands, alpha, vgrp)

    alpha0Plus  = alpha[2*0 + 0]; # coefficient of forward-traveling fundamental mode
    alpha0Minus = alpha[2*0 + 1]; # coefficient of backward-traveling fundamental mode

    print("refl:, {}, {:.8f}".format(Lt, abs(alpha0Minus)**2))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-Lt',  type=float, default=3.0, help='taper length (default: 3.0)')
    parser.add_argument('-w2',  type=float, default=2.0, help='width of outgoing waveguide (default: 2.0)')
    parser.add_argument('-res', type=int, default=20, help='resolution (default: 20)')
    args = parser.parse_args()
    main(args)
```

To investigate the scaling, we compute the reflection coefficient for a range of taper lengths chosen on a logarithmic scale (i.e., 1, 2, 4, ..., 64). The results are shown below in blue. For reference, a quadratic scaling is shown in black. Consistent with analytical results, the reflection coefficient of the linear waveguide taper decreases quadratically with the taper length.

<center>
![](../images/refl_coeff_vs_taper_length.png)
</center>
