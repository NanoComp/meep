# -*- coding: utf-8 -*-

# polarization grating from C. Oh and M.J. Escuti, Optics Letters, Vol. 33, No. 20, pp. 2287-9, 2008
# note: reference uses z as the propagation direction and y as the out-of-plane direction; this script uses x and z, respectively

import meep as mp
import math
import argparse

resolution = 50        # pixels/μm

dpml = 1.0             # PML thickness
dsub = 1.0             # substrate thickness
dpad = 1.0             # padding thickness

k_point = mp.Vector3(0,0,0)

n_0 = 1.55
delta_n = 0.159
epsilon_diag = mp.Matrix(mp.Vector3(n_0**2,0,0),mp.Vector3(0,n_0**2,0),mp.Vector3(0,0,(n_0+delta_n)**2))

wvl = 0.54             # center wavelength
fcen = 1/wvl           # center frequency
df = 0.05*fcen         # frequency width

def pol_grating(d,ph,gp,nmode):
    sx = dpml+dsub+d+d+dpad+dpml
    sy = gp

    cell_size = mp.Vector3(sx,sy,0)
    pml_layers = [mp.PML(thickness=dpml,direction=mp.X)]

    # twist angle of nematic director; from equation 1b
    def phi(p):
        xx  = p.x-(-0.5*sx+dpml+dsub)
        if (xx >= 0) and (xx <= d):
            return math.pi*p.y/gp + ph*xx/d
        else:
            return math.pi*p.y/gp - ph*xx/d + 2*ph

    # return the anisotropic permittivity tensor for a uniaxial, twisted nematic liquid crystal
    def lc_mat(p):
        # rotation matrix for rotation around x axis
        Rx = mp.Matrix(mp.Vector3(1,0,0),mp.Vector3(0,math.cos(phi(p)),math.sin(phi(p))),mp.Vector3(0,-math.sin(phi(p)),math.cos(phi(p))))
        lc_epsilon = Rx * epsilon_diag * Rx.transpose()
        lc_epsilon_diag = mp.Vector3(lc_epsilon[0].x,lc_epsilon[1].y,lc_epsilon[2].z)
        lc_epsilon_offdiag = mp.Vector3(lc_epsilon[1].x,lc_epsilon[2].x,lc_epsilon[2].y)
        return mp.Medium(epsilon_diag=lc_epsilon_diag,epsilon_offdiag=lc_epsilon_offdiag)

    geometry = [mp.Block(center=mp.Vector3(-0.5*sx+0.5*(dpml+dsub)),size=mp.Vector3(dpml+dsub,mp.inf,mp.inf),material=mp.Medium(index=n_0)),
                mp.Block(center=mp.Vector3(-0.5*sx+dpml+dsub+d),size=mp.Vector3(2*d,mp.inf,mp.inf),material=lc_mat)]

    # linear-polarized planewave pulse source
    src_pt = mp.Vector3(-0.5*sx+dpml+0.3*dsub,0,0)
    sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df), component=mp.Ez, center=src_pt, size=mp.Vector3(0,sy,0)),
               mp.Source(mp.GaussianSource(fcen,fwidth=df), component=mp.Ey, center=src_pt, size=mp.Vector3(0,sy,0))]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=pml_layers,
                        k_point=k_point,
                        sources=sources,
                        default_material=mp.Medium(index=n_0))

    refl_pt = mp.Vector3(-0.5*sx+dpml+0.5*dsub,0,0)
    refl_flux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=refl_pt, size=mp.Vector3(0,sy,0)))

    sim.run(until_after_sources=100)

    input_flux = mp.get_fluxes(refl_flux)
    input_flux_data = sim.get_flux_data(refl_flux)

    sim.reset_meep()

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=pml_layers,
                        k_point=k_point,
                        sources=sources,
                        geometry=geometry)

    refl_flux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=refl_pt, size=mp.Vector3(0,sy,0)))
    sim.load_minus_flux_data(refl_flux,input_flux_data)

    tran_pt = mp.Vector3(0.5*sx-dpml-0.5*dpad,0,0)
    tran_flux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=tran_pt, size=mp.Vector3(0,sy,0)))

    sim.run(until_after_sources=300)

    res1 = sim.get_eigenmode_coefficients(tran_flux, range(1,nmode+1), eig_parity=mp.ODD_Z+mp.EVEN_Y)
    res2 = sim.get_eigenmode_coefficients(tran_flux, range(1,nmode+1), eig_parity=mp.EVEN_Z+mp.ODD_Y)
    angles = [math.degrees(math.acos(kdom.x/fcen)) for kdom in res1.kdom]

    return input_flux[0], angles, res1.alpha[:,0,0], res2.alpha[:,0,0];


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-dd', type=float, default=1.7, help='chiral layer thickness (default: 1.7 μm)')
    parser.add_argument('-ph', type=float, default=70, help='chiral layer twist angle (default: 70°)')
    parser.add_argument('-gp', type=float, default=6.5, help='grating period (default: 6.5 μm)')
    parser.add_argument('-nmode', type=int, default=5, help='number of mode coefficients to compute (default: 5)')
    args = parser.parse_args()
    
    input_flux, angles, coeffs1, coeffs2 = pol_grating(args.dd,math.radians(args.ph),args.gp,args.nmode)

    tran1 = abs(coeffs1+1j*coeffs2)**2/input_flux
    tran2 = abs(coeffs1-1j*coeffs2)**2/input_flux
    print("tran:, 0, 0, {:.5f}".format(0.5*(tran1[0]+tran2[0])))
    for m in range(1,args.nmode):
        print("tran:, +{}, +{:.2f}, {:.5f}".format(m,angles[m],tran1[m]))
        print("tran:, -{}, -{:.2f}, {:.5f}".format(m,angles[m],tran2[m]))
