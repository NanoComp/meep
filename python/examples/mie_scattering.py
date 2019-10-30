import meep as mp
import numpy as np
import matplotlib.pyplot as plt

resolution = 35

wvl = 1.0
fcen = 1/wvl

dpml = 1.0
dair = 2.0

pml_layers = [mp.PML(thickness=dpml)]

symmetries = [mp.Mirror(mp.Y),
              mp.Mirror(mp.Z,phase=-1)]

def mie_scat(r):
    s = 2*(dpml+dair+r)
    cell_size = mp.Vector3(s,s,s)

    sources = [mp.Source(mp.GaussianSource(fcen,fwidth=0.2*fcen),
                         center=mp.Vector3(-0.5*s+dpml),
                         size=mp.Vector3(0,s,s),
                         component=mp.Ez)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=pml_layers,
                        sources=sources,
                        k_point=mp.Vector3(),
                        symmetries=symmetries)

    input = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(0.5*s-dpml),size=mp.Vector3(0,s-2*dpml,s-2*dpml)))

    box = sim.add_flux(fcen, 0, 1,
                       mp.FluxRegion(center=mp.Vector3(x=-r),size=mp.Vector3(0,2*r,2*r),weight=+1),
                       mp.FluxRegion(center=mp.Vector3(x=+r),size=mp.Vector3(0,2*r,2*r),weight=-1),
                       mp.FluxRegion(center=mp.Vector3(y=-r),size=mp.Vector3(2*r,0,2*r),weight=+1),
                       mp.FluxRegion(center=mp.Vector3(y=+r),size=mp.Vector3(2*r,0,2*r),weight=-1),
                       mp.FluxRegion(center=mp.Vector3(z=-r),size=mp.Vector3(2*r,2*r,0),weight=+1),
                       mp.FluxRegion(center=mp.Vector3(z=+r),size=mp.Vector3(2*r,2*r,0),weight=-1))

    sim.run(until_after_sources=20)

    input_flux = mp.get_fluxes(input)[0]
    box_flux_data = sim.get_flux_data(box)
    box_flux = mp.get_fluxes(box)[0]
    if mp.am_master():
        print("input-flux:, {:.6f}".format(input_flux))

    sim.reset_meep()

    geometry = [mp.Sphere(material=mp.metal,
                          center=mp.Vector3(),
                          radius=r)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=pml_layers,
                        sources=sources,
                        k_point=mp.Vector3(),
                        symmetries=symmetries,
                        geometry=geometry)

    box = sim.add_flux(fcen, 0, 1,
                       mp.FluxRegion(center=mp.Vector3(x=-r),size=mp.Vector3(0,2*r,2*r),weight=+1),
                       mp.FluxRegion(center=mp.Vector3(x=+r),size=mp.Vector3(0,2*r,2*r),weight=-1),
                       mp.FluxRegion(center=mp.Vector3(y=-r),size=mp.Vector3(2*r,0,2*r),weight=+1),
                       mp.FluxRegion(center=mp.Vector3(y=+r),size=mp.Vector3(2*r,0,2*r),weight=-1),
                       mp.FluxRegion(center=mp.Vector3(z=-r),size=mp.Vector3(2*r,2*r,0),weight=+1),
                       mp.FluxRegion(center=mp.Vector3(z=+r),size=mp.Vector3(2*r,2*r,0),weight=-1))
    
    sim.load_minus_flux_data(box, box_flux_data)

    sim.run(until_after_sources=40)

    scat_flux = mp.get_fluxes(box)[0]
    if mp.am_master():
        print("scat-flux:, {:.6f}".format(scat_flux))

    return scat_flux

if __name__ == '__main__':
    rs = np.logspace(-1,1,30) * wvl/(2*np.pi)
    rcs = np.zeros(len(rs))
    for m in range(len(rs)):
        rcs[m] = mie_scat(rs[m])/(4*np.pi*rs[m]**2)
        print("rcs:, {:.4f}, {:.6f}".format(rs[m],rcs[m]))

    non_zero_idx = np.nonzero(rcs)
    plt.figure(dpi=150)
    plt.loglog(2*np.pi/wvl*rs[non_zero_idx],abs(rcs[non_zero_idx]),'bo-')
    plt.loglog(2*np.pi/wvl*rs[non_zero_idx],np.ones(np.size(non_zero_idx)),'r-')    
    plt.grid(True,which="both",ls="-")
    plt.xlabel('(sphere circumference)/wavelength, 2πr/λ')
    plt.ylabel(r'(scattered power)/(sphere surface area), P/(4πr$^{2}$)')
    plt.title('Mie scattering of a lossless metallic sphere')
    plt.show()
