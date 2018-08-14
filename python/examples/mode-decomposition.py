# -*- coding: utf-8 -*-

import meep as mp # pixels/μm

resolution = 60

w1 = 1            # width of waveguide 1
w2 = 2            # width of waveguide 2
Lw = 10           # length of waveguide 1/2

Si = mp.Medium(epsilon=12.0)

dair = 3.0
dpml = 5.0

sx = dpml+Lw+dpml
sy = dpml+dair+w2+dair+dpml
cell_size = mp.Vector3(sx,sy,0)

prism_x = sx+1
half_w1 = 0.5*w1
half_w2 = 0.5*w2

# straight waveguide
vertices = [mp.Vector3(-prism_x,half_w1),
            mp.Vector3(prism_x,half_w1),
            mp.Vector3(prism_x,-half_w1),
            mp.Vector3(-prism_x,-half_w1)]
    
boundary_layers = [mp.PML(dpml)]

# mode wavelength
lcen = 6.67
# mode frequency
fcen = 1/lcen

sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen,fwidth=0.2*fcen),
                              component=mp.Ez,
                              center=mp.Vector3(-0.5*sx+dpml+0.2*Lw,0,0),
                              size=mp.Vector3(0,sy-2*dpml,0),
                              eig_match_freq=True,
                              eig_parity=mp.ODD_Z+mp.EVEN_Y)]

symmetries = [mp.Mirror(mp.Y)]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=boundary_layers,
                    geometry=[mp.Prism(vertices,height=mp.inf,material=Si)],
                    sources=sources,
                    symmetries=symmetries)

mon_pt = mp.Vector3(-0.5*sx+dpml+0.5*Lw,0,0)
flux = sim.add_flux(fcen,0,1,mp.FluxRegion(center=mon_pt,size=mp.Vector3(0,sy-2*dpml,0)))

sim.run(until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,mon_pt,1e-10))

incident_coeffs, vgrp, kpoints = sim.get_eigenmode_coefficients(flux,[1],eig_parity=mp.ODD_Z+mp.EVEN_Y)
incident_flux = mp.get_fluxes(flux)
incident_flux_data = sim.get_flux_data(flux)

for m in range(5):
    Lt = 2**m
    sx = dpml+Lw+Lt+Lw+dpml
    cell_size = mp.Vector3(sx,sy,0)

    prism_x = sx+1
    half_Lt = 0.5*Lt
    
    # linear taper
    vertices = [mp.Vector3(-prism_x,half_w1),
                mp.Vector3(-half_Lt,half_w1),
                mp.Vector3(half_Lt,half_w2),
                mp.Vector3(prism_x,half_w2),
                mp.Vector3(prism_x,-half_w2),
                mp.Vector3(half_Lt,-half_w2),
                mp.Vector3(-half_Lt,-half_w1),
                mp.Vector3(-prism_x,-half_w1)]

    sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen,fwidth=0.2*fcen),
                                  component=mp.Ez,
                                  center=mp.Vector3(-0.5*sx+dpml+0.2*Lw,0,0),
                                  size=mp.Vector3(0,sy-2*dpml,0),
                                  eig_match_freq=True,
                                  eig_parity=mp.ODD_Z+mp.EVEN_Y)]

    sim.reset_meep()
    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=boundary_layers,
                        geometry=[mp.Prism(vertices,height=mp.inf,material=Si)],
                        sources=sources,
                        symmetries=symmetries)

    mon_pt = mp.Vector3(-0.5*sx+dpml+0.5*Lw,0,0)
    refl_flux = sim.add_flux(fcen,0,1,mp.FluxRegion(center=mon_pt,size=mp.Vector3(0,sy-2*dpml,0)))

    sim.load_minus_flux_data(refl_flux,incident_flux_data)
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,mon_pt,1e-10))
        
    coeffs, vgrp, kpoints = sim.get_eigenmode_coefficients(refl_flux,[1],eig_parity=mp.ODD_Z+mp.EVEN_Y)
    taper_flux = mp.get_fluxes(refl_flux)
    print("refl:, {}, {:.8f}, {:.8f}".format(Lt,abs(coeffs[0,0,1])**2/abs(incident_coeffs[0,0,0])**2,-taper_flux[0]/incident_flux[0]))
