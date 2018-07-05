# -*- coding: utf-8 -*-

import meep as mp
import math

resolution = 40        # pixels/μm

dsub = 3.0             # substrate thickness
dpad = 3.0             # padding between grating and pml
gp = 10.0              # grating period
gh = 0.5               # grating height
gdc = 0.5              # grating duty cycle

dpml = 1.0             # PML thickness
sx = dpml+dsub+gh+dpad+dpml
sy = gp

cell_size = mp.Vector3(sx,sy,0)
pml_layers = [ mp.PML(thickness=dpml,direction=mp.X) ]

wvl_min = 0.4           # min wavelength
wvl_max = 0.6           # max wavelength
fmin = 1/wvl_max        # min frequency
fmax = 1/wvl_min        # max frequency
fcen = 0.5*(fmin+fmax)  # center frequency
df = fmax-fmin          # frequency width

src_pos = -0.5*sx+dpml+0.5*dsub
sources = [ mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ez, center=mp.Vector3(src_pos,0,0), size=mp.Vector3(0,sy,0)) ]

k_point = mp.Vector3(0,0,0)

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    k_point=k_point,
                    sources=sources)

nfreq = 21
xm = 0.5*sx-dpml-0.5*dpad
eig_mon = sim.add_eigenmode(fcen, df, nfreq, mp.FluxRegion(center=mp.Vector3(xm,0,0), size=mp.Vector3(0,sy,0)))

sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(xm,0,0), 1e-9))

coeffs0, vgrps0, kpoints0 = sim.get_eigenmode_coefficients(eig_mon, [1], eig_parity=mp.ODD_Z+mp.EVEN_Y)

sim.reset_meep()

glass = mp.Medium(index=1.5)

geometry = [ mp.Block(material=glass, size=mp.Vector3(dpml+dsub,mp.inf,mp.inf), center=mp.Vector3(-0.5*sx+0.5*(dpml+dsub),0,0)),
             mp.Block(material=glass, size=mp.Vector3(gh,gdc*gp,mp.inf), center=mp.Vector3(-0.5*sx+dpml+dsub+0.5*gh,0,0)) ]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    k_point=k_point,
                    sources=sources)

eig_mon = sim.add_eigenmode(fcen, df, nfreq, mp.FluxRegion(center=mp.Vector3(xm,0,0), size=mp.Vector3(0,sy,0)))

sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(xm,0,0), 1e-9))

freqs = mp.get_eigenmode_freqs(eig_mon)

kx = lambda m,freq: math.sqrt(freq**2 - (m/gp)**2)
theta_out = lambda m,freq: math.acos(kx(m,freq)/freq)

nmode = 10
for nm in range(nmode):
  coeffs, vgrps, kpoints = sim.get_eigenmode_coefficients(eig_mon, [nm+1], eig_parity=mp.ODD_Z+mp.EVEN_Y)
  for nf in range(nfreq):
    mode_wvl = 1/freqs[nf]
    mode_angle = math.degrees(theta_out(nm,freqs[nf]))
    mode_tran = abs(coeffs[0,nf,0])**2/abs(coeffs0[0,nf,0])**2
    print("grating{}:, {:.5f}, {:.2f}, {:.8f}".format(nm,mode_wvl,mode_angle,mode_tran))
