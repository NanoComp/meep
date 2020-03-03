import meep as mp
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import meep_adjoint as mpa

# --------------------------- #
# Design filter
# --------------------------- #

freq_min = 1/1.7
freq_max = 1/1.4
nfreq = 200
freqs = np.linspace(freq_min,freq_max,nfreq)
resolution = 10
run_time = 4000

# --------------------------- #
# Run normal simulation
# --------------------------- #

cell = [16,8,0]
geometry = [mp.Block(mp.Vector3(mp.inf,1,mp.inf),
                     center=mp.Vector3(),
                     material=mp.Medium(epsilon=12))]
fcen = 1/1.55
gauss_src = mp.GaussianSource(frequency=fcen,fwidth=0.1/1.55,is_integrated=False)
sources = [mp.EigenModeSource(gauss_src,
                     eig_band=1,
                     size=[0,8],
                     center=[-6,0])]
pml_layers = [mp.PML(1.0)]

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    resolution=resolution)
mon = sim.add_dft_fields([mp.Ez],freq_min,freq_max,nfreq,center=[0,0],size=[0,1])

sim.run(until=run_time)

fields = np.zeros((nfreq,),dtype=np.complex128)
# just store one spatial point for each freq
for f in range(nfreq):
    fields[f] = sim.get_dft_array(mon,mp.Ez,f)[10]

# build simple bandpass filter
dt = sim.fields.dt
fs = 1/dt
freqs_scipy = freqs * np.pi / (fs/2)
num_taps = 321
taps = signal.firwin(num_taps,[freq_min/(fs/2), freq_max/(fs/2)], pass_zero='bandstop',window='boxcar')

w,h = signal.freqz(taps,worN=freqs_scipy)

'''plt.figure()
plt.plot(w,np.abs(h))
plt.show()
quit()'''

# frequency domain calculation
desired_fields = h * fields

# --------------------------- #
# Run filtered simulation
# --------------------------- #

filtered_src = mp.FilteredSource(fcen,freqs,h,gauss_src)

sources = [mp.EigenModeSource(filtered_src,
                     eig_band=1,
                     size=[0,8],
                     center=[-6,0])]

sim.reset_meep()
sim.change_sources(sources)

mon = sim.add_dft_fields([mp.Ez],freq_min,freq_max,nfreq,center=[0,0],size=[0,1])
sim.run(until=run_time)

fields_filtered = np.zeros((nfreq,),dtype=np.complex128)
# just store one spatial point for each freq
for f in range(nfreq):
    fields_filtered[f] = sim.get_dft_array(mon,mp.Ez,f)[10]

# --------------------------- #
# Compare results
# --------------------------- #
print(np.abs(fields_filtered/desired_fields), np.angle(fields_filtered/desired_fields))

plt.figure()
plt.subplot(2,1,1)
plt.semilogy(freqs,np.abs(desired_fields),label='Frequency Domain')
plt.semilogy(freqs,np.abs(fields_filtered),'-.',label='Time Domain')
#plt.plot(freqs,np.abs(fields),'--',label='Gaussian Src')
#plt.plot(freqs,np.abs(h),'--',label='Window')
plt.grid()
plt.xlabel('Frequency')
plt.ylabel('Magnitude')
plt.legend()

plt.subplot(2,1,2)
plt.plot(freqs,np.unwrap(np.angle(desired_fields)),label='Frequency Domain')
plt.plot(freqs,np.unwrap(np.angle(fields_filtered)),'-.',label='Time Domain')
#plt.plot(freqs,np.unwrap(np.angle(fields)),'--',label='Gaussian Src')
#plt.plot(freqs,np.unwrap(np.angle(h)),'--',label='Window')
plt.grid()
plt.xlabel('Frequency')
plt.ylabel('Angle')
plt.legend()

plt.tight_layout()
plt.savefig('filtered.png')

plt.show()