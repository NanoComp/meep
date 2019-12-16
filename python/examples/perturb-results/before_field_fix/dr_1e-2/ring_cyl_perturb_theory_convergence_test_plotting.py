from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

resolutions = [100., 119., 141., 168., 200.]
dr=0.01

ez_harminv_freqs_at_r = [0.17578499188503546, 0.17578451594130998, 0.17578418683861435, 0.17578394599265096, 0.17578377702683184]
ez_harminv_freqs_at_r_plus_dr=[0.17493286732441568, 0.1749319032176484, 0.1749304377581377, 0.1749311300203282, 0.17493166203553098]
ez_perturb_theory_dw_dr= [-0.08544696870452952, -0.0854466133618931, -0.08544638217797605, -0.08544620705772135, -0.08544608457833697]
ez_center_diff_dw_dr=[-0.08521245606197825, -0.08526127236615777, -0.08537490804766401, -0.08528159723227569, -0.08521149913008619]
ez_predicted_freqs = [ez_perturb_theory_dw_dr[i] * dr + ez_harminv_freqs_at_r[i] for i in range(len(ez_harminv_freqs_at_r))]
ez_rel_err_freqs = [abs((ez_predicted_freqs[i] - ez_harminv_freqs_at_r_plus_dr[i]) / ez_harminv_freqs_at_r_plus_dr[i]) for i in range(len(ez_harminv_freqs_at_r))]
ez_rel_err_dw_dr = [abs((ez_perturb_theory_dw_dr[i] - ez_center_diff_dw_dr[i]) / ez_center_diff_dw_dr[i]) for i in range(len(ez_harminv_freqs_at_r))]


hz_harminv_freqs_at_r = [0.20831609666994255, 0.20831734643587063, 0.20831821060456387, 0.20831884302817685, 0.20831928670696853]
hz_harminv_freqs_at_r_plus_dr=[0.20751800980999938, 0.2075217880913565, 0.20752200005586988, 0.2075226964202726, 0.20752119916867948]
hz_perturb_theory_dw_dr= [-0.07491224947164242, -0.07491639745122385, -0.07492042959596894, -0.07492439294728007, -0.07492802777136065]
hz_center_diff_dw_dr=[-0.07980868599431756, -0.07955583445141412, -0.07962105486939985, -0.07961466079042523, -0.07980875382890551]
hz_predicted_freqs = [hz_perturb_theory_dw_dr[i] * dr + hz_harminv_freqs_at_r[i] for i in range(len(hz_harminv_freqs_at_r))]
hz_rel_err_freqs = [abs((hz_predicted_freqs[i] - hz_harminv_freqs_at_r_plus_dr[i]) / hz_harminv_freqs_at_r_plus_dr[i]) for i in range(len(hz_harminv_freqs_at_r))]
hz_rel_err_dw_dr = [abs((hz_perturb_theory_dw_dr[i] - hz_center_diff_dw_dr[i]) / hz_center_diff_dw_dr[i]) for i in range(len(hz_harminv_freqs_at_r))]


plt.figure(dpi=150)
plt.loglog(resolutions, ez_rel_err_dw_dr, 'bo-', label='Parallel Fields')
plt.loglog(resolutions, hz_rel_err_dw_dr, 'ro-', label='Perpendicular Fields')
plt.loglog(resolutions,np.divide(10**1.5,resolutions),'k--',label="linear reference")
plt.loglog(resolutions,10**1/np.power(resolutions,2),'k-',label="quadratic reference")
plt.grid(True, which='both', ls='-')
plt.xlabel('resolution')
plt.ylabel('relative error between $dω/dR$')
plt.legend(loc='lower right')
plt.title('Comparison of Perturbation Theory and \nCenter-Difference Calculations in Finding $dω/dR$\n$dR={}$'.format(dr))
plt.tight_layout()
# plt.show()
plt.savefig('ring_cyl_perturbation_theory.dw_dR_error.png')

# plt.figure(dpi=150)
# plt.loglog(resolutions, ez_rel_err_freqs, 'bo-', label='Parallel Fields')
# plt.loglog(resolutions, hz_rel_err_freqs, 'ro-', label='Perpendicular Fields')
# plt.grid(True, which='both', ls='-')
# plt.xlabel('resolution')
# plt.ylabel('relative error between $ω(R+dR)$')
# plt.legend(loc='upper left')
# plt.title('Comparison of resonance frequencies at $R+dR$ predicted by\nperturbation theory and found with Harminv')
# plt.tight_layout()
# # plt.show()
# plt.savefig('ring_cyl_perturbation_theory.freqs_error.png')
# # plt.clf()