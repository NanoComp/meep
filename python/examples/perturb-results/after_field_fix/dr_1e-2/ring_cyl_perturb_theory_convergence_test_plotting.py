from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

resolutions = [100., 108., 117., 126., 136., 147., 159., 171., 185., 200.]
dr=0.01

ez_harminv_freqs_at_r = [0.17578499188503546, 0.1757847608103432, 0.17578455537857748, 0.17578439237721472,
                         0.17578424785029104, 0.17578412168581883, 0.1757840128067264, 0.17578392603253604,
                         0.17578384536183025, 0.17578377702683184]
ez_harminv_freqs_at_r_plus_dr=[0.17493286732441568, 0.17493253208433043, 0.17493202684132708, 0.17493145106882438,
                               0.1749307775769501, 0.17493003234628496, 0.17493061361979617, 0.17493125525340514,
                               0.1749316030925222, 0.17493166203553098]
ez_perturb_theory_dw_dr= [-0.08544696870452953, -0.08544680156782197, -0.08544665318856819, -0.08544653049272684,
                          -0.08544642608708171, -0.08544633424219143, -0.08544625602889686, -0.08544618770170438,
                          -0.08544613110618413, -0.08544608457833698]
ez_center_diff_dw_dr=[-0.08521245606197825, -0.08522287260127603, -0.08525285372504021, -0.0852941308390337,
                      -0.08534702733409283, -0.08540893395338756, -0.0853399186930226, -0.08526707791308985,
                      -0.08522422693080511, -0.08521149913008619]
ez_predicted_freqs = [ez_perturb_theory_dw_dr[i] * dr + ez_harminv_freqs_at_r[i] for i in range(len(ez_harminv_freqs_at_r))]
ez_rel_err_freqs = [abs((ez_predicted_freqs[i] - ez_harminv_freqs_at_r_plus_dr[i]) / ez_harminv_freqs_at_r_plus_dr[i]) for i in range(len(ez_harminv_freqs_at_r))]
ez_rel_err_dw_dr = [abs((ez_perturb_theory_dw_dr[i] - ez_center_diff_dw_dr[i]) / ez_center_diff_dw_dr[i]) for i in range(len(ez_harminv_freqs_at_r))]


hz_harminv_freqs_at_r = [0.20831609666994255, 0.20831670344450806, 0.20831724286749728, 0.20831767089502762,
                         0.20831805040119272, 0.2083183816888786, 0.20831866759388476, 0.2083188954440263,
                         0.20831910727901737, 0.20831928670696853]
hz_harminv_freqs_at_r_plus_dr=[0.20751800980784021, 0.2075202825560755, 0.20752161500695757, 0.20752213954018445,
                               0.20752214528337207, 0.2075217367722817, 0.20752240164156374, 0.20752270181193921,
                               0.20752226495840775, 0.20752119916663955]
hz_perturb_theory_dw_dr= [-0.08050386030223725, -0.08045908261479044, -0.08041608995367247, -0.08037929760166317,
                          -0.08034414799917983, -0.08031105245314289, -0.08028020744599708, -0.08025370103709917,
                          -0.08022713578743966, -0.08020283538083425]
hz_center_diff_dw_dr=[-0.07980868621023374, -0.07964208884325696, -0.0795627860539716, -0.07955313548431708,
                      -0.07959051178206555, -0.07966449165968947, -0.0796265952321018, -0.07961936320870777,
                      -0.07968423206096142, -0.07980875403289789]
hz_predicted_freqs = [hz_perturb_theory_dw_dr[i] * dr + hz_harminv_freqs_at_r[i] for i in range(len(hz_harminv_freqs_at_r))]
hz_rel_err_freqs = [abs((hz_predicted_freqs[i] - hz_harminv_freqs_at_r_plus_dr[i]) / hz_harminv_freqs_at_r_plus_dr[i]) for i in range(len(hz_harminv_freqs_at_r))]
hz_rel_err_dw_dr = [abs((hz_perturb_theory_dw_dr[i] - hz_center_diff_dw_dr[i]) / hz_center_diff_dw_dr[i]) for i in range(len(hz_harminv_freqs_at_r))]


plt.figure(dpi=150)
plt.loglog(resolutions, ez_rel_err_dw_dr, 'bo-', label='Parallel Fields')
plt.loglog(resolutions, hz_rel_err_dw_dr, 'ro-', label='Perpendicular Fields')
plt.loglog(resolutions,np.divide(10**0.25,resolutions),'k--',label="linear reference")
plt.loglog(resolutions,10**1/np.power(resolutions,2),'k-',label="quadratic reference")
plt.grid(True, which='both', ls='-')
plt.xlabel('resolution')
plt.ylabel('relative error between $dω/dR$')
plt.legend(loc='lower left')
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