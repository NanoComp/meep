from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

resolutions = [100., 108., 117., 126., 136., 147., 159., 171., 185., 200.]
dr=0.0001

ez_harminv_freqs_at_r = [0.17578499188503546, 0.1757847608103432, 0.17578455537857748, 0.17578439237721472,
                         0.17578424785029104, 0.17578412168581883, 0.1757840128067264, 0.17578392603253604,
                         0.17578384536183025, 0.17578377702683184]
ez_harminv_freqs_at_r_plus_dr=[0.17577642308106872, 0.17577619201658898, 0.1757759866056893, 0.1757758236159623,
                               0.1757756791068255, 0.1757755529501322, 0.17577544409278273, 0.17577535732226707,
                               0.17577527665712633, 0.1757752083234335]
ez_perturb_theory_dw_dr= [-0.08544696870452953, -0.08544680156782197, -0.08544665318856819, -0.08544653049272684,
                          -0.08544642608708171, -0.08544633424219143, -0.08544625602889686, -0.08544618770170438,
                          -0.08544613110618413, -0.08544608457833698]
ez_center_diff_dw_dr=[-0.08568803966740868, -0.08568793754215598, -0.08568772888184473, -0.08568761252408796,
                      -0.08568743465525719, -0.08568735686637075, -0.08568713943668538, -0.08568710268969104,
                      -0.08568704703920682, -0.08568703398353916]
ez_predicted_freqs = [ez_perturb_theory_dw_dr[i] * dr + ez_harminv_freqs_at_r[i] for i in range(len(ez_harminv_freqs_at_r))]
ez_rel_err_freqs = [abs((ez_predicted_freqs[i] - ez_harminv_freqs_at_r_plus_dr[i]) / ez_harminv_freqs_at_r_plus_dr[i]) for i in range(len(ez_harminv_freqs_at_r))]
ez_rel_err_dw_dr = [abs((ez_perturb_theory_dw_dr[i] - ez_center_diff_dw_dr[i]) / ez_center_diff_dw_dr[i]) for i in range(len(ez_harminv_freqs_at_r))]


hz_harminv_freqs_at_r = [0.20831609666994255, 0.20831670344450806, 0.20831724286749728, 0.20831767089502762,
                         0.20831805040119272, 0.2083183816888786, 0.20831866759388476, 0.2083188954440263,
                         0.20831910727901737, 0.20831928670696853]
hz_harminv_freqs_at_r_plus_dr=[0.20830833969511953, 0.20830892347644797, 0.2083094410956956, 0.20830985065315252,
                               0.20831021270984773, 0.20831052771210723, 0.2083107985901013, 0.20831101362746496,
                               0.20831121272481284, 0.20831138058057572]
hz_perturb_theory_dw_dr= [-0.08050386030223725, -0.08045908261479044, -0.08041608995367247, -0.08037929760166317,
                          -0.08034414799917983, -0.08031105245314289, -0.08028020744599708, -0.08025370103709917,
                          -0.08022713578743966, -0.08020283538083425]
hz_center_diff_dw_dr=[-0.07756974823025509, -0.07779968060089848, -0.07801771801679847, -0.07820241875095002,
                      -0.07837691344991793, -0.07853976771360349, -0.07869003783467221, -0.07881816561333688,
                      -0.07894554204523896, -0.0790612639281485]
hz_predicted_freqs = [hz_perturb_theory_dw_dr[i] * dr + hz_harminv_freqs_at_r[i] for i in range(len(hz_harminv_freqs_at_r))]
hz_rel_err_freqs = [abs((hz_predicted_freqs[i] - hz_harminv_freqs_at_r_plus_dr[i]) / hz_harminv_freqs_at_r_plus_dr[i]) for i in range(len(hz_harminv_freqs_at_r))]
hz_rel_err_dw_dr = [abs((hz_perturb_theory_dw_dr[i] - hz_center_diff_dw_dr[i]) / hz_center_diff_dw_dr[i]) for i in range(len(hz_harminv_freqs_at_r))]


plt.figure(dpi=150)
plt.loglog(resolutions, ez_rel_err_dw_dr, 'bo-', label='Parallel Fields')
plt.loglog(resolutions, hz_rel_err_dw_dr, 'ro-', label='Perpendicular Fields')
plt.loglog(resolutions,np.divide(10**0.6,resolutions),'k--',label="linear reference")
plt.loglog(resolutions,10**1.25/np.power(resolutions,2),'k-',label="quadratic reference")
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