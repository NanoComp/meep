from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

resolutions = [100., 108., 117., 126., 136., 147., 159., 171., 185., 200.]
dr=0.001

ez_harminv_freqs_at_r = [0.17578499188503546, 0.1757847608103432, 0.17578455537857748, 0.17578439237721472,
                         0.17578424785029104, 0.17578412168581883, 0.1757840128067264, 0.17578392603253604,
                         0.17578384536183025, 0.17578377702683184]
ez_harminv_freqs_at_r_plus_dr=[0.17569917255321835, 0.17569894167519365, 0.17569873642999956, 0.17569857357192784,
                               0.17569842917935144, 0.17569830312460166, 0.17569819435506692, 0.17569810765454064,
                               0.17569802705459053, 0.17569795877622293 ]
ez_perturb_theory_dw_dr= [-0.08544696870452953, -0.08544680156782197, -0.08544665318856819, -0.08544653049272684,
                          -0.08544642608708171, -0.08544633424219143, -0.08544625602889686, -0.08544618770170438,
                          -0.08544613110618413, -0.08544608457833698]
ez_center_diff_dw_dr=[-0.08581933181711632, -0.0858191351495452, -0.08581894857792594, -0.08581880528688024,
                      -0.08581867093959694, -0.08581856121717135, -0.0858184516594751, -0.08581837799540026,
                      -0.08581830723972117, -0.08581825060891002]
ez_predicted_freqs = [ez_perturb_theory_dw_dr[i] * dr + ez_harminv_freqs_at_r[i] for i in range(len(ez_harminv_freqs_at_r))]
ez_rel_err_freqs = [abs((ez_predicted_freqs[i] - ez_harminv_freqs_at_r_plus_dr[i]) / ez_harminv_freqs_at_r_plus_dr[i]) for i in range(len(ez_harminv_freqs_at_r))]
ez_rel_err_dw_dr = [abs((ez_perturb_theory_dw_dr[i] - ez_center_diff_dw_dr[i]) / ez_center_diff_dw_dr[i]) for i in range(len(ez_harminv_freqs_at_r))]


hz_harminv_freqs_at_r = [0.20831609666994255, 0.20831670344450806, 0.20831724286749728, 0.20831767089502762,
                         0.20831805040119272, 0.2083183816888786, 0.20831866759388476, 0.2083188954440263,
                         0.20831910727901737, 0.20831928670696853]
hz_harminv_freqs_at_r_plus_dr=[0.2082381119328325, 0.2082384973208359, 0.208238826759217, 0.20823907703700542,
                               0.20823928858511062, 0.20823946324299591, 0.20823960453547732, 0.2082397092335266,
                               0.2082397986774414, 0.20823986694842345]
hz_perturb_theory_dw_dr= [-0.08050386030223725, -0.08045908261479044, -0.08041608995367247, -0.08037929760166317,
                          -0.08034414799917983, -0.08031105245314289, -0.08028020744599708, -0.08025370103709917,
                          -0.08022713578743966, -0.08020283538083425]
hz_center_diff_dw_dr=[-0.07798473711004283, -0.07820612367215318, -0.07841610828027146, -0.07859385802219676,
                      -0.07876181608210131, -0.07891844588267527, -0.07906305840743588, -0.07918621049968211,
                      -0.07930860157595587, -0.0794197585450851]
hz_predicted_freqs = [hz_perturb_theory_dw_dr[i] * dr + hz_harminv_freqs_at_r[i] for i in range(len(hz_harminv_freqs_at_r))]
hz_rel_err_freqs = [abs((hz_predicted_freqs[i] - hz_harminv_freqs_at_r_plus_dr[i]) / hz_harminv_freqs_at_r_plus_dr[i]) for i in range(len(hz_harminv_freqs_at_r))]
hz_rel_err_dw_dr = [abs((hz_perturb_theory_dw_dr[i] - hz_center_diff_dw_dr[i]) / hz_center_diff_dw_dr[i]) for i in range(len(hz_harminv_freqs_at_r))]


plt.figure(dpi=150)
plt.loglog(resolutions, ez_rel_err_dw_dr, 'bo-', label='Parallel Fields')
plt.loglog(resolutions, hz_rel_err_dw_dr, 'ro-', label='Perpendicular Fields')
plt.loglog(resolutions,np.divide(10**0.6,resolutions),'k--',label="linear reference")
plt.loglog(resolutions,10**1.5/np.power(resolutions,2),'k-',label="quadratic reference")
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