from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

resolutions = [10., 12., 15., 19., 24., 29., 36., 45., 55., 69., 85., 105., 130., 161., 200]
dr=0.005

ez_harminv_freqs_at_r = [0.1759448601220787, 0.17589545954904895, 0.17585517530855344, 0.1758281671751298,
                         0.17581146477408213, 0.1758026194227426, 0.1757958653854038, 0.17579136921227456,
                         0.17578872608481283, 0.17578677411197724, 0.17578561398088877, 0.17578484129648425,
                         0.17578433056032927, 0.1757839969847718, 0.17578377702683184]
ez_harminv_freqs_at_r_plus_dr=[0.17551207760849705, 0.17546290957473493, 0.17542280821158046, 0.17539592207279095,
                               0.17537929525287474, 0.17537048999006336, 0.17536376661394543, 0.17535929086427382,
                               0.17535665975135972, 0.1753547166544632, 0.17535356179926678, 0.17535367148118317,
                               0.17535572415056194, 0.17535651233115, 0.17535657748276823]
ez_perturb_theory_dw_dr= [-0.08556458559555223, -0.08552711862118463, -0.08549788470502898, -0.08547837070590848,
                          -0.08546626495549077, -0.08545983429512045, -0.08545491268408786, -0.0854516315647058,
                          -0.08544970006548072, -0.08544827179271082, -0.08544742494591485, -0.08544686014784159,
                          -0.08544647465295202, -0.08544623444751495, -0.08544608457833697]
ez_center_diff_dw_dr=[-0.08655650271632842, -0.08650999486280453, -0.08647341939459485, -0.08644902046777148,
                      -0.08643390424147857, -0.08642588653585137, -0.08641975429167226, -0.08641566960014835,
                      -0.08641326669062144, -0.08641149150280802, -0.08641043632439671, -0.08623396306021713,
                      -0.08572128195346584, -0.08549693072436026, -0.08543990881272334]
ez_predicted_freqs = [ez_perturb_theory_dw_dr[i] * dr + ez_harminv_freqs_at_r[i] for i in range(len(ez_harminv_freqs_at_r))]
ez_rel_err_freqs = [abs((ez_predicted_freqs[i] - ez_harminv_freqs_at_r_plus_dr[i]) / ez_harminv_freqs_at_r_plus_dr[i]) for i in range(len(ez_harminv_freqs_at_r))]
ez_rel_err_dw_dr = [abs((ez_perturb_theory_dw_dr[i] - ez_center_diff_dw_dr[i]) / ez_center_diff_dw_dr[i]) for i in range(len(ez_harminv_freqs_at_r))]


hz_harminv_freqs_at_r = [0.20789654295983773, 0.2080258186253275, 0.20813167370381228, 0.2082026698953792,
                         0.2082465614878404, 0.2082697997431897, 0.20828754064681493, 0.20829934933634858,
                         0.20830629066980802, 0.20831141665557598, 0.2083144630991991, 0.20831649209162345,
                         0.20831783323035333, 0.20831870913136713, 0.20831928670696853]
hz_harminv_freqs_at_r_plus_dr=[0.2077080711910781, 0.20778586115865166, 0.20784540385964567, 0.20788181733747052,
                               0.20790172421444852, 0.20791077007195918, 0.20791639390657132, 0.20791901206666977,
                               0.20791976375054255, 0.20791960830661763, 0.20791897935893577, 0.2079190570589992,
                               0.20792051990628624, 0.2079201644855126, 0.20791882018715632]
hz_perturb_theory_dw_dr= [-0.07600868841807405, -0.07561161259060968, -0.07530124517842114, -0.07510833477400632,
                          -0.0750016263230163, -0.07495277031573416, -0.0749222550188294, -0.07490807060134752,
                          -0.07490405512083827, -0.07490503543658283, -0.0749086049770053, -0.07491339596721278,
                          -0.07491851535123699, -0.0749234548641879, -0.07492802777136065]
hz_center_diff_dw_dr=[-0.03769435375192698, -0.047991493335169944, -0.05725396883332068, -0.06417051158173481,
                      -0.06896745467837584, -0.07180593424610526, -0.07422934804872106, -0.07606745393576309,
                      -0.07730538385309349, -0.07836166979167114, -0.07909674805266498, -0.07948700652484764,
                      -0.0794626648134178, -0.07970892917090744, -0.08009330396244185]
hz_predicted_freqs = [hz_perturb_theory_dw_dr[i] * dr + hz_harminv_freqs_at_r[i] for i in range(len(hz_harminv_freqs_at_r))]
hz_rel_err_freqs = [abs((hz_predicted_freqs[i] - hz_harminv_freqs_at_r_plus_dr[i]) / hz_harminv_freqs_at_r_plus_dr[i]) for i in range(len(hz_harminv_freqs_at_r))]
hz_rel_err_dw_dr = [abs((hz_perturb_theory_dw_dr[i] - hz_center_diff_dw_dr[i]) / hz_center_diff_dw_dr[i]) for i in range(len(hz_harminv_freqs_at_r))]


plt.figure(dpi=150)
plt.loglog(resolutions, ez_rel_err_dw_dr, 'bo-', label='Parallel Fields')
plt.loglog(resolutions, hz_rel_err_dw_dr, 'ro-', label='Perpendicular Fields')
plt.loglog(resolutions,np.divide(10**1.5,resolutions),'k--',label="linear reference")
plt.loglog(resolutions,10**0.5/np.power(resolutions,2),'k-',label="quadratic reference")
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