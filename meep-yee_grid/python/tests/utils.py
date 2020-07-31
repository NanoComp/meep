import numpy as np


def compare_arrays(test_instance, exp, res, tol=1e-3):
    exp_1d = exp.ravel()
    res_1d = res.ravel()

    norm_exp = np.linalg.norm(exp_1d)
    norm_res = np.linalg.norm(res_1d)

    if norm_exp == 0:
        test_instance.assertEqual(norm_res, 0)
    else:
        diff = np.linalg.norm(res_1d - exp_1d) / norm_exp
        test_instance.assertLess(diff, tol)
