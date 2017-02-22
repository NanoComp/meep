// This is a version of 2D_convergence.cpp using a
// material_function_for_polygons instead of eps function.

#include <meep.hpp>
using namespace meep;
#include <stdlib.h>
using namespace std;


void rotate_polygon(double* pts, const int numpts, const double angle)
{
    for (int i = 0; i < numpts; ++i) {
        double x = pts[2 * i];
        double y = pts[2 * i + 1];
        pts[2 * i] = x * cos(angle) - y * sin(angle);
        pts[2 * i + 1] = x * sin(angle) + y * cos(angle);
    }
}

void shift_polygon(double* pts, const int numpts, const vec &shift)
{
    for (int i = 0; i < numpts; ++i) {
        pts[2 * i] += shift.x();
        pts[2 * i + 1] += shift.y();
    }
}

// uniform pseudo-random number in [min,max]
static double urand(double min, double max)
{
    return (rand() * ((max - min) / RAND_MAX) + min);
}

double get_fill_ratio(const double angle, const double xshift, const double yshift)
{
    if (angle < -pi || angle > pi)
        abort("testing angle must be between -pi and +pi");
    double ta(tan(angle));
    double x_top, x_bot, y_lft, y_rgt;
    x_top = (yshift - 0.5) * ta + xshift;
    x_bot = (yshift + 0.5) * ta + xshift;
    y_lft = (xshift + 0.5) / ta + yshift;
    y_rgt = (xshift - 0.5) / ta + yshift;

    double pts[6][2];
    int numpts = 2;

    if (angle < -pi/2.0) {
        pts[0][0] = 0.5;
        pts[0][1] = 0.5;
        if (x_top >= -0.5) {
            pts[1][0] = x_top;
            pts[1][1] = 0.5;
        }
        else {
            numpts++;
            pts[1][0] = -0.5;
            pts[1][1] = 0.5;
            pts[2][0] = -0.5;
            pts[2][1] = y_lft;
        }
        if (x_bot <= 0.5) {
            pts[numpts][0] = x_bot;
            pts[numpts][1] = -0.5;
            numpts++;
            pts[numpts][0] = 0.5;
            pts[numpts][1] = -0.5;
        }
        else {
            pts[numpts][0] = 0.5;
            pts[numpts][1] = y_rgt;
        }
        pts[numpts+1][0] = 0.5;
        pts[numpts+1][1] = 0.5;
        numpts += 2;
    }
    else if (angle < 0.0) {
        pts[0][0] = -0.5;
        pts[0][1] = 0.5;
        if (y_lft >= -0.5) {
            pts[1][0] = -0.5;
            pts[1][1] = y_lft;
        }
        else {
            numpts++;
            pts[1][0] = -0.5;
            pts[1][1] = -0.5;
            pts[2][0] = x_bot;
            pts[2][1] = -0.5;
        }
        if (y_rgt <= 0.5) {
            pts[numpts][0] = 0.5;
            pts[numpts][1] = y_rgt;
            numpts++;
            pts[numpts][0] = 0.5;
            pts[numpts][1] = 0.5;
        }
        else {
            pts[numpts][0] = x_top;
            pts[numpts][1] = 0.5;
        }
        pts[numpts+1][0] = -0.5;
        pts[numpts+1][1] = 0.5;
        numpts += 2;
    }
    else if (angle < pi/2.0){
        pts[0][0] = -0.5;
        pts[0][1] = -0.5;
        if (x_bot <= 0.5) {
            pts[1][0] = x_bot;
            pts[1][1] = -0.5;
        }
        else {
            numpts++;
            pts[1][0] = 0.5;
            pts[1][1] = -0.5;
            pts[2][0] = 0.5;
            pts[2][1] = y_rgt;
        }
        if (x_top >= -0.5) {
            pts[numpts][0] = x_top;
            pts[numpts][1] = 0.5;
            numpts++;
            pts[numpts][0] = -0.5;
            pts[numpts][1] = 0.5;
        }
        else {
            pts[numpts][0] = -0.5;
            pts[numpts][1] = y_lft;
        }
        pts[numpts+1][0] = -0.5;
        pts[numpts+1][1] = -0.5;
        numpts += 2;
    }
    else {
        pts[0][0] = 0.5;
        pts[0][1] = -0.5;
        if (y_rgt <= 0.5) {
            pts[1][0] = 0.5;
            pts[1][1] = y_rgt;
        }
        else {
            numpts++;
            pts[1][0] = 0.5;
            pts[1][1] = 0.5;
            pts[2][0] = x_top;
            pts[2][1] = 0.5;
        }
        if (y_lft >= -0.5) {
            pts[numpts][0] = -0.5;
            pts[numpts][1] = y_lft;
            numpts++;
            pts[numpts][0] = -0.5;
            pts[numpts][1] = -0.5;
        }
        else {
            pts[numpts][0] = x_bot;
            pts[numpts][1] = -0.5;
        }
        pts[numpts+1][0] = 0.5;
        pts[numpts+1][1] = -0.5;
        numpts += 2;
    }
    polygon pol = polygon(pts, numpts);
    return pol.get_area();
};

/** Tests the normal direction, mean eps and mean 1/eps of a polygonal
 *  interface in 2D at multiple randomly chosen angles and offsets. */
void check_normals_2D(const int num_random_trials = 10000)
{
    master_printf("Checking interface of polygonal material with air (2D).\n");
    const int w = 6;
    grid_volume gv(voltwo(w, w, 1));
    volume vol(gv.surroundings());
    vec center = gv.center();
    ivec icenter = gv.icenter();

    for (int i = 0; i < num_random_trials; ++i) {
        double pts[8] = {
            -w, -w,
            0, -w,
            0,  w,
            -w,  w
        };
        double test_angle = urand(-pi, pi);
        vec test_shift = vec(urand(-0.5, 0.5), urand(-0.5, 0.5));

        rotate_polygon(pts, 4, test_angle);
        shift_polygon(pts, 4, center + test_shift);
        polygon pol = polygon(pts, 4, 2);
        material_function_for_polygons matfun(gv);
        matfun.add_polygon(pol, 12.0);

        vec grad(matfun.normal_vector(icenter));
        double grad_ang = atan2(grad.y(), grad.x());
        double fratio(get_fill_ratio(test_angle, test_shift.x(), test_shift.y()));
        double test_meps = fratio * 11 + 1;
        double meps = matfun.mean_eps(icenter);
        double test_minveps = 1 - fratio * 11./12.;
        double minveps = matfun.mean_inveps(icenter);


        if (fabs(test_angle - grad_ang) > 1e-8 ||
            fabs(test_minveps - minveps) > 1e-8 ||
            fabs(test_meps - meps) > 1e-8
        )
        {
            master_printf("error at #%i: test shift: (%.14f, %.14f) ; test angle: %.14f\n",
                          i, test_shift.x(), test_shift.y(), test_angle * 180 / pi);
            master_printf("expected angle: %.14f ; mean_eps: %.14f ; mean_inv_eps: %.14f\n",
                          test_angle * 180 / pi, test_meps, test_minveps);
            master_printf("     got angle: %.14f ; mean_eps: %.14f ; mean_inv_eps: %.14f\n\n",
                          grad_ang * 180 / pi, meps, minveps);
            abort("error");
        }
    }
}

/** Tests the normal direction, mean eps and mean 1/eps of a polygonal
 *  interface in 2D between 2 different materials at multiple randomly
 *  chosen angles and offsets. Also tests usage of inner polygons. */
void check_normals_2D_2materials(const int num_random_trials = 10000)
{
    master_printf("Checking interface between two polygonal materials (2D).\n");
    const int w = 6;
    grid_volume gv(voltwo(w, w, 1));
    volume vol(gv.surroundings());
    vec center = gv.center();
    ivec icenter = gv.icenter();

    for (int i = 0; i < num_random_trials; ++i) {
        double test_angle = urand(-pi, pi);
        vec test_shift = vec(urand(-0.5, 0.5), urand(-0.5, 0.5));

        double pts1[8] = {
            -w, -w,
            w, -w,
            w,  w,
            -w,  w
        };
        double pts2[8] = {
            w, w,
            0, w,
            0,  -w,
            w,  -w
        };
        rotate_polygon(pts2, 4, test_angle);
        shift_polygon(pts2, 4, center + test_shift);

        material_function_for_polygons matfun(gv);
        polygon pol(pts1, 4, 2);
        pol.add_inner_polygon(polygon(pts2, 4, 2));
        matfun.add_polygon(pol, 12.0);
        pol = polygon(pts2, 4, 2);
        matfun.add_polygon(pol, 4.0);

        vec grad(matfun.normal_vector(icenter));
        double grad_ang = atan2(grad.y(), grad.x());
        double fratio(get_fill_ratio(test_angle, test_shift.x(), test_shift.y()));
        double test_meps = fratio * 12 + (1-fratio) * 4;
        double meps = matfun.mean_eps(icenter);
        double test_minveps = fratio / 12. + (1 - fratio) / 4.;
        double minveps = matfun.mean_inveps(icenter);


        if (fabs(test_angle - grad_ang) > 1e-8 ||
            fabs(test_minveps - minveps) > 1e-8 ||
            fabs(test_meps - meps) > 1e-8
        )
        {
            master_printf("error at #%i: test shift: (%.14f, %.14f) ; test angle: %.14f\n",
                          i, test_shift.x(), test_shift.y(), test_angle * 180 / pi);
            master_printf("expected angle: %.14f ; mean_eps: %.14f ; mean_inv_eps: %.14f\n",
                          test_angle * 180 / pi, test_meps, test_minveps);
            master_printf("     got angle: %.14f ; mean_eps: %.14f ; mean_inv_eps: %.14f\n\n",
                          grad_ang * 180 / pi, meps, minveps);
            abort("error");
        }
    }
}

/** Tests the normal direction, mean eps and mean 1/eps of a polygonal
 *  interface in 3D between 2 different material stacks at multiple randomly
 *  chosen angles and offsets. Also tests usage of inner polygons. */
void check_normals_3D(const int num_random_trials = 10000)
{
    master_printf("Checking interface of polygonal material with air (3D) - transversal components\n");
    const int w = 6;
    const int sh = int(w/2);
    grid_volume gv(vol3d(w, w, 4*sh, 1));
    volume vol(gv.surroundings());
    vec center = gv.center();
    ivec icenter = gv.icenter();
    double maxpd = 0, mpi = 0, mps = 0;
    vec mpshift = zero_vec(D3);

    for (int i = 0; i < num_random_trials; ++i) {
        double ptsA[8] = {
            -w, -w,
            0, -w,
            0,  w,
            -w,  w
        };
        double ptsB[8] = {
            w, w,
            0, w,
            0,  -w,
            w,  -w
        };

        double test_angle = urand(-pi, pi);
        vec test_shift = vec(urand(-0.5, 0.5), urand(-0.5, 0.5));
        double zs = urand(-0.5, 0.5);

        double eps = 12.0;
        double layer_thicknessA[3] = {sh + zs, 2 * sh, sh - zs};
        double layer_epsA[3] = {1.0, eps, 1.0};
        double layer_thicknessB[3] = {sh + zs, sh,  2 * sh - zs};
        double layer_epsB[3] = {1.0, eps, 1.0};
        material_function_for_polygons matfun(gv);
        unsigned int matstack_idA = matfun.add_material_stack(layer_thicknessA, layer_epsA, 3);
        unsigned int matstack_idB = matfun.add_material_stack(layer_thicknessB, layer_epsB, 3);

        rotate_polygon(ptsA, 4, test_angle);
        shift_polygon(ptsA, 4, center + test_shift);
        rotate_polygon(ptsB, 4, test_angle);
        shift_polygon(ptsB, 4, center + test_shift);
        polygon pol = polygon(ptsA, 4, 2);
        matfun.add_polygon(pol, matstack_idA);
        pol = polygon(ptsB, 4, 2);
        matfun.add_polygon(pol, matstack_idB);

        /////////////////////////////////////////////////////////////////
        // position at lower material interface where the two stacks meet:
        ivec ibottom = icenter - unit_ivec(D3, Z) * (2*sh);
        vec grad(matfun.normal_vector(ibottom));
        double fratio = 0.5 - zs;
        double test_meps = fratio * (eps - 1) + 1;
        double meps = matfun.mean_eps(ibottom);
        double test_minveps = 1 - fratio * (eps - 1) / eps;
        double minveps = matfun.mean_inveps(ibottom);
        // check that gradient is always pointing downwards at ibottom:
        if (fabs(grad.x()) > 1e-8 || fabs(grad.y()) > 1e-8 || grad.z() > -1e-8) {
            master_printf("error at #%i: test shift: (%.14f, %.14f, %.14f) ; test angle: %.14f\n",
                          i, test_shift.x(), test_shift.y(), zs, test_angle * 180 / pi);
            master_printf("gradient at %i, %i, %i (bottom): (%.14f, %.14f, %.14f) does not point downwards\n",
                          ibottom.x(), ibottom.y(), ibottom.z(), grad.x(), grad.y(), grad.z());
            abort("error");
        }
        // check mean epsilon and mean inverse epsilon:
        if (fabs(test_minveps - minveps) > 1e-8 ||
            fabs(test_meps - meps) > 1e-8
        )
        {
            master_printf("error at #%i: test shift: (%.14f, %.14f, %.14f) ; test angle: %.14f\n",
                          i, test_shift.x(), test_shift.y(), zs, test_angle * 180 / pi);
            master_printf("at %i, %i, %i (bottom):\n",
                          ibottom.x(), ibottom.y(), ibottom.z());
            master_printf("expected mean_eps: %.14f ; mean_inv_eps: %.14f\n",
                          test_meps, test_minveps);
            master_printf("     got mean_eps: %.14f ; mean_inv_eps: %.14f\n\n",
                          meps, minveps);
            abort("error");
        }

        /////////////////////////////////////////////////////////////////
        // position at upper material edge where there is only material in one stack:
        ivec itop = icenter + unit_ivec(D3, Z) * (2*sh);
        grad = matfun.normal_vector(itop);

        double phi_ang = atan2(grad.y(), grad.x());
        fratio = (0.5 + zs) * get_fill_ratio(
            test_angle, test_shift.x(), test_shift.y());
        test_meps = fratio * (eps - 1) + 1;
        meps = matfun.mean_eps(itop);
        test_minveps = 1 - fratio * (eps - 1) / eps;
        minveps = matfun.mean_inveps(itop);

        if (fabs(test_angle - phi_ang) > 1e-8 ||
            grad.z() < 1e-8 || // must point upwards
            fabs(test_minveps - minveps) > 1e-8 ||
            fabs(test_meps - meps) > 1e-8
        )
        {
            master_printf("error at #%i: test shift: (%.14f, %.14f, %.14f) ; test angle: %.14f\n",
                          i, test_shift.x(), test_shift.y(), zs, test_angle * 180 / pi);
            master_printf("gradient at %i, %i, %i (top): (%.14f, %.14f, %.16f)\n",
                          itop.x(), itop.y(), itop.z(), grad.x(), grad.y(), grad.z());
            master_printf("expected phi: %.14f ; mean_eps: %.14f ; mean_inv_eps: %.14f ; gradient.z > 0\n",
                          test_angle * 180 / pi, test_meps, test_minveps);
            master_printf("     got phi: %.14f ; mean_eps: %.14f ; mean_inv_eps: %.14f ; gradient.z %s 0\n\n",
                          phi_ang * 180 / pi, meps, minveps, grad.z() < 1e-12 ? "<=" : ">");
            abort("error");
        }

        /////////////////////////////////////////////////////////////////
        // center material edge, where the layers in the two stacks have different heights:
        grad = matfun.normal_vector(icenter);

        phi_ang = atan2(grad.y(), grad.x());
        fratio = 0.5 + zs + (0.5 - zs) * get_fill_ratio(
            test_angle, test_shift.x(), test_shift.y());
        test_meps = fratio * (eps - 1) + 1;
        meps = matfun.mean_eps(icenter);
        test_minveps = 1 - fratio * (eps - 1) / eps;
        minveps = matfun.mean_inveps(icenter);

        // grad is expected to point in opposite direction:
        if (fabs(test_angle - phi_ang) - pi > 1e-8 ||
            grad.z() < 1e-8 || // must point upwards
            fabs(test_minveps - minveps) > 1e-8 ||
            fabs(test_meps - meps) > 1e-8
        )
        {
            master_printf("error at #%i: test shift: (%.14f, %.14f, %.14f) ; test angle: %.14f\n",
                          i, test_shift.x(), test_shift.y(), zs, test_angle * 180 / pi);
            master_printf("gradient at %i, %i, %i (center): (%.14f, %.14f, %.16f)\n",
                          icenter.x(), icenter.y(), icenter.z(), grad.x(), grad.y(), grad.z());
            master_printf("expected phi: %.14f ; mean_eps: %.14f ; mean_inv_eps: %.14f ; gradient.z > 0\n",
                          test_angle * 180 / pi, test_meps, test_minveps);
            master_printf("     got phi: %.14f ; mean_eps: %.14f ; mean_inv_eps: %.14f ; gradient.z %s 0\n\n",
                          phi_ang * 180 / pi, meps, minveps, grad.z() < 1e-12 ? "<=" : ">");
            abort("error");
        }

    }
}

/** Tests the z component of a polygonal interface's normal direction in
 *  3D between 2 different material stacks for multiple randomly
 *  chosen angles and z-offsets. */
void check_z_normals_3D(const int num_random_trials = 10000)
{
    master_printf("Checking interface of polygonal material with air (3D) - z component\n");
    const int w = 6;
    const int sh = int(w/2);
    grid_volume gv(vol3d(w, w, 4*sh, 1));
    volume vol(gv.surroundings());
    vec center = gv.center();
    ivec icenter = gv.icenter();
    double maxpd = 0, mpi = 0, mps = 0;
    vec mpshift = zero_vec(D3);

    for (int i = 0; i < num_random_trials; ++i) {
        double ptsA[8] = {
            -w, -w,
            0, -w,
            0,  w,
            -w,  w
        };
        double ptsB[8] = {
            w, w,
            0, w,
            0,  -w,
            w,  -w
        };
        double test_angle = urand(-pi, pi);
        // in this test case, the lateral shift cannot be random, since we
        // don't know the analytical solution in those cases:
        vec test_shift = vec(0, 0);
        double zs = urand(-0.5, 0.5);

        double eps = 12.0;
        double layer_thicknessA[3] = {sh + zs, 2 * sh, sh - zs};
        double layer_epsA[3] = {1.0, eps, 1.0};
        double layer_thicknessB[3] = {sh + zs, sh,  2 * sh - zs};
        double layer_epsB[3] = {1.0, eps, 1.0};
        material_function_for_polygons matfun(gv);
        unsigned int matstack_idA = matfun.add_material_stack(layer_thicknessA, layer_epsA, 3);
        unsigned int matstack_idB = matfun.add_material_stack(layer_thicknessB, layer_epsB, 3);

        rotate_polygon(ptsA, 4, test_angle);
        shift_polygon(ptsA, 4, center + test_shift);
        rotate_polygon(ptsB, 4, test_angle);
        shift_polygon(ptsB, 4, center + test_shift);
        polygon pol = polygon(ptsA, 4, 2);
        matfun.add_polygon(pol, matstack_idA);
        pol = polygon(ptsB, 4, 2);
        matfun.add_polygon(pol, matstack_idB);

        // position at upper material edge where there is only material in one stack:
        ivec itop = icenter + unit_ivec(D3, Z) * (2*sh);// - unit_ivec(D3, X) * 1;
        vec grad = matfun.normal_vector(itop);
        // Also test at center position, where the material from two stacks meet:
        vec grad2 = matfun.normal_vector(icenter);
        // expected grad.z(), from analytical integration of -z*chi over the unit sphere:
        double exp_gz = (eps - 1) / 8.0 * (1 - 4 * pow(zs / sqrt(3), 2));
        // error, normalized with chi, so this test is independent of changes of eps:
        double err = fabs(exp_gz - grad.z()) / (eps - 1);
        double err2 = fabs(exp_gz - grad2.z()) / (eps - 1);
        // The Lebedev quadrature scheme used in this case in matfun.normal_vector returns
        // only values at discrete positions, since it cannot distinquish small shifts
        // in the structure position. So we cannot be too strict on the error, but
        // still keep a tight upper bound:
        double maxerr = 0.0158;
        if (err > maxerr || err2 > maxerr)
        {
            master_printf("error at #%i: test shift: (%.14f, %.14f, %.14f) ; test angle: %.14f\n",
                          i, test_shift.x(), test_shift.y(), zs, test_angle * 180 / pi);
            master_printf("gradient at %i, %i, %i: (%.14f, %.14f, %.14f), error: %f\n",
                          itop.x(), itop.y(), itop.z(), grad.x(), grad.y(), grad.z(), err);
            master_printf("gradient at %i, %i, %i: (%.14f, %.14f, %.14f), error: %f\n",
                          icenter.x(), icenter.y(), icenter.z(), grad2.x(), grad2.y(), grad2.z(), err2);
            master_printf("z-component not within allowed margin of expected analytical result: %f\n",
                          exp_gz);
            abort("error");
        }
    }
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  quiet = true;

  srand(271828); // use fixed random sequence

  check_normals_2D();
  check_normals_2D_2materials();
  check_normals_3D();
  check_z_normals_3D();

  return 0;
}
