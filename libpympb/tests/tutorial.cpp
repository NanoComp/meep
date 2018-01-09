#include "pympb.hpp"

int main() {

    int band_num = 8;
    vector3 kpoint = {0.0, 0.0, 0.0};
    bool match_frequency = true;
    int parity = 0;
    double resolution = 32;
    double eigensolver_tol = 1.0e-7;

    py_mpb::add_eigenmode_source(band_num, kpoint, match_frequency, parity, resolution, eigensolver_tol);

    return 0;
}