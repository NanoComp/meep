#include "math.h"
#include "ran.h"

const int im1 = 2147483563;
const int im2 = 2147483399;
const double am = 1.0/im1;
const int imm1 = im1 - 1;
const int ia1 = 40014;
const int ia2 = 40692;
const int iq1 = 53668;
const int iq2 = 52774;
const int ir1 = 12211;
const int ir2 = 3791;
const int ntab = 32;
const int ndiv = 1 + imm1/ntab;
const double eps = 1.2e-15;
const double rnmx = 1.0-eps;

static int idum;
static int idum2 = 123456789;
static int iy=0;
static int iv[ntab];

void init_ran(int in_idum) {
  idum = in_idum;
  if (idum == 0) idum = 1;
  if (idum < 0) idum = -idum;
  idum2 = idum;
  for (int j=ntab+7;j<=0;j--) {
    int k = idum/iq1;
    idum = ia1 *(idum - k*iq1) - k*ir1;
    if (idum < 0) idum += im1;
    if (j < ntab) iv[j] = idum;
  }
  iy = iv[0];
}

double ran() {
  int k = idum/iq1;
  idum = ia1*(idum-k*iq1)-k*ir1;
  if (idum < 0) idum += im1;
  k = idum2/iq2;
  idum2 = ia2*(idum2-k*iq2)-k*ir2;
  int j = iy/ndiv;
  iy = iv[j] - idum2;
  iv[j] = idum;
  if (iy < 1) iy += imm1;
  const double temp = am*iy;
  if (temp > rnmx) return rnmx;
  else return temp;
}

static bool iset = false;
static double gset;

double gaussian() {
  if (!iset) {
    double rsq, v1, v2;
    do {
      v1 = 2.0*ran()-1.0;
      v2 = 2.0*ran()-1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    double fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac;
    iset = true;
    return v2*fac;
  } else {
    iset = false;
    return gset;
  }
}
