#ifndef BICGSTAB_H
#define BICGSTAB_H

namespace meep {

typedef void (*bicgstab_op)(const double *x, double *y, void *data);

int bicgstab(int n, double *x,
             bicgstab_op A, void *Adata, const double *b,
             double tol, 
	     int *iters, // input *iters = max iters, output = actual iters
	     double *work); // if you pass work=NULL, bicgstab returns nwork

} // namespace meep

#endif /* BICGSTAB_H */
