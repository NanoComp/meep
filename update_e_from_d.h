if (pol) {
  double *d_minus_p = new double[ntot];
  if (e_sources) {
    FOR_E_AND_D(ec,dc) if (f[ec][0]) {
      DOCMP {
        for (int i=0;i<ntot;i++) {
          d_minus_p[i] = f[dc][cmp][i];
          for (polarization *p = pol; p; p = p->next) {
            d_minus_p[i] -= p->P[ec][cmp][i];
          }
        }
        for (src *s = e_sources; s; s = s->next) {
          d_minus_p[s->i] -= (cmp==0) ? real(s->get_dipole_now()*s->A[ec])
             : imag(s->get_dipole_now()*s->A[ec]);
        }
        for (int i=0;i<ntot;i++) {
          f[ec][cmp][i] = (ma->inveps[ec][component_direction(ec)][i])*d_minus_p[i];
        }
      }
    }
  } else {
    FOR_E_AND_D(ec,dc) if (f[ec][0]) {
      DOCMP {
        for (int i=0;i<ntot;i++) {
          d_minus_p[i] = f[dc][cmp][i];
          for (polarization *p = pol; p; p = p->next) {
            d_minus_p[i] -= p->P[ec][cmp][i];
          }
        }
        for (int i=0;i<ntot;i++) {
          f[ec][cmp][i] = (ma->inveps[ec][component_direction(ec)][i])*d_minus_p[i];
        }
      }
    }
  }
  delete[] d_minus_p;
} else {
  if (e_sources) {
    double *d_minus_p = new double[ntot];
    FOR_E_AND_D(ec,dc) if (f[ec][0]) {
      DOCMP {
        for (int i=0;i<ntot;i++) {
          d_minus_p[i] = f[dc][cmp][i];
        }
        for (src *s = e_sources; s; s = s->next) {
          d_minus_p[s->i] -= (cmp==0) ? real(s->get_dipole_now()*s->A[ec])
             : imag(s->get_dipole_now()*s->A[ec]);
        }
        for (int i=0;i<ntot;i++) {
          f[ec][cmp][i] = (ma->inveps[ec][component_direction(ec)][i])*d_minus_p[i];
        }
      }
    }
    delete[] d_minus_p;
  } else {
    FOR_E_AND_D(ec,dc) if (f[ec][0]) {
      DOCMP {
        for (int i=0;i<ntot;i++) {
          f[ec][cmp][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][cmp][i];
        }
      }
    }
  }
}
