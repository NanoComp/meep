if (pol) {
  double *d_minus_p = new double[ntot];
  if (e_sources) {
    if (is_real) {
      FOR_E_AND_D(ec,dc) if (f[ec][0]) {
        for (int i=0;i<ntot;i++) {
          for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
            np->energy[ec][i] = op->energy[ec][i];
            np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])
               * f[ec][0][i];
          }
        }
        for (int i=0;i<ntot;i++) {
          d_minus_p[i] = f[dc][0][i];
          for (polarization *p = pol; p; p = p->next) {
            d_minus_p[i] -= p->P[ec][0][i];
          }
        }
        for (src *s = e_sources; s; s = s->next) {
          d_minus_p[s->i] -= real(s->get_dipole_now()*s->A[ec]);
        }
        for (int i=0;i<ntot;i++) {
          f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*d_minus_p[i];
        }
      }
    } else { // not is_real
      FOR_E_AND_D(ec,dc) if (f[ec][0]) {
        for (int i=0;i<ntot;i++) {
          for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
            np->energy[ec][i] = op->energy[ec][i];
            np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])
               * f[ec][0][i];
            np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])
               * f[ec][1][i];
          }
        }
        for (int i=0;i<ntot;i++) {
          d_minus_p[i] = f[dc][0][i];
          for (polarization *p = pol; p; p = p->next) {
            d_minus_p[i] -= p->P[ec][0][i];
          }
        }
        for (src *s = e_sources; s; s = s->next) {
          d_minus_p[s->i] -= real(s->get_dipole_now()*s->A[ec]);
        }
        for (int i=0;i<ntot;i++) {
          f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*d_minus_p[i];
        }
        for (int i=0;i<ntot;i++) {
          d_minus_p[i] = f[dc][1][i];
          for (polarization *p = pol; p; p = p->next) {
            d_minus_p[i] -= p->P[ec][1][i];
          }
        }
        for (src *s = e_sources; s; s = s->next) {
          d_minus_p[s->i] -= imag(s->get_dipole_now()*s->A[ec]);
        }
        for (int i=0;i<ntot;i++) {
          f[ec][1][i] = (ma->inveps[ec][component_direction(ec)][i])*d_minus_p[i];
        }
      }
    }
  } else { // not e_sources
    if (is_real) {
      FOR_E_AND_D(ec,dc) if (f[ec][0]) {
        for (int i=0;i<ntot;i++) {
          for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
            np->energy[ec][i] = op->energy[ec][i];
            np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])
               * f[ec][0][i];
          }
        }
        for (int i=0;i<ntot;i++) {
          d_minus_p[i] = f[dc][0][i];
          for (polarization *p = pol; p; p = p->next) {
            d_minus_p[i] -= p->P[ec][0][i];
          }
        }
        for (int i=0;i<ntot;i++) {
          f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*d_minus_p[i];
        }
      }
    } else { // not is_real
      FOR_E_AND_D(ec,dc) if (f[ec][0]) {
        for (int i=0;i<ntot;i++) {
          for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
            np->energy[ec][i] = op->energy[ec][i];
            np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])
               * f[ec][0][i];
            np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])
               * f[ec][1][i];
          }
        }
        for (int i=0;i<ntot;i++) {
          d_minus_p[i] = f[dc][0][i];
          for (polarization *p = pol; p; p = p->next) {
            d_minus_p[i] -= p->P[ec][0][i];
          }
        }
        for (int i=0;i<ntot;i++) {
          f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*d_minus_p[i];
        }
        for (int i=0;i<ntot;i++) {
          d_minus_p[i] = f[dc][1][i];
          for (polarization *p = pol; p; p = p->next) {
            d_minus_p[i] -= p->P[ec][1][i];
          }
        }
        for (int i=0;i<ntot;i++) {
          f[ec][1][i] = (ma->inveps[ec][component_direction(ec)][i])*d_minus_p[i];
        }
      }
    }
  }
  delete[] d_minus_p;
} else { // not pol
  if (e_sources) {
    double *d_minus_p = new double[ntot];
    if (is_real) {
      FOR_E_AND_D(ec,dc) if (f[ec][0]) {
        for (int i=0;i<ntot;i++) {
          d_minus_p[i] = f[dc][0][i];
        }
        for (src *s = e_sources; s; s = s->next) {
          d_minus_p[s->i] -= real(s->get_dipole_now()*s->A[ec]);
        }
        for (int i=0;i<ntot;i++) {
          f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*d_minus_p[i];
        }
      }
    } else { // not is_real
      FOR_E_AND_D(ec,dc) if (f[ec][0]) {
        for (int i=0;i<ntot;i++) {
          d_minus_p[i] = f[dc][0][i];
        }
        for (src *s = e_sources; s; s = s->next) {
          d_minus_p[s->i] -= real(s->get_dipole_now()*s->A[ec]);
        }
        for (int i=0;i<ntot;i++) {
          f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*d_minus_p[i];
        }
        for (int i=0;i<ntot;i++) {
          d_minus_p[i] = f[dc][1][i];
        }
        for (src *s = e_sources; s; s = s->next) {
          d_minus_p[s->i] -= imag(s->get_dipole_now()*s->A[ec]);
        }
        for (int i=0;i<ntot;i++) {
          f[ec][1][i] = (ma->inveps[ec][component_direction(ec)][i])*d_minus_p[i];
        }
      }
    }
    delete[] d_minus_p;
  } else { // not e_sources
    if (is_real) {
      FOR_E_AND_D(ec,dc) if (f[ec][0]) {
        for (int i=0;i<ntot;i++) {
          f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][0][i];
        }
      }
    } else { // not is_real
      FOR_E_AND_D(ec,dc) if (f[ec][0]) {
        for (int i=0;i<ntot;i++) {
          f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][0][i];
        }
        for (int i=0;i<ntot;i++) {
          f[ec][1][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][1][i];
        }
      }
    }
  }
}
