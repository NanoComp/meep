if (pol) {
  double *d_minus_p = new double[ntot];
  if (e_sources) {
    if (is_real) {
      if (f[Er][0]) {
        // Am in cylindrical coordinates.
        FOR_E_AND_D(ec,dc) if (f[ec][0]) {
          const int yee_idx = v.yee_index(ec);
          const int d_ec = component_direction(ec);
          const int s_ec = stride_any_direction[d_ec];
          const direction d_1 = (direction)(((d_ec-2)+1)%3+2);
          const component ec_1 = direction_component(ec,d_1);
          const int s_1 = stride_any_direction[d_1];
          const direction d_2 = (direction)(((d_ec-2)+2)%3+2);
          const component ec_2 = direction_component(ec,d_2);
          const int s_2 = stride_any_direction[d_2];
          for (int i=0;i<ntot;i++) {
            for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
              np->energy[ec][i] = op->energy[ec][i];
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
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
      } else { // not f[Er][0]
        if (f[Ey][0]) {
          if (f[Ez][0]) {
            if (stride_any_direction[Z]) {
              // Am in 3D
              FOR_E_AND_D(ec,dc) if (f[ec][0]) {
                const int yee_idx = v.yee_index(ec);
                const int d_ec = component_direction(ec);
                const int s_ec = stride_any_direction[d_ec];
                const direction d_1 = (direction)((d_ec+1)%3);
                const component ec_1 = direction_component(ec,d_1);
                const int s_1 = stride_any_direction[d_1];
                const direction d_2 = (direction)((d_ec+2)%3);
                const component ec_2 = direction_component(ec,d_2);
                const int s_2 = stride_any_direction[d_2];
                for (int i=0;i<ntot;i++) {
                  for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                    np->energy[ec][i] = op->energy[ec][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
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
            } else { // not stride_any_direction[Z]
              // Am in 2D
              FOR_E_AND_D(ec,dc) if (f[ec][0]) {
                const int yee_idx = v.yee_index(ec);
                const int d_ec = component_direction(ec);
                const int s_ec = stride_any_direction[d_ec];
                const direction d_1 = (direction)((d_ec+1)%3);
                const component ec_1 = direction_component(ec,d_1);
                const int s_1 = stride_any_direction[d_1];
                const direction d_2 = (direction)((d_ec+2)%3);
                const component ec_2 = direction_component(ec,d_2);
                const int s_2 = stride_any_direction[d_2];
                for (int i=0;i<ntot;i++) {
                  for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                    np->energy[ec][i] = op->energy[ec][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
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
            }
          } else { // not f[Ez][0]
            // Am in 2D TE
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              const direction d_1 = (direction)((d_ec+1)%2);
              const component ec_1 = direction_component(ec,d_1);
              const int s_1 = stride_any_direction[d_1];
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
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
          }
        } else { // not f[Ey][0]
          if (f[Ex][0]) {
            // Am in 1D
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
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
          } else { // not f[Ex][0]
            // Am in 2D TM
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
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
          }
        }
      }
    } else { // not is_real
      if (f[Er][0]) {
        // Am in cylindrical coordinates.
        FOR_E_AND_D(ec,dc) if (f[ec][0]) {
          const int yee_idx = v.yee_index(ec);
          const int d_ec = component_direction(ec);
          const int s_ec = stride_any_direction[d_ec];
          const direction d_1 = (direction)(((d_ec-2)+1)%3+2);
          const component ec_1 = direction_component(ec,d_1);
          const int s_1 = stride_any_direction[d_1];
          const direction d_2 = (direction)(((d_ec-2)+2)%3+2);
          const component ec_2 = direction_component(ec,d_2);
          const int s_2 = stride_any_direction[d_2];
          for (int i=0;i<ntot;i++) {
            for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
              np->energy[ec][i] = op->energy[ec][i];
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                 * f[ec][0][i];
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
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
      } else { // not f[Er][0]
        if (f[Ey][0]) {
          if (f[Ez][0]) {
            if (stride_any_direction[Z]) {
              // Am in 3D
              FOR_E_AND_D(ec,dc) if (f[ec][0]) {
                const int yee_idx = v.yee_index(ec);
                const int d_ec = component_direction(ec);
                const int s_ec = stride_any_direction[d_ec];
                const direction d_1 = (direction)((d_ec+1)%3);
                const component ec_1 = direction_component(ec,d_1);
                const int s_1 = stride_any_direction[d_1];
                const direction d_2 = (direction)((d_ec+2)%3);
                const component ec_2 = direction_component(ec,d_2);
                const int s_2 = stride_any_direction[d_2];
                for (int i=0;i<ntot;i++) {
                  for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                    np->energy[ec][i] = op->energy[ec][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                       * f[ec][0][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
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
            } else { // not stride_any_direction[Z]
              // Am in 2D
              FOR_E_AND_D(ec,dc) if (f[ec][0]) {
                const int yee_idx = v.yee_index(ec);
                const int d_ec = component_direction(ec);
                const int s_ec = stride_any_direction[d_ec];
                const direction d_1 = (direction)((d_ec+1)%3);
                const component ec_1 = direction_component(ec,d_1);
                const int s_1 = stride_any_direction[d_1];
                const direction d_2 = (direction)((d_ec+2)%3);
                const component ec_2 = direction_component(ec,d_2);
                const int s_2 = stride_any_direction[d_2];
                for (int i=0;i<ntot;i++) {
                  for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                    np->energy[ec][i] = op->energy[ec][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                       * f[ec][0][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
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
          } else { // not f[Ez][0]
            // Am in 2D TE
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              const direction d_1 = (direction)((d_ec+1)%2);
              const component ec_1 = direction_component(ec,d_1);
              const int s_1 = stride_any_direction[d_1];
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
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
        } else { // not f[Ey][0]
          if (f[Ex][0]) {
            // Am in 1D
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
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
          } else { // not f[Ex][0]
            // Am in 2D TM
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
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
        }
      }
    }
  } else { // not e_sources
    if (is_real) {
      if (f[Er][0]) {
        // Am in cylindrical coordinates.
        FOR_E_AND_D(ec,dc) if (f[ec][0]) {
          const int yee_idx = v.yee_index(ec);
          const int d_ec = component_direction(ec);
          const int s_ec = stride_any_direction[d_ec];
          const direction d_1 = (direction)(((d_ec-2)+1)%3+2);
          const component ec_1 = direction_component(ec,d_1);
          const int s_1 = stride_any_direction[d_1];
          const direction d_2 = (direction)(((d_ec-2)+2)%3+2);
          const component ec_2 = direction_component(ec,d_2);
          const int s_2 = stride_any_direction[d_2];
          for (int i=0;i<ntot;i++) {
            for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
              np->energy[ec][i] = op->energy[ec][i];
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
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
      } else { // not f[Er][0]
        if (f[Ey][0]) {
          if (f[Ez][0]) {
            if (stride_any_direction[Z]) {
              // Am in 3D
              FOR_E_AND_D(ec,dc) if (f[ec][0]) {
                const int yee_idx = v.yee_index(ec);
                const int d_ec = component_direction(ec);
                const int s_ec = stride_any_direction[d_ec];
                const direction d_1 = (direction)((d_ec+1)%3);
                const component ec_1 = direction_component(ec,d_1);
                const int s_1 = stride_any_direction[d_1];
                const direction d_2 = (direction)((d_ec+2)%3);
                const component ec_2 = direction_component(ec,d_2);
                const int s_2 = stride_any_direction[d_2];
                for (int i=0;i<ntot;i++) {
                  for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                    np->energy[ec][i] = op->energy[ec][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
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
            } else { // not stride_any_direction[Z]
              // Am in 2D
              FOR_E_AND_D(ec,dc) if (f[ec][0]) {
                const int yee_idx = v.yee_index(ec);
                const int d_ec = component_direction(ec);
                const int s_ec = stride_any_direction[d_ec];
                const direction d_1 = (direction)((d_ec+1)%3);
                const component ec_1 = direction_component(ec,d_1);
                const int s_1 = stride_any_direction[d_1];
                const direction d_2 = (direction)((d_ec+2)%3);
                const component ec_2 = direction_component(ec,d_2);
                const int s_2 = stride_any_direction[d_2];
                for (int i=0;i<ntot;i++) {
                  for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                    np->energy[ec][i] = op->energy[ec][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
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
            }
          } else { // not f[Ez][0]
            // Am in 2D TE
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              const direction d_1 = (direction)((d_ec+1)%2);
              const component ec_1 = direction_component(ec,d_1);
              const int s_1 = stride_any_direction[d_1];
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
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
          }
        } else { // not f[Ey][0]
          if (f[Ex][0]) {
            // Am in 1D
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
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
          } else { // not f[Ex][0]
            // Am in 2D TM
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
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
          }
        }
      }
    } else { // not is_real
      if (f[Er][0]) {
        // Am in cylindrical coordinates.
        FOR_E_AND_D(ec,dc) if (f[ec][0]) {
          const int yee_idx = v.yee_index(ec);
          const int d_ec = component_direction(ec);
          const int s_ec = stride_any_direction[d_ec];
          const direction d_1 = (direction)(((d_ec-2)+1)%3+2);
          const component ec_1 = direction_component(ec,d_1);
          const int s_1 = stride_any_direction[d_1];
          const direction d_2 = (direction)(((d_ec-2)+2)%3+2);
          const component ec_2 = direction_component(ec,d_2);
          const int s_2 = stride_any_direction[d_2];
          for (int i=0;i<ntot;i++) {
            for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
              np->energy[ec][i] = op->energy[ec][i];
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                 * f[ec][0][i];
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
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
      } else { // not f[Er][0]
        if (f[Ey][0]) {
          if (f[Ez][0]) {
            if (stride_any_direction[Z]) {
              // Am in 3D
              FOR_E_AND_D(ec,dc) if (f[ec][0]) {
                const int yee_idx = v.yee_index(ec);
                const int d_ec = component_direction(ec);
                const int s_ec = stride_any_direction[d_ec];
                const direction d_1 = (direction)((d_ec+1)%3);
                const component ec_1 = direction_component(ec,d_1);
                const int s_1 = stride_any_direction[d_1];
                const direction d_2 = (direction)((d_ec+2)%3);
                const component ec_2 = direction_component(ec,d_2);
                const int s_2 = stride_any_direction[d_2];
                for (int i=0;i<ntot;i++) {
                  for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                    np->energy[ec][i] = op->energy[ec][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                       * f[ec][0][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
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
            } else { // not stride_any_direction[Z]
              // Am in 2D
              FOR_E_AND_D(ec,dc) if (f[ec][0]) {
                const int yee_idx = v.yee_index(ec);
                const int d_ec = component_direction(ec);
                const int s_ec = stride_any_direction[d_ec];
                const direction d_1 = (direction)((d_ec+1)%3);
                const component ec_1 = direction_component(ec,d_1);
                const int s_1 = stride_any_direction[d_1];
                const direction d_2 = (direction)((d_ec+2)%3);
                const component ec_2 = direction_component(ec,d_2);
                const int s_2 = stride_any_direction[d_2];
                for (int i=0;i<ntot;i++) {
                  for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                    np->energy[ec][i] = op->energy[ec][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                       * f[ec][0][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
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
          } else { // not f[Ez][0]
            // Am in 2D TE
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              const direction d_1 = (direction)((d_ec+1)%2);
              const component ec_1 = direction_component(ec,d_1);
              const int s_1 = stride_any_direction[d_1];
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
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
        } else { // not f[Ey][0]
          if (f[Ex][0]) {
            // Am in 1D
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
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
          } else { // not f[Ex][0]
            // Am in 2D TM
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
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
      }
    }
  }
  delete[] d_minus_p;
} else { // not pol
  if (e_sources) {
    double *d_minus_p = new double[ntot];
    if (is_real) {
      if (f[Er][0]) {
        // Am in cylindrical coordinates.
        FOR_E_AND_D(ec,dc) if (f[ec][0]) {
          const int yee_idx = v.yee_index(ec);
          const int d_ec = component_direction(ec);
          const int s_ec = stride_any_direction[d_ec];
          const direction d_1 = (direction)(((d_ec-2)+1)%3+2);
          const component ec_1 = direction_component(ec,d_1);
          const int s_1 = stride_any_direction[d_1];
          const direction d_2 = (direction)(((d_ec-2)+2)%3+2);
          const component ec_2 = direction_component(ec,d_2);
          const int s_2 = stride_any_direction[d_2];
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
      } else { // not f[Er][0]
        if (f[Ey][0]) {
          if (f[Ez][0]) {
            if (stride_any_direction[Z]) {
              // Am in 3D
              FOR_E_AND_D(ec,dc) if (f[ec][0]) {
                const int yee_idx = v.yee_index(ec);
                const int d_ec = component_direction(ec);
                const int s_ec = stride_any_direction[d_ec];
                const direction d_1 = (direction)((d_ec+1)%3);
                const component ec_1 = direction_component(ec,d_1);
                const int s_1 = stride_any_direction[d_1];
                const direction d_2 = (direction)((d_ec+2)%3);
                const component ec_2 = direction_component(ec,d_2);
                const int s_2 = stride_any_direction[d_2];
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
            } else { // not stride_any_direction[Z]
              // Am in 2D
              FOR_E_AND_D(ec,dc) if (f[ec][0]) {
                const int yee_idx = v.yee_index(ec);
                const int d_ec = component_direction(ec);
                const int s_ec = stride_any_direction[d_ec];
                const direction d_1 = (direction)((d_ec+1)%3);
                const component ec_1 = direction_component(ec,d_1);
                const int s_1 = stride_any_direction[d_1];
                const direction d_2 = (direction)((d_ec+2)%3);
                const component ec_2 = direction_component(ec,d_2);
                const int s_2 = stride_any_direction[d_2];
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
            }
          } else { // not f[Ez][0]
            // Am in 2D TE
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              const direction d_1 = (direction)((d_ec+1)%2);
              const component ec_1 = direction_component(ec,d_1);
              const int s_1 = stride_any_direction[d_1];
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
          }
        } else { // not f[Ey][0]
          if (f[Ex][0]) {
            // Am in 1D
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
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
          } else { // not f[Ex][0]
            // Am in 2D TM
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
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
          }
        }
      }
    } else { // not is_real
      if (f[Er][0]) {
        // Am in cylindrical coordinates.
        FOR_E_AND_D(ec,dc) if (f[ec][0]) {
          const int yee_idx = v.yee_index(ec);
          const int d_ec = component_direction(ec);
          const int s_ec = stride_any_direction[d_ec];
          const direction d_1 = (direction)(((d_ec-2)+1)%3+2);
          const component ec_1 = direction_component(ec,d_1);
          const int s_1 = stride_any_direction[d_1];
          const direction d_2 = (direction)(((d_ec-2)+2)%3+2);
          const component ec_2 = direction_component(ec,d_2);
          const int s_2 = stride_any_direction[d_2];
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
      } else { // not f[Er][0]
        if (f[Ey][0]) {
          if (f[Ez][0]) {
            if (stride_any_direction[Z]) {
              // Am in 3D
              FOR_E_AND_D(ec,dc) if (f[ec][0]) {
                const int yee_idx = v.yee_index(ec);
                const int d_ec = component_direction(ec);
                const int s_ec = stride_any_direction[d_ec];
                const direction d_1 = (direction)((d_ec+1)%3);
                const component ec_1 = direction_component(ec,d_1);
                const int s_1 = stride_any_direction[d_1];
                const direction d_2 = (direction)((d_ec+2)%3);
                const component ec_2 = direction_component(ec,d_2);
                const int s_2 = stride_any_direction[d_2];
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
            } else { // not stride_any_direction[Z]
              // Am in 2D
              FOR_E_AND_D(ec,dc) if (f[ec][0]) {
                const int yee_idx = v.yee_index(ec);
                const int d_ec = component_direction(ec);
                const int s_ec = stride_any_direction[d_ec];
                const direction d_1 = (direction)((d_ec+1)%3);
                const component ec_1 = direction_component(ec,d_1);
                const int s_1 = stride_any_direction[d_1];
                const direction d_2 = (direction)((d_ec+2)%3);
                const component ec_2 = direction_component(ec,d_2);
                const int s_2 = stride_any_direction[d_2];
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
          } else { // not f[Ez][0]
            // Am in 2D TE
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              const direction d_1 = (direction)((d_ec+1)%2);
              const component ec_1 = direction_component(ec,d_1);
              const int s_1 = stride_any_direction[d_1];
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
        } else { // not f[Ey][0]
          if (f[Ex][0]) {
            // Am in 1D
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
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
          } else { // not f[Ex][0]
            // Am in 2D TM
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
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
        }
      }
    }
    delete[] d_minus_p;
  } else { // not e_sources
    if (is_real) {
      if (f[Er][0]) {
        // Am in cylindrical coordinates.
        FOR_E_AND_D(ec,dc) if (f[ec][0]) {
          const int yee_idx = v.yee_index(ec);
          const int d_ec = component_direction(ec);
          const int s_ec = stride_any_direction[d_ec];
          const direction d_1 = (direction)(((d_ec-2)+1)%3+2);
          const component ec_1 = direction_component(ec,d_1);
          const int s_1 = stride_any_direction[d_1];
          const direction d_2 = (direction)(((d_ec-2)+2)%3+2);
          const component ec_2 = direction_component(ec,d_2);
          const int s_2 = stride_any_direction[d_2];
          for (int i=0;i<ntot;i++) {
            f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][0][i];
          }
        }
      } else { // not f[Er][0]
        if (f[Ey][0]) {
          if (f[Ez][0]) {
            if (stride_any_direction[Z]) {
              // Am in 3D
              FOR_E_AND_D(ec,dc) if (f[ec][0]) {
                const int yee_idx = v.yee_index(ec);
                const int d_ec = component_direction(ec);
                const int s_ec = stride_any_direction[d_ec];
                const direction d_1 = (direction)((d_ec+1)%3);
                const component ec_1 = direction_component(ec,d_1);
                const int s_1 = stride_any_direction[d_1];
                const direction d_2 = (direction)((d_ec+2)%3);
                const component ec_2 = direction_component(ec,d_2);
                const int s_2 = stride_any_direction[d_2];
                for (int i=0;i<ntot;i++) {
                  f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][0][i];
                }
              }
            } else { // not stride_any_direction[Z]
              // Am in 2D
              FOR_E_AND_D(ec,dc) if (f[ec][0]) {
                const int yee_idx = v.yee_index(ec);
                const int d_ec = component_direction(ec);
                const int s_ec = stride_any_direction[d_ec];
                const direction d_1 = (direction)((d_ec+1)%3);
                const component ec_1 = direction_component(ec,d_1);
                const int s_1 = stride_any_direction[d_1];
                const direction d_2 = (direction)((d_ec+2)%3);
                const component ec_2 = direction_component(ec,d_2);
                const int s_2 = stride_any_direction[d_2];
                for (int i=0;i<ntot;i++) {
                  f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][0][i];
                }
              }
            }
          } else { // not f[Ez][0]
            // Am in 2D TE
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              const direction d_1 = (direction)((d_ec+1)%2);
              const component ec_1 = direction_component(ec,d_1);
              const int s_1 = stride_any_direction[d_1];
              for (int i=0;i<ntot;i++) {
                f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][0][i];
              }
            }
          }
        } else { // not f[Ey][0]
          if (f[Ex][0]) {
            // Am in 1D
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              for (int i=0;i<ntot;i++) {
                f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][0][i];
              }
            }
          } else { // not f[Ex][0]
            // Am in 2D TM
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              for (int i=0;i<ntot;i++) {
                f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][0][i];
              }
            }
          }
        }
      }
    } else { // not is_real
      if (f[Er][0]) {
        // Am in cylindrical coordinates.
        FOR_E_AND_D(ec,dc) if (f[ec][0]) {
          const int yee_idx = v.yee_index(ec);
          const int d_ec = component_direction(ec);
          const int s_ec = stride_any_direction[d_ec];
          const direction d_1 = (direction)(((d_ec-2)+1)%3+2);
          const component ec_1 = direction_component(ec,d_1);
          const int s_1 = stride_any_direction[d_1];
          const direction d_2 = (direction)(((d_ec-2)+2)%3+2);
          const component ec_2 = direction_component(ec,d_2);
          const int s_2 = stride_any_direction[d_2];
          for (int i=0;i<ntot;i++) {
            f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][0][i];
          }
          for (int i=0;i<ntot;i++) {
            f[ec][1][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][1][i];
          }
        }
      } else { // not f[Er][0]
        if (f[Ey][0]) {
          if (f[Ez][0]) {
            if (stride_any_direction[Z]) {
              // Am in 3D
              FOR_E_AND_D(ec,dc) if (f[ec][0]) {
                const int yee_idx = v.yee_index(ec);
                const int d_ec = component_direction(ec);
                const int s_ec = stride_any_direction[d_ec];
                const direction d_1 = (direction)((d_ec+1)%3);
                const component ec_1 = direction_component(ec,d_1);
                const int s_1 = stride_any_direction[d_1];
                const direction d_2 = (direction)((d_ec+2)%3);
                const component ec_2 = direction_component(ec,d_2);
                const int s_2 = stride_any_direction[d_2];
                for (int i=0;i<ntot;i++) {
                  f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][0][i];
                }
                for (int i=0;i<ntot;i++) {
                  f[ec][1][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][1][i];
                }
              }
            } else { // not stride_any_direction[Z]
              // Am in 2D
              FOR_E_AND_D(ec,dc) if (f[ec][0]) {
                const int yee_idx = v.yee_index(ec);
                const int d_ec = component_direction(ec);
                const int s_ec = stride_any_direction[d_ec];
                const direction d_1 = (direction)((d_ec+1)%3);
                const component ec_1 = direction_component(ec,d_1);
                const int s_1 = stride_any_direction[d_1];
                const direction d_2 = (direction)((d_ec+2)%3);
                const component ec_2 = direction_component(ec,d_2);
                const int s_2 = stride_any_direction[d_2];
                for (int i=0;i<ntot;i++) {
                  f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][0][i];
                }
                for (int i=0;i<ntot;i++) {
                  f[ec][1][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][1][i];
                }
              }
            }
          } else { // not f[Ez][0]
            // Am in 2D TE
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              const direction d_1 = (direction)((d_ec+1)%2);
              const component ec_1 = direction_component(ec,d_1);
              const int s_1 = stride_any_direction[d_1];
              for (int i=0;i<ntot;i++) {
                f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][0][i];
              }
              for (int i=0;i<ntot;i++) {
                f[ec][1][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][1][i];
              }
            }
          }
        } else { // not f[Ey][0]
          if (f[Ex][0]) {
            // Am in 1D
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
              for (int i=0;i<ntot;i++) {
                f[ec][0][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][0][i];
              }
              for (int i=0;i<ntot;i++) {
                f[ec][1][i] = (ma->inveps[ec][component_direction(ec)][i])*f[dc][1][i];
              }
            }
          } else { // not f[Ex][0]
            // Am in 2D TM
            FOR_E_AND_D(ec,dc) if (f[ec][0]) {
              const int yee_idx = v.yee_index(ec);
              const int d_ec = component_direction(ec);
              const int s_ec = stride_any_direction[d_ec];
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
    }
  }
}
