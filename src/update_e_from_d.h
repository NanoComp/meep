if (pol) {
  double *d_minus_p[5][2];
  DOCMP {
    FOR_ELECTRIC_COMPONENTS(ec) {
      if (f[ec]) {
        d_minus_p[ec][cmp] = new double[v.ntot()];
      } else { // not f[ec]
        d_minus_p[ec][cmp] = 0;
      }
    }
  }
  if (e_sources) {
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
        if (is_real) {
          for (int i=0;i<ntot;i++) {
            for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
              np->energy[ec][i] = op->energy[ec][i];
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                 * f[ec][0][i];
            }
          }
        } else { // not is_real
          for (int i=0;i<ntot;i++) {
            for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
              np->energy[ec][i] = op->energy[ec][i];
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                 * f[ec][0][i];
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
                 * f[ec][1][i];
            }
          }
        }
        DOCMP {
          for (int i=0;i<ntot;i++) {
            d_minus_p[ec][cmp][i] = f[dc][cmp][i];
            for (polarization *p = pol; p; p = p->next) {
              d_minus_p[ec][cmp][i] -= p->P[ec][cmp][i];
            }
          }
          for (src *s = e_sources; s; s = s->next) {
            d_minus_p[ec][cmp][s->i] -= (cmp) ? imag(s->get_dipole_now()*s->A[ec])
               : real(s->get_dipole_now()*s->A[ec]);
          }
        }
      }
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
        if (ma->inveps[ec][d_1]) {
          if (ma->inveps[ec][d_2]) {
            DOCMP {
              if (num_any_direction[Z]==1) {
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                     + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                        + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                           + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                              + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                 + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                       + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                          + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                             + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                   + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR++) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                         + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                            + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                               + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                  + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                     + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                    }
                  }
                }
              }
            }
          } else { // not ma->inveps[ec][d_2]
            DOCMP {
              if (num_any_direction[Z]==1) {
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                     + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                        + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                           + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                       + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                          + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                             + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR++) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                         + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                            + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                               + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                    }
                  }
                }
              }
            }
          }
        } else { // not ma->inveps[ec][d_1]
          if (ma->inveps[ec][d_2]) {
            DOCMP {
              if (num_any_direction[Z]==1) {
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                     + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                        + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                           + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                       + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                          + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                             + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR++) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                         + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                            + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                               + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                    }
                  }
                }
              }
            }
          } else { // not ma->inveps[ec][d_2]
            DOCMP {
              if (num_any_direction[Z]==1) {
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR++) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                    }
                  }
                }
              }
            }
          }
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
              if (is_real) {
                for (int i=0;i<ntot;i++) {
                  for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                    np->energy[ec][i] = op->energy[ec][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                       * f[ec][0][i];
                  }
                }
              } else { // not is_real
                for (int i=0;i<ntot;i++) {
                  for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                    np->energy[ec][i] = op->energy[ec][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                       * f[ec][0][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
                       * f[ec][1][i];
                  }
                }
              }
              DOCMP {
                for (int i=0;i<ntot;i++) {
                  d_minus_p[ec][cmp][i] = f[dc][cmp][i];
                  for (polarization *p = pol; p; p = p->next) {
                    d_minus_p[ec][cmp][i] -= p->P[ec][cmp][i];
                  }
                }
                for (src *s = e_sources; s; s = s->next) {
                  d_minus_p[ec][cmp][s->i] -= (cmp) ? imag(s->get_dipole_now()*s->A[ec])
                     : real(s->get_dipole_now()*s->A[ec]);
                }
              }
            }
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
              if (ma->inveps[ec][d_1]) {
                if (ma->inveps[ec][d_2]) {
                  DOCMP {
                    if (num_any_direction[Z]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                   + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                      + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                         + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                  + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                     + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                        + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                           + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                 + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                    + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                       + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                          + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                             + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iY=0; iY<num_any_direction[Y]; iY++) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                   + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                      + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                         + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                            + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                               + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                } else { // not ma->inveps[ec][d_2]
                  DOCMP {
                    if (num_any_direction[Z]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                   + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                  + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                     + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                 + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                    + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                       + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iY=0; iY<num_any_direction[Y]; iY++) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                   + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                      + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                         + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              } else { // not ma->inveps[ec][d_1]
                if (ma->inveps[ec][d_2]) {
                  DOCMP {
                    if (num_any_direction[Z]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                   + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                  + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                     + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                 + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                    + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                       + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iY=0; iY<num_any_direction[Y]; iY++) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                   + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                      + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                         + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                } else { // not ma->inveps[ec][d_2]
                  DOCMP {
                    if (num_any_direction[Z]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                          f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iY=0; iY<num_any_direction[Y]; iY++) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
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
              if (is_real) {
                for (int i=0;i<ntot;i++) {
                  for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                    np->energy[ec][i] = op->energy[ec][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                       * f[ec][0][i];
                  }
                }
              } else { // not is_real
                for (int i=0;i<ntot;i++) {
                  for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                    np->energy[ec][i] = op->energy[ec][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                       * f[ec][0][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
                       * f[ec][1][i];
                  }
                }
              }
              DOCMP {
                for (int i=0;i<ntot;i++) {
                  d_minus_p[ec][cmp][i] = f[dc][cmp][i];
                  for (polarization *p = pol; p; p = p->next) {
                    d_minus_p[ec][cmp][i] -= p->P[ec][cmp][i];
                  }
                }
                for (src *s = e_sources; s; s = s->next) {
                  d_minus_p[ec][cmp][s->i] -= (cmp) ? imag(s->get_dipole_now()*s->A[ec])
                     : real(s->get_dipole_now()*s->A[ec]);
                }
              }
            }
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
              if (ma->inveps[ec][d_1]) {
                if (ma->inveps[ec][d_2]) {
                  DOCMP {
                    if (num_any_direction[X]==1) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                           + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                              + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                 + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                    + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                       + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                   + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                      + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                         + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                  + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                     + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                        + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                           + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                          }
                        }
                      }
                    }
                  }
                } else { // not ma->inveps[ec][d_2]
                  DOCMP {
                    if (num_any_direction[X]==1) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                           + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                              + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                 + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                   + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                  + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                     + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                          }
                        }
                      }
                    }
                  }
                }
              } else { // not ma->inveps[ec][d_1]
                if (ma->inveps[ec][d_2]) {
                  DOCMP {
                    if (num_any_direction[X]==1) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                           + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                              + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                 + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                   + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                  + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                     + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                          }
                        }
                      }
                    }
                  }
                } else { // not ma->inveps[ec][d_2]
                  DOCMP {
                    if (num_any_direction[X]==1) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                          }
                        }
                      }
                    }
                  }
                }
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
            if (is_real) {
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                }
              }
            } else { // not is_real
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
                     * f[ec][1][i];
                }
              }
            }
            DOCMP {
              for (int i=0;i<ntot;i++) {
                d_minus_p[ec][cmp][i] = f[dc][cmp][i];
                for (polarization *p = pol; p; p = p->next) {
                  d_minus_p[ec][cmp][i] -= p->P[ec][cmp][i];
                }
              }
              for (src *s = e_sources; s; s = s->next) {
                d_minus_p[ec][cmp][s->i] -= (cmp) ? imag(s->get_dipole_now()*s->A[ec])
                   : real(s->get_dipole_now()*s->A[ec]);
              }
            }
          }
          FOR_E_AND_D(ec,dc) if (f[ec][0]) {
            const int yee_idx = v.yee_index(ec);
            const int d_ec = component_direction(ec);
            const int s_ec = stride_any_direction[d_ec];
            const direction d_1 = (direction)((d_ec+1)%2);
            const component ec_1 = direction_component(ec,d_1);
            const int s_1 = stride_any_direction[d_1];
            if (ma->inveps[ec][d_1]) {
              DOCMP {
                if (num_any_direction[X]==1) {
                  for (int iY=0; iY<num_any_direction[Y]; iY++) {
                    const int i = yee_idx + iY;
                    f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                       + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                          + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                             + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                  }
                } else { // not num_any_direction[X]==1
                  if (num_any_direction[Y]==1) {
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      const int i = yee_idx + iX*stride_any_direction[X];
                      f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                         + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                            + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                               + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                    }
                  } else { // not num_any_direction[Y]==1
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iX*stride_any_direction[X] + iY;
                        f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                           + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                              + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                 + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                      }
                    }
                  }
                }
              }
            } else { // not ma->inveps[ec][d_1]
              DOCMP {
                if (num_any_direction[X]==1) {
                  for (int iY=0; iY<num_any_direction[Y]; iY++) {
                    const int i = yee_idx + iY;
                    f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                  }
                } else { // not num_any_direction[X]==1
                  if (num_any_direction[Y]==1) {
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      const int i = yee_idx + iX*stride_any_direction[X];
                      f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                    }
                  } else { // not num_any_direction[Y]==1
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iX*stride_any_direction[X] + iY;
                        f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                      }
                    }
                  }
                }
              }
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
            if (is_real) {
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                }
              }
            } else { // not is_real
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
                     * f[ec][1][i];
                }
              }
            }
            DOCMP {
              for (int i=0;i<ntot;i++) {
                d_minus_p[ec][cmp][i] = f[dc][cmp][i];
                for (polarization *p = pol; p; p = p->next) {
                  d_minus_p[ec][cmp][i] -= p->P[ec][cmp][i];
                }
              }
              for (src *s = e_sources; s; s = s->next) {
                d_minus_p[ec][cmp][s->i] -= (cmp) ? imag(s->get_dipole_now()*s->A[ec])
                   : real(s->get_dipole_now()*s->A[ec]);
              }
            }
          }
          FOR_E_AND_D(ec,dc) if (f[ec][0]) {
            const int yee_idx = v.yee_index(ec);
            const int d_ec = component_direction(ec);
            const int s_ec = stride_any_direction[d_ec];
            DOCMP {
              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                const int i = yee_idx + iZ;
                f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
              }
            }
          }
        } else { // not f[Ex][0]
          // Am in 2D TM
          FOR_E_AND_D(ec,dc) if (f[ec][0]) {
            const int yee_idx = v.yee_index(ec);
            const int d_ec = component_direction(ec);
            const int s_ec = stride_any_direction[d_ec];
            if (is_real) {
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                }
              }
            } else { // not is_real
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
                     * f[ec][1][i];
                }
              }
            }
            DOCMP {
              for (int i=0;i<ntot;i++) {
                d_minus_p[ec][cmp][i] = f[dc][cmp][i];
                for (polarization *p = pol; p; p = p->next) {
                  d_minus_p[ec][cmp][i] -= p->P[ec][cmp][i];
                }
              }
              for (src *s = e_sources; s; s = s->next) {
                d_minus_p[ec][cmp][s->i] -= (cmp) ? imag(s->get_dipole_now()*s->A[ec])
                   : real(s->get_dipole_now()*s->A[ec]);
              }
            }
          }
          FOR_E_AND_D(ec,dc) if (f[ec][0]) {
            const int yee_idx = v.yee_index(ec);
            const int d_ec = component_direction(ec);
            const int s_ec = stride_any_direction[d_ec];
            DOCMP {
              if (num_any_direction[X]==1) {
                for (int iY=0; iY<num_any_direction[Y]; iY++) {
                  const int i = yee_idx + iY;
                  f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                }
              } else { // not num_any_direction[X]==1
                if (num_any_direction[Y]==1) {
                  for (int iX=0; iX<num_any_direction[X]; iX++) {
                    const int i = yee_idx + iX*stride_any_direction[X];
                    f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                  }
                } else { // not num_any_direction[Y]==1
                  for (int iX=0; iX<num_any_direction[X]; iX++) {
                    for (int iY=0; iY<num_any_direction[Y]; iY++) {
                      const int i = yee_idx + iX*stride_any_direction[X] + iY;
                      f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  } else { // not e_sources
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
        if (is_real) {
          for (int i=0;i<ntot;i++) {
            for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
              np->energy[ec][i] = op->energy[ec][i];
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                 * f[ec][0][i];
            }
          }
        } else { // not is_real
          for (int i=0;i<ntot;i++) {
            for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
              np->energy[ec][i] = op->energy[ec][i];
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                 * f[ec][0][i];
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
                 * f[ec][1][i];
            }
          }
        }
        DOCMP {
          for (int i=0;i<ntot;i++) {
            d_minus_p[ec][cmp][i] = f[dc][cmp][i];
            for (polarization *p = pol; p; p = p->next) {
              d_minus_p[ec][cmp][i] -= p->P[ec][cmp][i];
            }
          }
        }
      }
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
        if (ma->inveps[ec][d_1]) {
          if (ma->inveps[ec][d_2]) {
            DOCMP {
              if (num_any_direction[Z]==1) {
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                     + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                        + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                           + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                              + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                 + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                       + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                          + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                             + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                   + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR++) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                         + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                            + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                               + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                  + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                     + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                    }
                  }
                }
              }
            }
          } else { // not ma->inveps[ec][d_2]
            DOCMP {
              if (num_any_direction[Z]==1) {
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                     + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                        + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                           + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                       + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                          + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                             + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR++) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                         + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                            + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                               + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                    }
                  }
                }
              }
            }
          }
        } else { // not ma->inveps[ec][d_1]
          if (ma->inveps[ec][d_2]) {
            DOCMP {
              if (num_any_direction[Z]==1) {
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                     + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                        + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                           + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                       + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                          + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                             + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR++) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                         + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                            + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                               + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                    }
                  }
                }
              }
            }
          } else { // not ma->inveps[ec][d_2]
            DOCMP {
              if (num_any_direction[Z]==1) {
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR++) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                    }
                  }
                }
              }
            }
          }
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
              if (is_real) {
                for (int i=0;i<ntot;i++) {
                  for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                    np->energy[ec][i] = op->energy[ec][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                       * f[ec][0][i];
                  }
                }
              } else { // not is_real
                for (int i=0;i<ntot;i++) {
                  for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                    np->energy[ec][i] = op->energy[ec][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                       * f[ec][0][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
                       * f[ec][1][i];
                  }
                }
              }
              DOCMP {
                for (int i=0;i<ntot;i++) {
                  d_minus_p[ec][cmp][i] = f[dc][cmp][i];
                  for (polarization *p = pol; p; p = p->next) {
                    d_minus_p[ec][cmp][i] -= p->P[ec][cmp][i];
                  }
                }
              }
            }
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
              if (ma->inveps[ec][d_1]) {
                if (ma->inveps[ec][d_2]) {
                  DOCMP {
                    if (num_any_direction[Z]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                   + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                      + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                         + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                  + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                     + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                        + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                           + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                 + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                    + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                       + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                          + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                             + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iY=0; iY<num_any_direction[Y]; iY++) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                   + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                      + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                         + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                            + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                               + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                } else { // not ma->inveps[ec][d_2]
                  DOCMP {
                    if (num_any_direction[Z]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                   + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                  + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                     + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                 + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                    + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                       + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iY=0; iY<num_any_direction[Y]; iY++) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                   + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                      + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                         + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              } else { // not ma->inveps[ec][d_1]
                if (ma->inveps[ec][d_2]) {
                  DOCMP {
                    if (num_any_direction[Z]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                   + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                  + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                     + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                 + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                    + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                       + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iY=0; iY<num_any_direction[Y]; iY++) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                   + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                      + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                         + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                } else { // not ma->inveps[ec][d_2]
                  DOCMP {
                    if (num_any_direction[Z]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                          f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iY=0; iY<num_any_direction[Y]; iY++) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
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
              if (is_real) {
                for (int i=0;i<ntot;i++) {
                  for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                    np->energy[ec][i] = op->energy[ec][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                       * f[ec][0][i];
                  }
                }
              } else { // not is_real
                for (int i=0;i<ntot;i++) {
                  for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                    np->energy[ec][i] = op->energy[ec][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                       * f[ec][0][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
                       * f[ec][1][i];
                  }
                }
              }
              DOCMP {
                for (int i=0;i<ntot;i++) {
                  d_minus_p[ec][cmp][i] = f[dc][cmp][i];
                  for (polarization *p = pol; p; p = p->next) {
                    d_minus_p[ec][cmp][i] -= p->P[ec][cmp][i];
                  }
                }
              }
            }
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
              if (ma->inveps[ec][d_1]) {
                if (ma->inveps[ec][d_2]) {
                  DOCMP {
                    if (num_any_direction[X]==1) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                           + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                              + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                 + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                    + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                       + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                   + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                      + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                         + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                  + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                     + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                        + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                           + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                          }
                        }
                      }
                    }
                  }
                } else { // not ma->inveps[ec][d_2]
                  DOCMP {
                    if (num_any_direction[X]==1) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                           + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                              + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                 + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                   + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                  + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                     + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                          }
                        }
                      }
                    }
                  }
                }
              } else { // not ma->inveps[ec][d_1]
                if (ma->inveps[ec][d_2]) {
                  DOCMP {
                    if (num_any_direction[X]==1) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                           + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                              + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                 + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                   + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                  + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                     + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                          }
                        }
                      }
                    }
                  }
                } else { // not ma->inveps[ec][d_2]
                  DOCMP {
                    if (num_any_direction[X]==1) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                          }
                        }
                      }
                    }
                  }
                }
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
            if (is_real) {
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                }
              }
            } else { // not is_real
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
                     * f[ec][1][i];
                }
              }
            }
            DOCMP {
              for (int i=0;i<ntot;i++) {
                d_minus_p[ec][cmp][i] = f[dc][cmp][i];
                for (polarization *p = pol; p; p = p->next) {
                  d_minus_p[ec][cmp][i] -= p->P[ec][cmp][i];
                }
              }
            }
          }
          FOR_E_AND_D(ec,dc) if (f[ec][0]) {
            const int yee_idx = v.yee_index(ec);
            const int d_ec = component_direction(ec);
            const int s_ec = stride_any_direction[d_ec];
            const direction d_1 = (direction)((d_ec+1)%2);
            const component ec_1 = direction_component(ec,d_1);
            const int s_1 = stride_any_direction[d_1];
            if (ma->inveps[ec][d_1]) {
              DOCMP {
                if (num_any_direction[X]==1) {
                  for (int iY=0; iY<num_any_direction[Y]; iY++) {
                    const int i = yee_idx + iY;
                    f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                       + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                          + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                             + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                  }
                } else { // not num_any_direction[X]==1
                  if (num_any_direction[Y]==1) {
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      const int i = yee_idx + iX*stride_any_direction[X];
                      f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                         + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                            + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                               + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                    }
                  } else { // not num_any_direction[Y]==1
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iX*stride_any_direction[X] + iY;
                        f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                           + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                              + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                 + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                      }
                    }
                  }
                }
              }
            } else { // not ma->inveps[ec][d_1]
              DOCMP {
                if (num_any_direction[X]==1) {
                  for (int iY=0; iY<num_any_direction[Y]; iY++) {
                    const int i = yee_idx + iY;
                    f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                  }
                } else { // not num_any_direction[X]==1
                  if (num_any_direction[Y]==1) {
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      const int i = yee_idx + iX*stride_any_direction[X];
                      f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                    }
                  } else { // not num_any_direction[Y]==1
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iX*stride_any_direction[X] + iY;
                        f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                      }
                    }
                  }
                }
              }
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
            if (is_real) {
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                }
              }
            } else { // not is_real
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
                     * f[ec][1][i];
                }
              }
            }
            DOCMP {
              for (int i=0;i<ntot;i++) {
                d_minus_p[ec][cmp][i] = f[dc][cmp][i];
                for (polarization *p = pol; p; p = p->next) {
                  d_minus_p[ec][cmp][i] -= p->P[ec][cmp][i];
                }
              }
            }
          }
          FOR_E_AND_D(ec,dc) if (f[ec][0]) {
            const int yee_idx = v.yee_index(ec);
            const int d_ec = component_direction(ec);
            const int s_ec = stride_any_direction[d_ec];
            DOCMP {
              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                const int i = yee_idx + iZ;
                f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
              }
            }
          }
        } else { // not f[Ex][0]
          // Am in 2D TM
          FOR_E_AND_D(ec,dc) if (f[ec][0]) {
            const int yee_idx = v.yee_index(ec);
            const int d_ec = component_direction(ec);
            const int s_ec = stride_any_direction[d_ec];
            if (is_real) {
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                }
              }
            } else { // not is_real
              for (int i=0;i<ntot;i++) {
                for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                  np->energy[ec][i] = op->energy[ec][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])
                     * f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])
                     * f[ec][1][i];
                }
              }
            }
            DOCMP {
              for (int i=0;i<ntot;i++) {
                d_minus_p[ec][cmp][i] = f[dc][cmp][i];
                for (polarization *p = pol; p; p = p->next) {
                  d_minus_p[ec][cmp][i] -= p->P[ec][cmp][i];
                }
              }
            }
          }
          FOR_E_AND_D(ec,dc) if (f[ec][0]) {
            const int yee_idx = v.yee_index(ec);
            const int d_ec = component_direction(ec);
            const int s_ec = stride_any_direction[d_ec];
            DOCMP {
              if (num_any_direction[X]==1) {
                for (int iY=0; iY<num_any_direction[Y]; iY++) {
                  const int i = yee_idx + iY;
                  f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                }
              } else { // not num_any_direction[X]==1
                if (num_any_direction[Y]==1) {
                  for (int iX=0; iX<num_any_direction[X]; iX++) {
                    const int i = yee_idx + iX*stride_any_direction[X];
                    f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                  }
                } else { // not num_any_direction[Y]==1
                  for (int iX=0; iX<num_any_direction[X]; iX++) {
                    for (int iY=0; iY<num_any_direction[Y]; iY++) {
                      const int i = yee_idx + iX*stride_any_direction[X] + iY;
                      f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  DOCMP {
    FOR_ELECTRIC_COMPONENTS(ec) {
      delete[] d_minus_p[ec][cmp];
    }
  }
} else { // not pol
  if (e_sources) {
    double *d_minus_p[5][2];
    DOCMP {
      FOR_ELECTRIC_COMPONENTS(ec) {
        if (f[ec]) {
          d_minus_p[ec][cmp] = new double[v.ntot()];
        } else { // not f[ec]
          d_minus_p[ec][cmp] = 0;
        }
      }
    }
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
        DOCMP {
          for (int i=0;i<ntot;i++) {
            d_minus_p[ec][cmp][i] = f[dc][cmp][i];
          }
          for (src *s = e_sources; s; s = s->next) {
            d_minus_p[ec][cmp][s->i] -= (cmp) ? imag(s->get_dipole_now()*s->A[ec])
               : real(s->get_dipole_now()*s->A[ec]);
          }
        }
      }
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
        if (ma->inveps[ec][d_1]) {
          if (ma->inveps[ec][d_2]) {
            DOCMP {
              if (num_any_direction[Z]==1) {
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                     + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                        + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                           + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                              + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                 + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                       + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                          + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                             + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                   + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR++) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                         + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                            + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                               + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                  + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                     + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                    }
                  }
                }
              }
            }
          } else { // not ma->inveps[ec][d_2]
            DOCMP {
              if (num_any_direction[Z]==1) {
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                     + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                        + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                           + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                       + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                          + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                             + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR++) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                         + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                            + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                               + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                    }
                  }
                }
              }
            }
          }
        } else { // not ma->inveps[ec][d_1]
          if (ma->inveps[ec][d_2]) {
            DOCMP {
              if (num_any_direction[Z]==1) {
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                     + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                        + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                           + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                       + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                          + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                             + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR++) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                         + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                            + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                               + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                    }
                  }
                }
              }
            }
          } else { // not ma->inveps[ec][d_2]
            DOCMP {
              if (num_any_direction[Z]==1) {
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR++) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                    }
                  }
                }
              }
            }
          }
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
              DOCMP {
                for (int i=0;i<ntot;i++) {
                  d_minus_p[ec][cmp][i] = f[dc][cmp][i];
                }
                for (src *s = e_sources; s; s = s->next) {
                  d_minus_p[ec][cmp][s->i] -= (cmp) ? imag(s->get_dipole_now()*s->A[ec])
                     : real(s->get_dipole_now()*s->A[ec]);
                }
              }
            }
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
              if (ma->inveps[ec][d_1]) {
                if (ma->inveps[ec][d_2]) {
                  DOCMP {
                    if (num_any_direction[Z]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                   + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                      + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                         + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                  + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                     + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                        + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                           + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                 + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                    + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                       + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                          + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                             + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iY=0; iY<num_any_direction[Y]; iY++) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                   + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                      + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                         + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                            + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                               + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                } else { // not ma->inveps[ec][d_2]
                  DOCMP {
                    if (num_any_direction[Z]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                   + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                  + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                     + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                 + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                    + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                       + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iY=0; iY<num_any_direction[Y]; iY++) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                   + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                      + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                         + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              } else { // not ma->inveps[ec][d_1]
                if (ma->inveps[ec][d_2]) {
                  DOCMP {
                    if (num_any_direction[Z]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                   + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                  + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                     + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                 + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                    + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                       + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iY=0; iY<num_any_direction[Y]; iY++) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                                   + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                      + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                         + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                } else { // not ma->inveps[ec][d_2]
                  DOCMP {
                    if (num_any_direction[Z]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                          f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iY=0; iY<num_any_direction[Y]; iY++) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
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
              DOCMP {
                for (int i=0;i<ntot;i++) {
                  d_minus_p[ec][cmp][i] = f[dc][cmp][i];
                }
                for (src *s = e_sources; s; s = s->next) {
                  d_minus_p[ec][cmp][s->i] -= (cmp) ? imag(s->get_dipole_now()*s->A[ec])
                     : real(s->get_dipole_now()*s->A[ec]);
                }
              }
            }
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
              if (ma->inveps[ec][d_1]) {
                if (ma->inveps[ec][d_2]) {
                  DOCMP {
                    if (num_any_direction[X]==1) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                           + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                              + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                 + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                    + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                       + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                   + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                      + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                         + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + (0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                  + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                     + d_minus_p[ec_1][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                        + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                           + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                          }
                        }
                      }
                    }
                  }
                } else { // not ma->inveps[ec][d_2]
                  DOCMP {
                    if (num_any_direction[X]==1) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                           + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                              + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                 + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                   + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                                  + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                     + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                          }
                        }
                      }
                    }
                  }
                }
              } else { // not ma->inveps[ec][d_1]
                if (ma->inveps[ec][d_2]) {
                  DOCMP {
                    if (num_any_direction[X]==1) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                           + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                              + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                 + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                             + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                   + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                               + 0.25*((ma->inveps[ec][d_2][i])*(d_minus_p[ec_2][cmp][i]
                                  + (d_minus_p[ec_2][cmp][i-s_2]) + (d_minus_p[ec_2][cmp][i+s_ec])
                                     + d_minus_p[ec_2][cmp][i-s_2+s_ec]));
                          }
                        }
                      }
                    }
                  }
                } else { // not ma->inveps[ec][d_2]
                  DOCMP {
                    if (num_any_direction[X]==1) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                          }
                        }
                      }
                    }
                  }
                }
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
            DOCMP {
              for (int i=0;i<ntot;i++) {
                d_minus_p[ec][cmp][i] = f[dc][cmp][i];
              }
              for (src *s = e_sources; s; s = s->next) {
                d_minus_p[ec][cmp][s->i] -= (cmp) ? imag(s->get_dipole_now()*s->A[ec])
                   : real(s->get_dipole_now()*s->A[ec]);
              }
            }
          }
          FOR_E_AND_D(ec,dc) if (f[ec][0]) {
            const int yee_idx = v.yee_index(ec);
            const int d_ec = component_direction(ec);
            const int s_ec = stride_any_direction[d_ec];
            const direction d_1 = (direction)((d_ec+1)%2);
            const component ec_1 = direction_component(ec,d_1);
            const int s_1 = stride_any_direction[d_1];
            if (ma->inveps[ec][d_1]) {
              DOCMP {
                if (num_any_direction[X]==1) {
                  for (int iY=0; iY<num_any_direction[Y]; iY++) {
                    const int i = yee_idx + iY;
                    f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                       + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                          + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                             + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                  }
                } else { // not num_any_direction[X]==1
                  if (num_any_direction[Y]==1) {
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      const int i = yee_idx + iX*stride_any_direction[X];
                      f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                         + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                            + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                               + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                    }
                  } else { // not num_any_direction[Y]==1
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iX*stride_any_direction[X] + iY;
                        f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i])
                           + 0.25*((ma->inveps[ec][d_1][i])*(d_minus_p[ec_1][cmp][i]
                              + (d_minus_p[ec_1][cmp][i-s_1]) + (d_minus_p[ec_1][cmp][i+s_ec])
                                 + d_minus_p[ec_1][cmp][i-s_1+s_ec]));
                      }
                    }
                  }
                }
              }
            } else { // not ma->inveps[ec][d_1]
              DOCMP {
                if (num_any_direction[X]==1) {
                  for (int iY=0; iY<num_any_direction[Y]; iY++) {
                    const int i = yee_idx + iY;
                    f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                  }
                } else { // not num_any_direction[X]==1
                  if (num_any_direction[Y]==1) {
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      const int i = yee_idx + iX*stride_any_direction[X];
                      f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                    }
                  } else { // not num_any_direction[Y]==1
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iX*stride_any_direction[X] + iY;
                        f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                      }
                    }
                  }
                }
              }
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
            DOCMP {
              for (int i=0;i<ntot;i++) {
                d_minus_p[ec][cmp][i] = f[dc][cmp][i];
              }
              for (src *s = e_sources; s; s = s->next) {
                d_minus_p[ec][cmp][s->i] -= (cmp) ? imag(s->get_dipole_now()*s->A[ec])
                   : real(s->get_dipole_now()*s->A[ec]);
              }
            }
          }
          FOR_E_AND_D(ec,dc) if (f[ec][0]) {
            const int yee_idx = v.yee_index(ec);
            const int d_ec = component_direction(ec);
            const int s_ec = stride_any_direction[d_ec];
            DOCMP {
              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                const int i = yee_idx + iZ;
                f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
              }
            }
          }
        } else { // not f[Ex][0]
          // Am in 2D TM
          FOR_E_AND_D(ec,dc) if (f[ec][0]) {
            const int yee_idx = v.yee_index(ec);
            const int d_ec = component_direction(ec);
            const int s_ec = stride_any_direction[d_ec];
            DOCMP {
              for (int i=0;i<ntot;i++) {
                d_minus_p[ec][cmp][i] = f[dc][cmp][i];
              }
              for (src *s = e_sources; s; s = s->next) {
                d_minus_p[ec][cmp][s->i] -= (cmp) ? imag(s->get_dipole_now()*s->A[ec])
                   : real(s->get_dipole_now()*s->A[ec]);
              }
            }
          }
          FOR_E_AND_D(ec,dc) if (f[ec][0]) {
            const int yee_idx = v.yee_index(ec);
            const int d_ec = component_direction(ec);
            const int s_ec = stride_any_direction[d_ec];
            DOCMP {
              if (num_any_direction[X]==1) {
                for (int iY=0; iY<num_any_direction[Y]; iY++) {
                  const int i = yee_idx + iY;
                  f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                }
              } else { // not num_any_direction[X]==1
                if (num_any_direction[Y]==1) {
                  for (int iX=0; iX<num_any_direction[X]; iX++) {
                    const int i = yee_idx + iX*stride_any_direction[X];
                    f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                  }
                } else { // not num_any_direction[Y]==1
                  for (int iX=0; iX<num_any_direction[X]; iX++) {
                    for (int iY=0; iY<num_any_direction[Y]; iY++) {
                      const int i = yee_idx + iX*stride_any_direction[X] + iY;
                      f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*d_minus_p[ec][cmp][i];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    DOCMP {
      FOR_ELECTRIC_COMPONENTS(ec) {
        delete[] d_minus_p[ec][cmp];
      }
    }
  } else { // not e_sources
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
      }
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
        if (ma->inveps[ec][d_1]) {
          if (ma->inveps[ec][d_2]) {
            DOCMP {
              if (num_any_direction[Z]==1) {
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                     + (0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                        + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                           + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                              + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                 + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                       + (0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                          + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                             + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                                + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                   + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR++) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                         + (0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                            + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                               + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                                  + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                     + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                    }
                  }
                }
              }
            }
          } else { // not ma->inveps[ec][d_2]
            DOCMP {
              if (num_any_direction[Z]==1) {
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                     + 0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                        + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                           + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]));
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                       + 0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                          + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                             + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]));
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR++) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                         + 0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                            + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                               + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]));
                    }
                  }
                }
              }
            }
          }
        } else { // not ma->inveps[ec][d_1]
          if (ma->inveps[ec][d_2]) {
            DOCMP {
              if (num_any_direction[Z]==1) {
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                     + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                        + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                           + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                       + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                          + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                             + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR++) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                         + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                            + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                               + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                    }
                  }
                }
              }
            }
          } else { // not ma->inveps[ec][d_2]
            DOCMP {
              if (num_any_direction[Z]==1) {
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR++) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
                    }
                  }
                }
              }
            }
          }
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
            }
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
              if (ma->inveps[ec][d_1]) {
                if (ma->inveps[ec][d_2]) {
                  DOCMP {
                    if (num_any_direction[Z]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                             + (0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                                + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                                   + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                                      + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                         + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                               + (0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                                  + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                                     + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                                        + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                           + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                                 + (0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                                    + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                                       + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                                          + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                             + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iY=0; iY<num_any_direction[Y]; iY++) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                                   + (0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                                      + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                                         + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                                            + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                               + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                } else { // not ma->inveps[ec][d_2]
                  DOCMP {
                    if (num_any_direction[Z]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                             + 0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                                + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                                   + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]));
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                               + 0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                                  + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                                     + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]));
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                                 + 0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                                    + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                                       + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]));
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iY=0; iY<num_any_direction[Y]; iY++) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                                   + 0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                                      + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                                         + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]));
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              } else { // not ma->inveps[ec][d_1]
                if (ma->inveps[ec][d_2]) {
                  DOCMP {
                    if (num_any_direction[Z]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                             + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                                + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                   + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                               + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                                  + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                     + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                                 + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                                    + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                       + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iY=0; iY<num_any_direction[Y]; iY++) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                                   + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                                      + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                         + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                } else { // not ma->inveps[ec][d_2]
                  DOCMP {
                    if (num_any_direction[Z]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                          f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX++) {
                            for (int iY=0; iY<num_any_direction[Y]; iY++) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
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
            }
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
              if (ma->inveps[ec][d_1]) {
                if (ma->inveps[ec][d_2]) {
                  DOCMP {
                    if (num_any_direction[X]==1) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                           + (0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                              + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                                 + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                                    + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                       + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                             + (0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                                + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                                   + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                                      + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                         + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                               + (0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                                  + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                                     + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]))) + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                                        + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                           + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                          }
                        }
                      }
                    }
                  }
                } else { // not ma->inveps[ec][d_2]
                  DOCMP {
                    if (num_any_direction[X]==1) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                           + 0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                              + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                                 + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]));
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                             + 0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                                + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                                   + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]));
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                               + 0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                                  + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                                     + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]));
                          }
                        }
                      }
                    }
                  }
                }
              } else { // not ma->inveps[ec][d_1]
                if (ma->inveps[ec][d_2]) {
                  DOCMP {
                    if (num_any_direction[X]==1) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                           + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                              + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                 + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                             + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                                + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                   + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                               + 0.25*((ma->inveps[ec][d_2][i])*((f[(component)(ec_2+10)][cmp][i])
                                  + (f[(component)(ec_2+10)][cmp][i-s_2]) + (f[(component)(ec_2+10)][cmp][i+s_ec])
                                     + f[(component)(ec_2+10)][cmp][i-s_2+s_ec]));
                          }
                        }
                      }
                    }
                  }
                } else { // not ma->inveps[ec][d_2]
                  DOCMP {
                    if (num_any_direction[X]==1) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
                          }
                        }
                      }
                    }
                  }
                }
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
          }
          FOR_E_AND_D(ec,dc) if (f[ec][0]) {
            const int yee_idx = v.yee_index(ec);
            const int d_ec = component_direction(ec);
            const int s_ec = stride_any_direction[d_ec];
            const direction d_1 = (direction)((d_ec+1)%2);
            const component ec_1 = direction_component(ec,d_1);
            const int s_1 = stride_any_direction[d_1];
            if (ma->inveps[ec][d_1]) {
              DOCMP {
                if (num_any_direction[X]==1) {
                  for (int iY=0; iY<num_any_direction[Y]; iY++) {
                    const int i = yee_idx + iY;
                    f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                       + 0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                          + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                             + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]));
                  }
                } else { // not num_any_direction[X]==1
                  if (num_any_direction[Y]==1) {
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      const int i = yee_idx + iX*stride_any_direction[X];
                      f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                         + 0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                            + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                               + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]));
                    }
                  } else { // not num_any_direction[Y]==1
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iX*stride_any_direction[X] + iY;
                        f[ec][cmp][i] = ((ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]))
                           + 0.25*((ma->inveps[ec][d_1][i])*((f[(component)(ec_1+10)][cmp][i])
                              + (f[(component)(ec_1+10)][cmp][i-s_1]) + (f[(component)(ec_1+10)][cmp][i+s_ec])
                                 + f[(component)(ec_1+10)][cmp][i-s_1+s_ec]));
                      }
                    }
                  }
                }
              }
            } else { // not ma->inveps[ec][d_1]
              DOCMP {
                if (num_any_direction[X]==1) {
                  for (int iY=0; iY<num_any_direction[Y]; iY++) {
                    const int i = yee_idx + iY;
                    f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
                  }
                } else { // not num_any_direction[X]==1
                  if (num_any_direction[Y]==1) {
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      const int i = yee_idx + iX*stride_any_direction[X];
                      f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
                    }
                  } else { // not num_any_direction[Y]==1
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iX*stride_any_direction[X] + iY;
                        f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
                      }
                    }
                  }
                }
              }
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
          }
          FOR_E_AND_D(ec,dc) if (f[ec][0]) {
            const int yee_idx = v.yee_index(ec);
            const int d_ec = component_direction(ec);
            const int s_ec = stride_any_direction[d_ec];
            DOCMP {
              for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                const int i = yee_idx + iZ;
                f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
              }
            }
          }
        } else { // not f[Ex][0]
          // Am in 2D TM
          FOR_E_AND_D(ec,dc) if (f[ec][0]) {
            const int yee_idx = v.yee_index(ec);
            const int d_ec = component_direction(ec);
            const int s_ec = stride_any_direction[d_ec];
          }
          FOR_E_AND_D(ec,dc) if (f[ec][0]) {
            const int yee_idx = v.yee_index(ec);
            const int d_ec = component_direction(ec);
            const int s_ec = stride_any_direction[d_ec];
            DOCMP {
              if (num_any_direction[X]==1) {
                for (int iY=0; iY<num_any_direction[Y]; iY++) {
                  const int i = yee_idx + iY;
                  f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
                }
              } else { // not num_any_direction[X]==1
                if (num_any_direction[Y]==1) {
                  for (int iX=0; iX<num_any_direction[X]; iX++) {
                    const int i = yee_idx + iX*stride_any_direction[X];
                    f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
                  }
                } else { // not num_any_direction[Y]==1
                  for (int iX=0; iX<num_any_direction[X]; iX++) {
                    for (int iY=0; iY<num_any_direction[Y]; iY++) {
                      const int i = yee_idx + iX*stride_any_direction[X] + iY;
                      f[ec][cmp][i] = (ma->inveps[ec][d_ec][i])*(f[(component)(ec+10)][cmp][i]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
