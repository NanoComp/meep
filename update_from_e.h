if (f[Er][0]) {
  FOR_ELECTRIC_COMPONENTS(ec) if (f[ec][0]) {
    const int yee_idx = v.yee_index(ec);
    if (is_real) {
      if (pol) {
        for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
          const double fac = np->saturation_factor;
          const double g = op->pb->gamma;
          const double om = op->pb->omeganot;
          const double invomsqr = 1.0/(om*om);
          const double funinv = 1.0/(1+0.5*g);
          if (fac) {
            if (fac > 0) {
              for (int i=0;i<ntot;i++) {
                np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
              }
              if (num_any_direction[Z]==1) {
                if (num_any_direction[R]==1) {
                  const int i = yee_idx;
                  const double energy_here = (np->energy[Ez][i]) + (np->energy[Ep][i])
                     + np->energy[Er][i];
                  np->s[ec][i] = max(-energy_here*fac, 0.0);
                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                } else {
                  for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                    const int i = yee_idx + iR*stride_any_direction[R];
                    const double energy_here = (np->energy[Ez][i]) + (np->energy[Ep][i])
                       + np->energy[Er][i];
                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                }
              } else {
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                    const int i = yee_idx + iZ;
                    const double energy_here = (np->energy[Ez][i]) + (np->energy[Ep][i])
                       + np->energy[Er][i];
                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                } else {
                  for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      const double energy_here = (np->energy[Ez][i]) + (np->energy[Ep][i])
                         + np->energy[Er][i];
                      np->s[ec][i] = max(-energy_here*fac, 0.0);
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    }
                  }
                }
              }
            } else {
              for (int i=0;i<ntot;i++) {
                np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
              }
              if (num_any_direction[Z]==1) {
                if (num_any_direction[R]==1) {
                  const int i = yee_idx;
                  const double energy_here = (np->energy[Ez][i]) + (np->energy[Ep][i])
                     + np->energy[Er][i];
                  np->s[ec][i] = energy_here*fac;
                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                } else {
                  for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                    const int i = yee_idx + iR*stride_any_direction[R];
                    const double energy_here = (np->energy[Ez][i]) + (np->energy[Ep][i])
                       + np->energy[Er][i];
                    np->s[ec][i] = energy_here*fac;
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                }
              } else {
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                    const int i = yee_idx + iZ;
                    const double energy_here = (np->energy[Ez][i]) + (np->energy[Ep][i])
                       + np->energy[Er][i];
                    np->s[ec][i] = energy_here*fac;
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                } else {
                  for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      const double energy_here = (np->energy[Ez][i]) + (np->energy[Ep][i])
                         + np->energy[Er][i];
                      np->s[ec][i] = energy_here*fac;
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    }
                  }
                }
              }
            }
          } else {
            for (int i=0;i<ntot;i++) {
              np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                 + np->s[ec][i]*f[ec][0][i]);
            }
          }
        }
      }
    } else {
      if (pol) {
        for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
          const double fac = np->saturation_factor;
          const double g = op->pb->gamma;
          const double om = op->pb->omeganot;
          const double invomsqr = 1.0/(om*om);
          const double funinv = 1.0/(1+0.5*g);
          if (fac) {
            if (fac > 0) {
              for (int i=0;i<ntot;i++) {
                np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
              }
              if (num_any_direction[Z]==1) {
                if (num_any_direction[R]==1) {
                  const int i = yee_idx;
                  const double energy_here = (np->energy[Ez][i]) + (np->energy[Ep][i])
                     + np->energy[Er][i];
                  np->s[ec][i] = max(-energy_here*fac, 0.0);
                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                     + np->s[ec][i]*f[ec][1][i]);
                } else {
                  for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                    const int i = yee_idx + iR*stride_any_direction[R];
                    const double energy_here = (np->energy[Ez][i]) + (np->energy[Ep][i])
                       + np->energy[Er][i];
                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  }
                }
              } else {
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                    const int i = yee_idx + iZ;
                    const double energy_here = (np->energy[Ez][i]) + (np->energy[Ep][i])
                       + np->energy[Er][i];
                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  }
                } else {
                  for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      const double energy_here = (np->energy[Ez][i]) + (np->energy[Ep][i])
                         + np->energy[Er][i];
                      np->s[ec][i] = max(-energy_here*fac, 0.0);
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                      op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                         + np->s[ec][i]*f[ec][1][i]);
                    }
                  }
                }
              }
            } else {
              for (int i=0;i<ntot;i++) {
                np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
              }
              if (num_any_direction[Z]==1) {
                if (num_any_direction[R]==1) {
                  const int i = yee_idx;
                  const double energy_here = (np->energy[Ez][i]) + (np->energy[Ep][i])
                     + np->energy[Er][i];
                  np->s[ec][i] = energy_here*fac;
                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                     + np->s[ec][i]*f[ec][1][i]);
                } else {
                  for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                    const int i = yee_idx + iR*stride_any_direction[R];
                    const double energy_here = (np->energy[Ez][i]) + (np->energy[Ep][i])
                       + np->energy[Er][i];
                    np->s[ec][i] = energy_here*fac;
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  }
                }
              } else {
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                    const int i = yee_idx + iZ;
                    const double energy_here = (np->energy[Ez][i]) + (np->energy[Ep][i])
                       + np->energy[Er][i];
                    np->s[ec][i] = energy_here*fac;
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  }
                } else {
                  for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      const double energy_here = (np->energy[Ez][i]) + (np->energy[Ep][i])
                         + np->energy[Er][i];
                      np->s[ec][i] = energy_here*fac;
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                      op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                         + np->s[ec][i]*f[ec][1][i]);
                    }
                  }
                }
              }
            }
          } else {
            for (int i=0;i<ntot;i++) {
              np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
              np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                 + np->s[ec][i]*f[ec][0][i]);
              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                 + np->s[ec][i]*f[ec][1][i]);
            }
          }
        }
      }
    }
  }
  // The polarizations got switched...
  polarization *temp = olpol;
  olpol = pol;
  pol = temp;
} else {
  if (stride_any_direction[Z]) {
    if (f[Ey][0]) {
      if (f[Ez][0]) {
        FOR_ELECTRIC_COMPONENTS(ec) if (f[ec][0]) {
          const int yee_idx = v.yee_index(ec);
          if (is_real) {
            if (pol) {
              for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                const double fac = np->saturation_factor;
                const double g = op->pb->gamma;
                const double om = op->pb->omeganot;
                const double invomsqr = 1.0/(om*om);
                const double funinv = 1.0/(1+0.5*g);
                if (fac) {
                  if (fac > 0) {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    }
                    if (num_any_direction[Z]==1) {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          const int i = yee_idx;
                          const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                             + np->energy[Ez][i];
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      } else {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            const int i = yee_idx + iX*stride_any_direction[X];
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    } else {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                            const int i = yee_idx + iZ;
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      } else {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                   + np->energy[Ez][i];
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            }
                          }
                        }
                      }
                    }
                  } else {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    }
                    if (num_any_direction[Z]==1) {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          const int i = yee_idx;
                          const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                             + np->energy[Ez][i];
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      } else {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            const int i = yee_idx + iX*stride_any_direction[X];
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    } else {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                            const int i = yee_idx + iZ;
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      } else {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                   + np->energy[Ez][i];
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                } else {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                }
              }
            }
          } else {
            if (pol) {
              for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                const double fac = np->saturation_factor;
                const double g = op->pb->gamma;
                const double om = op->pb->omeganot;
                const double invomsqr = 1.0/(om*om);
                const double funinv = 1.0/(1+0.5*g);
                if (fac) {
                  if (fac > 0) {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                      np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    }
                    if (num_any_direction[Z]==1) {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          const int i = yee_idx;
                          const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                             + np->energy[Ez][i];
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      } else {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            const int i = yee_idx + iX*stride_any_direction[X];
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    } else {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                            const int i = yee_idx + iZ;
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      } else {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                   + np->energy[Ez][i];
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                                op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                   + np->s[ec][i]*f[ec][1][i]);
                              }
                            }
                          }
                        }
                      }
                    }
                  } else {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                      np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    }
                    if (num_any_direction[Z]==1) {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          const int i = yee_idx;
                          const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                             + np->energy[Ez][i];
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      } else {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            const int i = yee_idx + iX*stride_any_direction[X];
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    } else {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                            const int i = yee_idx + iZ;
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      } else {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                   + np->energy[Ez][i];
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                                op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                   + np->s[ec][i]*f[ec][1][i]);
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                } else {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  }
                }
              }
            }
          }
        }
        // The polarizations got switched...
        polarization *temp = olpol;
        olpol = pol;
        pol = temp;
      } else {
        FOR_ELECTRIC_COMPONENTS(ec) if (f[ec][0]) {
          const int yee_idx = v.yee_index(ec);
          if (is_real) {
            if (pol) {
              for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                const double fac = np->saturation_factor;
                const double g = op->pb->gamma;
                const double om = op->pb->omeganot;
                const double invomsqr = 1.0/(om*om);
                const double funinv = 1.0/(1+0.5*g);
                if (fac) {
                  if (fac > 0) {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    }
                    if (num_any_direction[Z]==1) {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          const int i = yee_idx;
                          const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      } else {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            const int i = yee_idx + iX*stride_any_direction[X];
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    } else {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                            const int i = yee_idx + iZ;
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      } else {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            }
                          }
                        }
                      }
                    }
                  } else {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    }
                    if (num_any_direction[Z]==1) {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          const int i = yee_idx;
                          const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      } else {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            const int i = yee_idx + iX*stride_any_direction[X];
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    } else {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                            const int i = yee_idx + iZ;
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      } else {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                } else {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                }
              }
            }
          } else {
            if (pol) {
              for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                const double fac = np->saturation_factor;
                const double g = op->pb->gamma;
                const double om = op->pb->omeganot;
                const double invomsqr = 1.0/(om*om);
                const double funinv = 1.0/(1+0.5*g);
                if (fac) {
                  if (fac > 0) {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                      np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    }
                    if (num_any_direction[Z]==1) {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          const int i = yee_idx;
                          const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      } else {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            const int i = yee_idx + iX*stride_any_direction[X];
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    } else {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                            const int i = yee_idx + iZ;
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      } else {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                                op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                   + np->s[ec][i]*f[ec][1][i]);
                              }
                            }
                          }
                        }
                      }
                    }
                  } else {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                      np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    }
                    if (num_any_direction[Z]==1) {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          const int i = yee_idx;
                          const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      } else {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            const int i = yee_idx + iX*stride_any_direction[X];
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    } else {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                            const int i = yee_idx + iZ;
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      } else {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                                op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                   + np->s[ec][i]*f[ec][1][i]);
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                } else {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  }
                }
              }
            }
          }
        }
        // The polarizations got switched...
        polarization *temp = olpol;
        olpol = pol;
        pol = temp;
      }
    } else {
      FOR_ELECTRIC_COMPONENTS(ec) if (f[ec][0]) {
        const int yee_idx = v.yee_index(ec);
        if (is_real) {
          if (pol) {
            for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
              const double fac = np->saturation_factor;
              const double g = op->pb->gamma;
              const double om = op->pb->omeganot;
              const double invomsqr = 1.0/(om*om);
              const double funinv = 1.0/(1+0.5*g);
              if (fac) {
                if (fac > 0) {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  }
                  if (num_any_direction[Z]==1) {
                    const int i = yee_idx;
                    const double energy_here = np->energy[Ex][i];
                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  } else {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                      const int i = yee_idx + iZ;
                      const double energy_here = np->energy[Ex][i];
                      np->s[ec][i] = max(-energy_here*fac, 0.0);
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    }
                  }
                } else {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  }
                  if (num_any_direction[Z]==1) {
                    const int i = yee_idx;
                    const double energy_here = np->energy[Ex][i];
                    np->s[ec][i] = energy_here*fac;
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  } else {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                      const int i = yee_idx + iZ;
                      const double energy_here = np->energy[Ex][i];
                      np->s[ec][i] = energy_here*fac;
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    }
                  }
                }
              } else {
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                }
              }
            }
          }
        } else {
          if (pol) {
            for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
              const double fac = np->saturation_factor;
              const double g = op->pb->gamma;
              const double om = op->pb->omeganot;
              const double invomsqr = 1.0/(om*om);
              const double funinv = 1.0/(1+0.5*g);
              if (fac) {
                if (fac > 0) {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                  }
                  if (num_any_direction[Z]==1) {
                    const int i = yee_idx;
                    const double energy_here = np->energy[Ex][i];
                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  } else {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                      const int i = yee_idx + iZ;
                      const double energy_here = np->energy[Ex][i];
                      np->s[ec][i] = max(-energy_here*fac, 0.0);
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                      op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                         + np->s[ec][i]*f[ec][1][i]);
                    }
                  }
                } else {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                  }
                  if (num_any_direction[Z]==1) {
                    const int i = yee_idx;
                    const double energy_here = np->energy[Ex][i];
                    np->s[ec][i] = energy_here*fac;
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  } else {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                      const int i = yee_idx + iZ;
                      const double energy_here = np->energy[Ex][i];
                      np->s[ec][i] = energy_here*fac;
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                      op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                         + np->s[ec][i]*f[ec][1][i]);
                    }
                  }
                }
              } else {
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                     + np->s[ec][i]*f[ec][1][i]);
                }
              }
            }
          }
        }
      }
      // The polarizations got switched...
      polarization *temp = olpol;
      olpol = pol;
      pol = temp;
    }
  } else {
    if (f[Ey][0]) {
      if (f[Ez][0]) {
        FOR_ELECTRIC_COMPONENTS(ec) if (f[ec][0]) {
          const int yee_idx = v.yee_index(ec);
          if (is_real) {
            if (pol) {
              for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                const double fac = np->saturation_factor;
                const double g = op->pb->gamma;
                const double om = op->pb->omeganot;
                const double invomsqr = 1.0/(om*om);
                const double funinv = 1.0/(1+0.5*g);
                if (fac) {
                  if (fac > 0) {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    }
                    if (num_any_direction[X]==1) {
                      if (num_any_direction[Y]==1) {
                        const int i = yee_idx;
                        const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                           + np->energy[Ez][i];
                        np->s[ec][i] = max(-energy_here*fac, 0.0);
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      } else {
                        if (stride_any_direction[Y]==1) {
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iY;
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      }
                    } else {
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                             + np->energy[Ez][i];
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      } else {
                        if (stride_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY;
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    }
                  } else {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    }
                    if (num_any_direction[X]==1) {
                      if (num_any_direction[Y]==1) {
                        const int i = yee_idx;
                        const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                           + np->energy[Ez][i];
                        np->s[ec][i] = energy_here*fac;
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      } else {
                        if (stride_any_direction[Y]==1) {
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iY;
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      }
                    } else {
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                             + np->energy[Ez][i];
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      } else {
                        if (stride_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY;
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    }
                  }
                } else {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                }
              }
            }
          } else {
            if (pol) {
              for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                const double fac = np->saturation_factor;
                const double g = op->pb->gamma;
                const double om = op->pb->omeganot;
                const double invomsqr = 1.0/(om*om);
                const double funinv = 1.0/(1+0.5*g);
                if (fac) {
                  if (fac > 0) {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                      np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    }
                    if (num_any_direction[X]==1) {
                      if (num_any_direction[Y]==1) {
                        const int i = yee_idx;
                        const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                           + np->energy[Ez][i];
                        np->s[ec][i] = max(-energy_here*fac, 0.0);
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                        op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                           + np->s[ec][i]*f[ec][1][i]);
                      } else {
                        if (stride_any_direction[Y]==1) {
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iY;
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      }
                    } else {
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                             + np->energy[Ez][i];
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        }
                      } else {
                        if (stride_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY;
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    }
                  } else {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                      np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    }
                    if (num_any_direction[X]==1) {
                      if (num_any_direction[Y]==1) {
                        const int i = yee_idx;
                        const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                           + np->energy[Ez][i];
                        np->s[ec][i] = energy_here*fac;
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                        op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                           + np->s[ec][i]*f[ec][1][i]);
                      } else {
                        if (stride_any_direction[Y]==1) {
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iY;
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                               + np->energy[Ez][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      }
                    } else {
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                             + np->energy[Ez][i];
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        }
                      } else {
                        if (stride_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY;
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[Ex][i]) + (np->energy[Ey][i])
                                 + np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    }
                  }
                } else {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  }
                }
              }
            }
          }
        }
        // The polarizations got switched...
        polarization *temp = olpol;
        olpol = pol;
        pol = temp;
      } else {
        FOR_ELECTRIC_COMPONENTS(ec) if (f[ec][0]) {
          const int yee_idx = v.yee_index(ec);
          if (is_real) {
            if (pol) {
              for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                const double fac = np->saturation_factor;
                const double g = op->pb->gamma;
                const double om = op->pb->omeganot;
                const double invomsqr = 1.0/(om*om);
                const double funinv = 1.0/(1+0.5*g);
                if (fac) {
                  if (fac > 0) {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    }
                    if (num_any_direction[X]==1) {
                      if (num_any_direction[Y]==1) {
                        const int i = yee_idx;
                        const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                        np->s[ec][i] = max(-energy_here*fac, 0.0);
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      } else {
                        if (stride_any_direction[Y]==1) {
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iY;
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      }
                    } else {
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      } else {
                        if (stride_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY;
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    }
                  } else {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    }
                    if (num_any_direction[X]==1) {
                      if (num_any_direction[Y]==1) {
                        const int i = yee_idx;
                        const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                        np->s[ec][i] = energy_here*fac;
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      } else {
                        if (stride_any_direction[Y]==1) {
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iY;
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      }
                    } else {
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      } else {
                        if (stride_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY;
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    }
                  }
                } else {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                }
              }
            }
          } else {
            if (pol) {
              for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                const double fac = np->saturation_factor;
                const double g = op->pb->gamma;
                const double om = op->pb->omeganot;
                const double invomsqr = 1.0/(om*om);
                const double funinv = 1.0/(1+0.5*g);
                if (fac) {
                  if (fac > 0) {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                      np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    }
                    if (num_any_direction[X]==1) {
                      if (num_any_direction[Y]==1) {
                        const int i = yee_idx;
                        const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                        np->s[ec][i] = max(-energy_here*fac, 0.0);
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                        op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                           + np->s[ec][i]*f[ec][1][i]);
                      } else {
                        if (stride_any_direction[Y]==1) {
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iY;
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      }
                    } else {
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        }
                      } else {
                        if (stride_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY;
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    }
                  } else {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                      np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    }
                    if (num_any_direction[X]==1) {
                      if (num_any_direction[Y]==1) {
                        const int i = yee_idx;
                        const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                        np->s[ec][i] = energy_here*fac;
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                        op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                           + np->s[ec][i]*f[ec][1][i]);
                      } else {
                        if (stride_any_direction[Y]==1) {
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iY;
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      }
                    } else {
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        }
                      } else {
                        if (stride_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY;
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[Ex][i]) + np->energy[Ey][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    }
                  }
                } else {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  }
                }
              }
            }
          }
        }
        // The polarizations got switched...
        polarization *temp = olpol;
        olpol = pol;
        pol = temp;
      }
    } else {
      if (f[Ez][0]) {
        FOR_ELECTRIC_COMPONENTS(ec) if (f[ec][0]) {
          const int yee_idx = v.yee_index(ec);
          if (is_real) {
            if (pol) {
              for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                const double fac = np->saturation_factor;
                const double g = op->pb->gamma;
                const double om = op->pb->omeganot;
                const double invomsqr = 1.0/(om*om);
                const double funinv = 1.0/(1+0.5*g);
                if (fac) {
                  if (fac > 0) {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    }
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              const int i = yee_idx;
                              const double energy_here = np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            } else {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                const int i = yee_idx + iY*stride_any_direction[Y];
                                const double energy_here = np->energy[Ez][i];
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            }
                          } else {
                            const int i = yee_idx;
                            const double energy_here = np->energy[Ez][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              const int i = yee_idx;
                              const double energy_here = np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                  const int i = yee_idx + iY;
                                  const double energy_here = np->energy[Ez][i];
                                  np->s[ec][i] = max(-energy_here*fac, 0.0);
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                }
                              } else {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  const int i = yee_idx + iY*stride_any_direction[Y];
                                  const double energy_here = np->energy[Ez][i];
                                  np->s[ec][i] = max(-energy_here*fac, 0.0);
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                }
                              }
                            }
                          } else {
                            const int i = yee_idx;
                            const double energy_here = np->energy[Ez][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      } else {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                const int i = yee_idx + iX;
                                const double energy_here = np->energy[Ez][i];
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            } else {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  const int i = yee_idx + iX + iY*stride_any_direction[Y];
                                  const double energy_here = np->energy[Ez][i];
                                  np->s[ec][i] = max(-energy_here*fac, 0.0);
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                              const int i = yee_idx + iX;
                              const double energy_here = np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                const int i = yee_idx + iX*stride_any_direction[X];
                                const double energy_here = np->energy[Ez][i];
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                    const int i = yee_idx + iX*stride_any_direction[X] + iY;
                                    const double energy_here = np->energy[Ez][i];
                                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                       + np->s[ec][i]*f[ec][0][i]);
                                  }
                                }
                              } else {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                    const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                                    const double energy_here = np->energy[Ez][i];
                                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                       + np->s[ec][i]*f[ec][0][i]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              const int i = yee_idx + iX*stride_any_direction[X];
                              const double energy_here = np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    } else {
                      const int i = yee_idx;
                      const double energy_here = np->energy[Ez][i];
                      np->s[ec][i] = max(-energy_here*fac, 0.0);
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    }
                  } else {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    }
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              const int i = yee_idx;
                              const double energy_here = np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            } else {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                const int i = yee_idx + iY*stride_any_direction[Y];
                                const double energy_here = np->energy[Ez][i];
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            }
                          } else {
                            const int i = yee_idx;
                            const double energy_here = np->energy[Ez][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              const int i = yee_idx;
                              const double energy_here = np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                  const int i = yee_idx + iY;
                                  const double energy_here = np->energy[Ez][i];
                                  np->s[ec][i] = energy_here*fac;
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                }
                              } else {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  const int i = yee_idx + iY*stride_any_direction[Y];
                                  const double energy_here = np->energy[Ez][i];
                                  np->s[ec][i] = energy_here*fac;
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                }
                              }
                            }
                          } else {
                            const int i = yee_idx;
                            const double energy_here = np->energy[Ez][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      } else {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                const int i = yee_idx + iX;
                                const double energy_here = np->energy[Ez][i];
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            } else {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  const int i = yee_idx + iX + iY*stride_any_direction[Y];
                                  const double energy_here = np->energy[Ez][i];
                                  np->s[ec][i] = energy_here*fac;
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                              const int i = yee_idx + iX;
                              const double energy_here = np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                const int i = yee_idx + iX*stride_any_direction[X];
                                const double energy_here = np->energy[Ez][i];
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                    const int i = yee_idx + iX*stride_any_direction[X] + iY;
                                    const double energy_here = np->energy[Ez][i];
                                    np->s[ec][i] = energy_here*fac;
                                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                       + np->s[ec][i]*f[ec][0][i]);
                                  }
                                }
                              } else {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                    const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                                    const double energy_here = np->energy[Ez][i];
                                    np->s[ec][i] = energy_here*fac;
                                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                       + np->s[ec][i]*f[ec][0][i]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              const int i = yee_idx + iX*stride_any_direction[X];
                              const double energy_here = np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    } else {
                      const int i = yee_idx;
                      const double energy_here = np->energy[Ez][i];
                      np->s[ec][i] = energy_here*fac;
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    }
                  }
                } else {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                }
              }
            }
          } else {
            if (pol) {
              for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                const double fac = np->saturation_factor;
                const double g = op->pb->gamma;
                const double om = op->pb->omeganot;
                const double invomsqr = 1.0/(om*om);
                const double funinv = 1.0/(1+0.5*g);
                if (fac) {
                  if (fac > 0) {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                      np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    }
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              const int i = yee_idx;
                              const double energy_here = np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            } else {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                const int i = yee_idx + iY*stride_any_direction[Y];
                                const double energy_here = np->energy[Ez][i];
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                                op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                   + np->s[ec][i]*f[ec][1][i]);
                              }
                            }
                          } else {
                            const int i = yee_idx;
                            const double energy_here = np->energy[Ez][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              const int i = yee_idx;
                              const double energy_here = np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                  const int i = yee_idx + iY;
                                  const double energy_here = np->energy[Ez][i];
                                  np->s[ec][i] = max(-energy_here*fac, 0.0);
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                     + np->s[ec][i]*f[ec][1][i]);
                                }
                              } else {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  const int i = yee_idx + iY*stride_any_direction[Y];
                                  const double energy_here = np->energy[Ez][i];
                                  np->s[ec][i] = max(-energy_here*fac, 0.0);
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                     + np->s[ec][i]*f[ec][1][i]);
                                }
                              }
                            }
                          } else {
                            const int i = yee_idx;
                            const double energy_here = np->energy[Ez][i];
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      } else {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                const int i = yee_idx + iX;
                                const double energy_here = np->energy[Ez][i];
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                                op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                   + np->s[ec][i]*f[ec][1][i]);
                              }
                            } else {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  const int i = yee_idx + iX + iY*stride_any_direction[Y];
                                  const double energy_here = np->energy[Ez][i];
                                  np->s[ec][i] = max(-energy_here*fac, 0.0);
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                     + np->s[ec][i]*f[ec][1][i]);
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                              const int i = yee_idx + iX;
                              const double energy_here = np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                const int i = yee_idx + iX*stride_any_direction[X];
                                const double energy_here = np->energy[Ez][i];
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                                op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                   + np->s[ec][i]*f[ec][1][i]);
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                    const int i = yee_idx + iX*stride_any_direction[X] + iY;
                                    const double energy_here = np->energy[Ez][i];
                                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                       + np->s[ec][i]*f[ec][0][i]);
                                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                       + np->s[ec][i]*f[ec][1][i]);
                                  }
                                }
                              } else {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                    const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                                    const double energy_here = np->energy[Ez][i];
                                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                       + np->s[ec][i]*f[ec][0][i]);
                                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                       + np->s[ec][i]*f[ec][1][i]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              const int i = yee_idx + iX*stride_any_direction[X];
                              const double energy_here = np->energy[Ez][i];
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    } else {
                      const int i = yee_idx;
                      const double energy_here = np->energy[Ez][i];
                      np->s[ec][i] = max(-energy_here*fac, 0.0);
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                      op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                         + np->s[ec][i]*f[ec][1][i]);
                    }
                  } else {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                      np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    }
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              const int i = yee_idx;
                              const double energy_here = np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            } else {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                const int i = yee_idx + iY*stride_any_direction[Y];
                                const double energy_here = np->energy[Ez][i];
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                                op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                   + np->s[ec][i]*f[ec][1][i]);
                              }
                            }
                          } else {
                            const int i = yee_idx;
                            const double energy_here = np->energy[Ez][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              const int i = yee_idx;
                              const double energy_here = np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                  const int i = yee_idx + iY;
                                  const double energy_here = np->energy[Ez][i];
                                  np->s[ec][i] = energy_here*fac;
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                     + np->s[ec][i]*f[ec][1][i]);
                                }
                              } else {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  const int i = yee_idx + iY*stride_any_direction[Y];
                                  const double energy_here = np->energy[Ez][i];
                                  np->s[ec][i] = energy_here*fac;
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                     + np->s[ec][i]*f[ec][1][i]);
                                }
                              }
                            }
                          } else {
                            const int i = yee_idx;
                            const double energy_here = np->energy[Ez][i];
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      } else {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                const int i = yee_idx + iX;
                                const double energy_here = np->energy[Ez][i];
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                                op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                   + np->s[ec][i]*f[ec][1][i]);
                              }
                            } else {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  const int i = yee_idx + iX + iY*stride_any_direction[Y];
                                  const double energy_here = np->energy[Ez][i];
                                  np->s[ec][i] = energy_here*fac;
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                     + np->s[ec][i]*f[ec][1][i]);
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                              const int i = yee_idx + iX;
                              const double energy_here = np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                const int i = yee_idx + iX*stride_any_direction[X];
                                const double energy_here = np->energy[Ez][i];
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                                op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                   + np->s[ec][i]*f[ec][1][i]);
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                    const int i = yee_idx + iX*stride_any_direction[X] + iY;
                                    const double energy_here = np->energy[Ez][i];
                                    np->s[ec][i] = energy_here*fac;
                                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                       + np->s[ec][i]*f[ec][0][i]);
                                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                       + np->s[ec][i]*f[ec][1][i]);
                                  }
                                }
                              } else {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                    const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                                    const double energy_here = np->energy[Ez][i];
                                    np->s[ec][i] = energy_here*fac;
                                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                       + np->s[ec][i]*f[ec][0][i]);
                                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                       + np->s[ec][i]*f[ec][1][i]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              const int i = yee_idx + iX*stride_any_direction[X];
                              const double energy_here = np->energy[Ez][i];
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    } else {
                      const int i = yee_idx;
                      const double energy_here = np->energy[Ez][i];
                      np->s[ec][i] = energy_here*fac;
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                      op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                         + np->s[ec][i]*f[ec][1][i]);
                    }
                  }
                } else {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  }
                }
              }
            }
          }
        }
        // The polarizations got switched...
        polarization *temp = olpol;
        olpol = pol;
        pol = temp;
      } else {
        FOR_ELECTRIC_COMPONENTS(ec) if (f[ec][0]) {
          const int yee_idx = v.yee_index(ec);
          if (is_real) {
            if (pol) {
              for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                const double fac = np->saturation_factor;
                const double g = op->pb->gamma;
                const double om = op->pb->omeganot;
                const double invomsqr = 1.0/(om*om);
                const double funinv = 1.0/(1+0.5*g);
                if (fac) {
                  if (fac > 0) {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    }
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              const int i = yee_idx;
                              const double energy_here = 0;
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            } else {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                const int i = yee_idx + iY*stride_any_direction[Y];
                                const double energy_here = 0;
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            }
                          } else {
                            const int i = yee_idx;
                            const double energy_here = 0;
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              const int i = yee_idx;
                              const double energy_here = 0;
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                  const int i = yee_idx + iY;
                                  const double energy_here = 0;
                                  np->s[ec][i] = max(-energy_here*fac, 0.0);
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                }
                              } else {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  const int i = yee_idx + iY*stride_any_direction[Y];
                                  const double energy_here = 0;
                                  np->s[ec][i] = max(-energy_here*fac, 0.0);
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                }
                              }
                            }
                          } else {
                            const int i = yee_idx;
                            const double energy_here = 0;
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      } else {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                const int i = yee_idx + iX;
                                const double energy_here = 0;
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            } else {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  const int i = yee_idx + iX + iY*stride_any_direction[Y];
                                  const double energy_here = 0;
                                  np->s[ec][i] = max(-energy_here*fac, 0.0);
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                              const int i = yee_idx + iX;
                              const double energy_here = 0;
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                const int i = yee_idx + iX*stride_any_direction[X];
                                const double energy_here = 0;
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                    const int i = yee_idx + iX*stride_any_direction[X] + iY;
                                    const double energy_here = 0;
                                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                       + np->s[ec][i]*f[ec][0][i]);
                                  }
                                }
                              } else {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                    const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                                    const double energy_here = 0;
                                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                       + np->s[ec][i]*f[ec][0][i]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              const int i = yee_idx + iX*stride_any_direction[X];
                              const double energy_here = 0;
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    } else {
                      const int i = yee_idx;
                      const double energy_here = 0;
                      np->s[ec][i] = max(-energy_here*fac, 0.0);
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    }
                  } else {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    }
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              const int i = yee_idx;
                              const double energy_here = 0;
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            } else {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                const int i = yee_idx + iY*stride_any_direction[Y];
                                const double energy_here = 0;
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            }
                          } else {
                            const int i = yee_idx;
                            const double energy_here = 0;
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              const int i = yee_idx;
                              const double energy_here = 0;
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                  const int i = yee_idx + iY;
                                  const double energy_here = 0;
                                  np->s[ec][i] = energy_here*fac;
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                }
                              } else {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  const int i = yee_idx + iY*stride_any_direction[Y];
                                  const double energy_here = 0;
                                  np->s[ec][i] = energy_here*fac;
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                }
                              }
                            }
                          } else {
                            const int i = yee_idx;
                            const double energy_here = 0;
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      } else {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                const int i = yee_idx + iX;
                                const double energy_here = 0;
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            } else {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  const int i = yee_idx + iX + iY*stride_any_direction[Y];
                                  const double energy_here = 0;
                                  np->s[ec][i] = energy_here*fac;
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                              const int i = yee_idx + iX;
                              const double energy_here = 0;
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                const int i = yee_idx + iX*stride_any_direction[X];
                                const double energy_here = 0;
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                    const int i = yee_idx + iX*stride_any_direction[X] + iY;
                                    const double energy_here = 0;
                                    np->s[ec][i] = energy_here*fac;
                                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                       + np->s[ec][i]*f[ec][0][i]);
                                  }
                                }
                              } else {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                    const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                                    const double energy_here = 0;
                                    np->s[ec][i] = energy_here*fac;
                                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                       + np->s[ec][i]*f[ec][0][i]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              const int i = yee_idx + iX*stride_any_direction[X];
                              const double energy_here = 0;
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    } else {
                      const int i = yee_idx;
                      const double energy_here = 0;
                      np->s[ec][i] = energy_here*fac;
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    }
                  }
                } else {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                }
              }
            }
          } else {
            if (pol) {
              for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                const double fac = np->saturation_factor;
                const double g = op->pb->gamma;
                const double om = op->pb->omeganot;
                const double invomsqr = 1.0/(om*om);
                const double funinv = 1.0/(1+0.5*g);
                if (fac) {
                  if (fac > 0) {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                      np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    }
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              const int i = yee_idx;
                              const double energy_here = 0;
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            } else {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                const int i = yee_idx + iY*stride_any_direction[Y];
                                const double energy_here = 0;
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                                op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                   + np->s[ec][i]*f[ec][1][i]);
                              }
                            }
                          } else {
                            const int i = yee_idx;
                            const double energy_here = 0;
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              const int i = yee_idx;
                              const double energy_here = 0;
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                  const int i = yee_idx + iY;
                                  const double energy_here = 0;
                                  np->s[ec][i] = max(-energy_here*fac, 0.0);
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                     + np->s[ec][i]*f[ec][1][i]);
                                }
                              } else {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  const int i = yee_idx + iY*stride_any_direction[Y];
                                  const double energy_here = 0;
                                  np->s[ec][i] = max(-energy_here*fac, 0.0);
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                     + np->s[ec][i]*f[ec][1][i]);
                                }
                              }
                            }
                          } else {
                            const int i = yee_idx;
                            const double energy_here = 0;
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      } else {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                const int i = yee_idx + iX;
                                const double energy_here = 0;
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                                op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                   + np->s[ec][i]*f[ec][1][i]);
                              }
                            } else {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  const int i = yee_idx + iX + iY*stride_any_direction[Y];
                                  const double energy_here = 0;
                                  np->s[ec][i] = max(-energy_here*fac, 0.0);
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                     + np->s[ec][i]*f[ec][1][i]);
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                              const int i = yee_idx + iX;
                              const double energy_here = 0;
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                const int i = yee_idx + iX*stride_any_direction[X];
                                const double energy_here = 0;
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                                op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                   + np->s[ec][i]*f[ec][1][i]);
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                    const int i = yee_idx + iX*stride_any_direction[X] + iY;
                                    const double energy_here = 0;
                                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                       + np->s[ec][i]*f[ec][0][i]);
                                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                       + np->s[ec][i]*f[ec][1][i]);
                                  }
                                }
                              } else {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                    const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                                    const double energy_here = 0;
                                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                       + np->s[ec][i]*f[ec][0][i]);
                                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                       + np->s[ec][i]*f[ec][1][i]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              const int i = yee_idx + iX*stride_any_direction[X];
                              const double energy_here = 0;
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    } else {
                      const int i = yee_idx;
                      const double energy_here = 0;
                      np->s[ec][i] = max(-energy_here*fac, 0.0);
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                      op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                         + np->s[ec][i]*f[ec][1][i]);
                    }
                  } else {
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                      np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    }
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              const int i = yee_idx;
                              const double energy_here = 0;
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            } else {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                const int i = yee_idx + iY*stride_any_direction[Y];
                                const double energy_here = 0;
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                                op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                   + np->s[ec][i]*f[ec][1][i]);
                              }
                            }
                          } else {
                            const int i = yee_idx;
                            const double energy_here = 0;
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              const int i = yee_idx;
                              const double energy_here = 0;
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                  const int i = yee_idx + iY;
                                  const double energy_here = 0;
                                  np->s[ec][i] = energy_here*fac;
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                     + np->s[ec][i]*f[ec][1][i]);
                                }
                              } else {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  const int i = yee_idx + iY*stride_any_direction[Y];
                                  const double energy_here = 0;
                                  np->s[ec][i] = energy_here*fac;
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                     + np->s[ec][i]*f[ec][1][i]);
                                }
                              }
                            }
                          } else {
                            const int i = yee_idx;
                            const double energy_here = 0;
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      } else {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                const int i = yee_idx + iX;
                                const double energy_here = 0;
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                                op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                   + np->s[ec][i]*f[ec][1][i]);
                              }
                            } else {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  const int i = yee_idx + iX + iY*stride_any_direction[Y];
                                  const double energy_here = 0;
                                  np->s[ec][i] = energy_here*fac;
                                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                     + np->s[ec][i]*f[ec][0][i]);
                                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                     + np->s[ec][i]*f[ec][1][i]);
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                              const int i = yee_idx + iX;
                              const double energy_here = 0;
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                const int i = yee_idx + iX*stride_any_direction[X];
                                const double energy_here = 0;
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                                op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                   + np->s[ec][i]*f[ec][1][i]);
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                    const int i = yee_idx + iX*stride_any_direction[X] + iY;
                                    const double energy_here = 0;
                                    np->s[ec][i] = energy_here*fac;
                                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                       + np->s[ec][i]*f[ec][0][i]);
                                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                       + np->s[ec][i]*f[ec][1][i]);
                                  }
                                }
                              } else {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                    const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                                    const double energy_here = 0;
                                    np->s[ec][i] = energy_here*fac;
                                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                       + np->s[ec][i]*f[ec][0][i]);
                                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                       + np->s[ec][i]*f[ec][1][i]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              const int i = yee_idx + iX*stride_any_direction[X];
                              const double energy_here = 0;
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    } else {
                      const int i = yee_idx;
                      const double energy_here = 0;
                      np->s[ec][i] = energy_here*fac;
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                      op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                         + np->s[ec][i]*f[ec][1][i]);
                    }
                  }
                } else {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  }
                }
              }
            }
          }
        }
        // The polarizations got switched...
        polarization *temp = olpol;
        olpol = pol;
        pol = temp;
      }
    }
  }
}
