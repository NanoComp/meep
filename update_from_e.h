FOR_E_AND_D(ec,dc) if (f[ec][0]) {
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
            if (stride_any_direction[Z]) {
              if (num_any_direction[Z]==1) {
                if (stride_any_direction[Z]==1) {
                  if (stride_any_direction[R]) {
                    if (num_any_direction[R]==1) {
                      np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                      op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                         + np->s[ec][0]*f[ec][0][0]);
                    } else {
                      for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                        np->s[ec][iR*stride_any_direction[R]] = max(-np->energy[ec][iR*stride_any_direction[R]]*fac,
                           0.0);
                        op->P[ec][0][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][0][iR*stride_any_direction[R]]
                           + (0.5*g-1)*op->P[ec][0][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][0][iR*stride_any_direction[R]]);
                      }
                    }
                  } else {
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[Y]) {
                          if (num_any_direction[Y]==1) {
                            np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                            op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                               + np->s[ec][0]*f[ec][0][0]);
                          } else {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              np->s[ec][iY*stride_any_direction[Y]] = max(-np->energy[ec][iY*stride_any_direction[Y]]*fac,
                                 0.0);
                              op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                                 + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                            }
                          }
                        } else {
                          np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                          op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                             + np->s[ec][0]*f[ec][0][0]);
                        }
                      } else {
                        if (stride_any_direction[Y]) {
                          if (num_any_direction[Y]==1) {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              np->s[ec][iX*stride_any_direction[X]] = max(-np->energy[ec][iX*stride_any_direction[X]]*fac,
                                 0.0);
                              op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                                 + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = max(-np->energy[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]*fac,
                                      0.0);
                                op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                      + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                         + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y]]*f[ec][0][iX*stride_any_direction[X]
                                               + iY*stride_any_direction[Y]]);
                              }
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            np->s[ec][iX*stride_any_direction[X]] = max(-np->energy[ec][iX*stride_any_direction[X]]*fac,
                               0.0);
                            op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                               + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                          }
                        }
                      }
                    } else {
                      np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                      op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                         + np->s[ec][0]*f[ec][0][0]);
                    }
                  }
                } else {
                  if (stride_any_direction[R]) {
                    if (num_any_direction[R]==1) {
                      np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                      op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                         + np->s[ec][0]*f[ec][0][0]);
                    } else {
                      if (stride_any_direction[R]==1) {
                        for (int iR=0; iR<num_any_direction[R]; iR += 1) {
                          np->s[ec][iR] = max(-np->energy[ec][iR]*fac, 0.0);
                          op->P[ec][0][iR] = funinv*((2-om*om)*np->P[ec][0][iR] + (0.5*g-1)*op->P[ec][0][iR] + np->s[ec][iR]*f[ec][0][iR]);
                        }
                      } else {
                        for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                          np->s[ec][iR*stride_any_direction[R]] = max(-np->energy[ec][iR*stride_any_direction[R]]*fac,
                             0.0);
                          op->P[ec][0][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][0][iR*stride_any_direction[R]]
                             + (0.5*g-1)*op->P[ec][0][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][0][iR*stride_any_direction[R]]);
                        }
                      }
                    }
                  } else {
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                              op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                                 + np->s[ec][0]*f[ec][0][0]);
                            } else {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                np->s[ec][iY*stride_any_direction[Y]] = max(-np->energy[ec][iY*stride_any_direction[Y]]*fac,
                                   0.0);
                                op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                                   + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                              }
                            }
                          } else {
                            np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                            op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                               + np->s[ec][0]*f[ec][0][0]);
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                              op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                                 + np->s[ec][0]*f[ec][0][0]);
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                  np->s[ec][iY] = max(-np->energy[ec][iY]*fac, 0.0);
                                  op->P[ec][0][iY] = funinv*((2-om*om)*np->P[ec][0][iY] + (0.5*g-1)*op->P[ec][0][iY] + np->s[ec][iY]*f[ec][0][iY]);
                                }
                              } else {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  np->s[ec][iY*stride_any_direction[Y]] = max(-np->energy[ec][iY*stride_any_direction[Y]]*fac,
                                     0.0);
                                  op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                                     + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                                }
                              }
                            }
                          } else {
                            np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                            op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                               + np->s[ec][0]*f[ec][0][0]);
                          }
                        }
                      } else {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                np->s[ec][iX] = max(-np->energy[ec][iX]*fac, 0.0);
                                op->P[ec][0][iX] = funinv*((2-om*om)*np->P[ec][0][iX] + (0.5*g-1)*op->P[ec][0][iX] + np->s[ec][iX]*f[ec][0][iX]);
                              }
                            } else {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  np->s[ec][iX + iY*stride_any_direction[Y]] = max(-np->energy[ec][iX
                                     + iY*stride_any_direction[Y]]*fac, 0.0);
                                  op->P[ec][0][iX + iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iX
                                     + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX
                                        + iY*stride_any_direction[Y]] + np->s[ec][iX + iY*stride_any_direction[Y]]*f[ec][0][iX
                                           + iY*stride_any_direction[Y]]);
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                              np->s[ec][iX] = max(-np->energy[ec][iX]*fac, 0.0);
                              op->P[ec][0][iX] = funinv*((2-om*om)*np->P[ec][0][iX] + (0.5*g-1)*op->P[ec][0][iX] + np->s[ec][iX]*f[ec][0][iX]);
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                np->s[ec][iX*stride_any_direction[X]] = max(-np->energy[ec][iX*stride_any_direction[X]]*fac,
                                   0.0);
                                op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                                   + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                    np->s[ec][iX*stride_any_direction[X] + iY] = max(-np->energy[ec][iX*stride_any_direction[X]
                                       + iY]*fac, 0.0);
                                    op->P[ec][0][iX*stride_any_direction[X] + iY] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                       + iY] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                          + iY] + np->s[ec][iX*stride_any_direction[X] + iY]*f[ec][0][iX*stride_any_direction[X]
                                             + iY]);
                                  }
                                }
                              } else {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                    np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                       = max(-np->energy[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]*fac,
                                          0.0);
                                    op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                       = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                          + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                             + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                                + iY*stride_any_direction[Y]]*f[ec][0][iX*stride_any_direction[X]
                                                   + iY*stride_any_direction[Y]]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              np->s[ec][iX*stride_any_direction[X]] = max(-np->energy[ec][iX*stride_any_direction[X]]*fac,
                                 0.0);
                              op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                                 + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                            }
                          }
                        }
                      }
                    } else {
                      np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                      op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                         + np->s[ec][0]*f[ec][0][0]);
                    }
                  }
                }
              } else {
                if (stride_any_direction[Z]==1) {
                  if (stride_any_direction[R]) {
                    if (num_any_direction[R]==1) {
                      for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                        np->s[ec][iZ] = max(-np->energy[ec][iZ]*fac, 0.0);
                        op->P[ec][0][iZ] = funinv*((2-om*om)*np->P[ec][0][iZ] + (0.5*g-1)*op->P[ec][0][iZ] + np->s[ec][iZ]*f[ec][0][iZ]);
                      }
                    } else {
                      for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                        for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                          np->s[ec][iZ + iR*stride_any_direction[R]] = max(-np->energy[ec][iZ
                             + iR*stride_any_direction[R]]*fac, 0.0);
                          op->P[ec][0][iZ + iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][0][iZ
                             + iR*stride_any_direction[R]] + (0.5*g-1)*op->P[ec][0][iZ
                                + iR*stride_any_direction[R]] + np->s[ec][iZ + iR*stride_any_direction[R]]*f[ec][0][iZ
                                   + iR*stride_any_direction[R]]);
                        }
                      }
                    }
                  } else {
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[Y]) {
                          if (num_any_direction[Y]==1) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              np->s[ec][iZ] = max(-np->energy[ec][iZ]*fac, 0.0);
                              op->P[ec][0][iZ] = funinv*((2-om*om)*np->P[ec][0][iZ] + (0.5*g-1)*op->P[ec][0][iZ] + np->s[ec][iZ]*f[ec][0][iZ]);
                            }
                          } else {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                np->s[ec][iY*stride_any_direction[Y] + iZ] = max(-np->energy[ec][iY*stride_any_direction[Y]
                                   + iZ]*fac, 0.0);
                                op->P[ec][0][iY*stride_any_direction[Y] + iZ] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]
                                   + iZ] + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]
                                      + iZ] + np->s[ec][iY*stride_any_direction[Y] + iZ]*f[ec][0][iY*stride_any_direction[Y]
                                         + iZ]);
                              }
                            }
                          }
                        } else {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                            np->s[ec][iZ] = max(-np->energy[ec][iZ]*fac, 0.0);
                            op->P[ec][0][iZ] = funinv*((2-om*om)*np->P[ec][0][iZ] + (0.5*g-1)*op->P[ec][0][iZ] + np->s[ec][iZ]*f[ec][0][iZ]);
                          }
                        }
                      } else {
                        if (stride_any_direction[Y]) {
                          if (num_any_direction[Y]==1) {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                np->s[ec][iX*stride_any_direction[X] + iZ] = max(-np->energy[ec][iX*stride_any_direction[X]
                                   + iZ]*fac, 0.0);
                                op->P[ec][0][iX*stride_any_direction[X] + iZ] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                   + iZ] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                      + iZ] + np->s[ec][iX*stride_any_direction[X] + iZ]*f[ec][0][iX*stride_any_direction[X]
                                         + iZ]);
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                  np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                     + iZ] = max(-np->energy[ec][iX*stride_any_direction[X] +
                                        iY*stride_any_direction[Y] + iZ]*fac, 0.0);
                                  op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                     + iZ] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                        + iY*stride_any_direction[Y] + iZ] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                           + iY*stride_any_direction[Y] + iZ] + np->s[ec][iX*stride_any_direction[X]
                                              + iY*stride_any_direction[Y] + iZ]*f[ec][0][iX*stride_any_direction[X]
                                                 + iY*stride_any_direction[Y] + iZ]);
                                }
                              }
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              np->s[ec][iX*stride_any_direction[X] + iZ] = max(-np->energy[ec][iX*stride_any_direction[X]
                                 + iZ]*fac, 0.0);
                              op->P[ec][0][iX*stride_any_direction[X] + iZ] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                 + iZ] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                    + iZ] + np->s[ec][iX*stride_any_direction[X] + iZ]*f[ec][0][iX*stride_any_direction[X]
                                       + iZ]);
                            }
                          }
                        }
                      }
                    } else {
                      for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                        np->s[ec][iZ] = max(-np->energy[ec][iZ]*fac, 0.0);
                        op->P[ec][0][iZ] = funinv*((2-om*om)*np->P[ec][0][iZ] + (0.5*g-1)*op->P[ec][0][iZ] + np->s[ec][iZ]*f[ec][0][iZ]);
                      }
                    }
                  }
                } else {
                  if (stride_any_direction[R]) {
                    if (num_any_direction[R]==1) {
                      for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                        np->s[ec][iZ*stride_any_direction[Z]] = max(-np->energy[ec][iZ*stride_any_direction[Z]]*fac,
                           0.0);
                        op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                           + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                      }
                    } else {
                      if (stride_any_direction[R]==1) {
                        for (int iR=0; iR<num_any_direction[R]; iR += 1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                            np->s[ec][iZ*stride_any_direction[Z] + iR] = max(-np->energy[ec][iZ*stride_any_direction[Z]
                               + iR]*fac, 0.0);
                            op->P[ec][0][iZ*stride_any_direction[Z] + iR] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]
                               + iR] + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]
                                  + iR] + np->s[ec][iZ*stride_any_direction[Z] + iR]*f[ec][0][iZ*stride_any_direction[Z]
                                     + iR]);
                          }
                        }
                      } else {
                        for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                            np->s[ec][iZ*stride_any_direction[Z] + iR*stride_any_direction[R]]
                               = max(-np->energy[ec][iZ*stride_any_direction[Z] + iR*stride_any_direction[R]]*fac,
                                  0.0);
                            op->P[ec][0][iZ*stride_any_direction[Z] + iR*stride_any_direction[R]]
                               = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]
                                  + iR*stride_any_direction[R]] + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]
                                     + iR*stride_any_direction[R]] + np->s[ec][iZ*stride_any_direction[Z]
                                        + iR*stride_any_direction[R]]*f[ec][0][iZ*stride_any_direction[Z]
                                           + iR*stride_any_direction[R]]);
                          }
                        }
                      }
                    }
                  } else {
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                np->s[ec][iZ*stride_any_direction[Z]] = max(-np->energy[ec][iZ*stride_any_direction[Z]]*fac,
                                   0.0);
                                op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                                   + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                              }
                            } else {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                  np->s[ec][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                     = max(-np->energy[ec][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*fac,
                                        0.0);
                                  op->P[ec][0][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                     = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]
                                        + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]
                                           + iZ*stride_any_direction[Z]] + np->s[ec][iY*stride_any_direction[Y]
                                              + iZ*stride_any_direction[Z]]*f[ec][0][iY*stride_any_direction[Y]
                                                 + iZ*stride_any_direction[Z]]);
                                }
                              }
                            }
                          } else {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                              np->s[ec][iZ*stride_any_direction[Z]] = max(-np->energy[ec][iZ*stride_any_direction[Z]]*fac,
                                 0.0);
                              op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                                 + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                np->s[ec][iZ*stride_any_direction[Z]] = max(-np->energy[ec][iZ*stride_any_direction[Z]]*fac,
                                   0.0);
                                op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                                   + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                    np->s[ec][iY + iZ*stride_any_direction[Z]] = max(-np->energy[ec][iY
                                       + iZ*stride_any_direction[Z]]*fac, 0.0);
                                    op->P[ec][0][iY + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iY
                                       + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iY
                                          + iZ*stride_any_direction[Z]] + np->s[ec][iY + iZ*stride_any_direction[Z]]*f[ec][0][iY
                                             + iZ*stride_any_direction[Z]]);
                                  }
                                }
                              } else {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                    np->s[ec][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = max(-np->energy[ec][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*fac,
                                          0.0);
                                    op->P[ec][0][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]
                                          + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]
                                             + iZ*stride_any_direction[Z]] + np->s[ec][iY*stride_any_direction[Y]
                                                + iZ*stride_any_direction[Z]]*f[ec][0][iY*stride_any_direction[Y]
                                                   + iZ*stride_any_direction[Z]]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                              np->s[ec][iZ*stride_any_direction[Z]] = max(-np->energy[ec][iZ*stride_any_direction[Z]]*fac,
                                 0.0);
                              op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                                 + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                            }
                          }
                        }
                      } else {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                  np->s[ec][iX + iZ*stride_any_direction[Z]] = max(-np->energy[ec][iX
                                     + iZ*stride_any_direction[Z]]*fac, 0.0);
                                  op->P[ec][0][iX + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iX
                                     + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX
                                        + iZ*stride_any_direction[Z]] + np->s[ec][iX + iZ*stride_any_direction[Z]]*f[ec][0][iX
                                           + iZ*stride_any_direction[Z]]);
                                }
                              }
                            } else {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                    np->s[ec][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = max(-np->energy[ec][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*fac, 0.0);
                                    op->P[ec][0][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = funinv*((2-om*om)*np->P[ec][0][iX + iY*stride_any_direction[Y]
                                          + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX
                                             + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                                + np->s[ec][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*f[ec][0][iX
                                                   + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                np->s[ec][iX + iZ*stride_any_direction[Z]] = max(-np->energy[ec][iX
                                   + iZ*stride_any_direction[Z]]*fac, 0.0);
                                op->P[ec][0][iX + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iX
                                   + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX
                                      + iZ*stride_any_direction[Z]] + np->s[ec][iX + iZ*stride_any_direction[Z]]*f[ec][0][iX
                                         + iZ*stride_any_direction[Z]]);
                              }
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                  np->s[ec][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                     = max(-np->energy[ec][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]*fac,
                                        0.0);
                                  op->P[ec][0][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                     = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                        + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                           + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                              + iZ*stride_any_direction[Z]]*f[ec][0][iX*stride_any_direction[X]
                                                 + iZ*stride_any_direction[Z]]);
                                }
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                      np->s[ec][iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z]]
                                         = max(-np->energy[ec][iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z]]*fac, 0.0);
                                      op->P[ec][0][iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z]]
                                         = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                            + iY + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                               + iY + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                                  + iY + iZ*stride_any_direction[Z]]*f[ec][0][iX*stride_any_direction[X]
                                                     + iY + iZ*stride_any_direction[Z]]);
                                    }
                                  }
                                }
                              } else {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                      np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                         + iZ*stride_any_direction[Z]] = max(-np->energy[ec][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*fac,
                                               0.0);
                                      op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                         + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                               + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                                  + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                                     + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*f[ec][0][iX*stride_any_direction[X]
                                                        + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]);
                                    }
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                np->s[ec][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                   = max(-np->energy[ec][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]*fac,
                                      0.0);
                                op->P[ec][0][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                   = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                      + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                         + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                            + iZ*stride_any_direction[Z]]*f[ec][0][iX*stride_any_direction[X]
                                               + iZ*stride_any_direction[Z]]);
                              }
                            }
                          }
                        }
                      }
                    } else {
                      for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                        np->s[ec][iZ*stride_any_direction[Z]] = max(-np->energy[ec][iZ*stride_any_direction[Z]]*fac,
                           0.0);
                        op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                           + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                      }
                    }
                  }
                }
              }
            } else {
              if (stride_any_direction[R]) {
                if (num_any_direction[R]==1) {
                  np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                  op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                     + np->s[ec][0]*f[ec][0][0]);
                } else {
                  if (stride_any_direction[R]==1) {
                    for (int iR=0; iR<num_any_direction[R]; iR += 1) {
                      np->s[ec][iR] = max(-np->energy[ec][iR]*fac, 0.0);
                      op->P[ec][0][iR] = funinv*((2-om*om)*np->P[ec][0][iR] + (0.5*g-1)*op->P[ec][0][iR] + np->s[ec][iR]*f[ec][0][iR]);
                    }
                  } else {
                    for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                      np->s[ec][iR*stride_any_direction[R]] = max(-np->energy[ec][iR*stride_any_direction[R]]*fac,
                         0.0);
                      op->P[ec][0][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][0][iR*stride_any_direction[R]]
                         + (0.5*g-1)*op->P[ec][0][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][0][iR*stride_any_direction[R]]);
                    }
                  }
                }
              } else {
                if (stride_any_direction[X]) {
                  if (num_any_direction[X]==1) {
                    if (stride_any_direction[X]==1) {
                      if (stride_any_direction[Y]) {
                        if (num_any_direction[Y]==1) {
                          np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                          op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                             + np->s[ec][0]*f[ec][0][0]);
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            np->s[ec][iY*stride_any_direction[Y]] = max(-np->energy[ec][iY*stride_any_direction[Y]]*fac,
                               0.0);
                            op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                               + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                          }
                        }
                      } else {
                        np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                        op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                           + np->s[ec][0]*f[ec][0][0]);
                      }
                    } else {
                      if (stride_any_direction[Y]) {
                        if (num_any_direction[Y]==1) {
                          np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                          op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                             + np->s[ec][0]*f[ec][0][0]);
                        } else {
                          if (stride_any_direction[Y]==1) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              np->s[ec][iY] = max(-np->energy[ec][iY]*fac, 0.0);
                              op->P[ec][0][iY] = funinv*((2-om*om)*np->P[ec][0][iY] + (0.5*g-1)*op->P[ec][0][iY] + np->s[ec][iY]*f[ec][0][iY]);
                            }
                          } else {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              np->s[ec][iY*stride_any_direction[Y]] = max(-np->energy[ec][iY*stride_any_direction[Y]]*fac,
                                 0.0);
                              op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                                 + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                            }
                          }
                        }
                      } else {
                        np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                        op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                           + np->s[ec][0]*f[ec][0][0]);
                      }
                    }
                  } else {
                    if (stride_any_direction[X]==1) {
                      if (stride_any_direction[Y]) {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                            np->s[ec][iX] = max(-np->energy[ec][iX]*fac, 0.0);
                            op->P[ec][0][iX] = funinv*((2-om*om)*np->P[ec][0][iX] + (0.5*g-1)*op->P[ec][0][iX] + np->s[ec][iX]*f[ec][0][iX]);
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              np->s[ec][iX + iY*stride_any_direction[Y]] = max(-np->energy[ec][iX
                                 + iY*stride_any_direction[Y]]*fac, 0.0);
                              op->P[ec][0][iX + iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iX
                                 + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX
                                    + iY*stride_any_direction[Y]] + np->s[ec][iX + iY*stride_any_direction[Y]]*f[ec][0][iX
                                       + iY*stride_any_direction[Y]]);
                            }
                          }
                        }
                      } else {
                        for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                          np->s[ec][iX] = max(-np->energy[ec][iX]*fac, 0.0);
                          op->P[ec][0][iX] = funinv*((2-om*om)*np->P[ec][0][iX] + (0.5*g-1)*op->P[ec][0][iX] + np->s[ec][iX]*f[ec][0][iX]);
                        }
                      }
                    } else {
                      if (stride_any_direction[Y]) {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            np->s[ec][iX*stride_any_direction[X]] = max(-np->energy[ec][iX*stride_any_direction[X]]*fac,
                               0.0);
                            op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                               + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                          }
                        } else {
                          if (stride_any_direction[Y]==1) {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                np->s[ec][iX*stride_any_direction[X] + iY] = max(-np->energy[ec][iX*stride_any_direction[X]
                                   + iY]*fac, 0.0);
                                op->P[ec][0][iX*stride_any_direction[X] + iY] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                   + iY] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                      + iY] + np->s[ec][iX*stride_any_direction[X] + iY]*f[ec][0][iX*stride_any_direction[X]
                                         + iY]);
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = max(-np->energy[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]*fac,
                                      0.0);
                                op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                      + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                         + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y]]*f[ec][0][iX*stride_any_direction[X]
                                               + iY*stride_any_direction[Y]]);
                              }
                            }
                          }
                        }
                      } else {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          np->s[ec][iX*stride_any_direction[X]] = max(-np->energy[ec][iX*stride_any_direction[X]]*fac,
                             0.0);
                          op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                             + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                        }
                      }
                    }
                  }
                } else {
                  np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                  op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                     + np->s[ec][0]*f[ec][0][0]);
                }
              }
            }
          } else {
            for (int i=0;i<ntot;i++) {
              np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
            }
            if (stride_any_direction[Z]) {
              if (num_any_direction[Z]==1) {
                if (stride_any_direction[Z]==1) {
                  if (stride_any_direction[R]) {
                    if (num_any_direction[R]==1) {
                      np->s[ec][0] = np->energy[ec][0]*fac;
                      op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                         + np->s[ec][0]*f[ec][0][0]);
                    } else {
                      for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                        np->s[ec][iR*stride_any_direction[R]] = np->energy[ec][iR*stride_any_direction[R]]*fac;
                        op->P[ec][0][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][0][iR*stride_any_direction[R]]
                           + (0.5*g-1)*op->P[ec][0][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][0][iR*stride_any_direction[R]]);
                      }
                    }
                  } else {
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[Y]) {
                          if (num_any_direction[Y]==1) {
                            np->s[ec][0] = np->energy[ec][0]*fac;
                            op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                               + np->s[ec][0]*f[ec][0][0]);
                          } else {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              np->s[ec][iY*stride_any_direction[Y]] = np->energy[ec][iY*stride_any_direction[Y]]*fac;
                              op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                                 + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                            }
                          }
                        } else {
                          np->s[ec][0] = np->energy[ec][0]*fac;
                          op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                             + np->s[ec][0]*f[ec][0][0]);
                        }
                      } else {
                        if (stride_any_direction[Y]) {
                          if (num_any_direction[Y]==1) {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              np->s[ec][iX*stride_any_direction[X]] = np->energy[ec][iX*stride_any_direction[X]]*fac;
                              op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                                 + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = np->energy[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]*fac;
                                op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                      + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                         + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y]]*f[ec][0][iX*stride_any_direction[X]
                                               + iY*stride_any_direction[Y]]);
                              }
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            np->s[ec][iX*stride_any_direction[X]] = np->energy[ec][iX*stride_any_direction[X]]*fac;
                            op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                               + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                          }
                        }
                      }
                    } else {
                      np->s[ec][0] = np->energy[ec][0]*fac;
                      op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                         + np->s[ec][0]*f[ec][0][0]);
                    }
                  }
                } else {
                  if (stride_any_direction[R]) {
                    if (num_any_direction[R]==1) {
                      np->s[ec][0] = np->energy[ec][0]*fac;
                      op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                         + np->s[ec][0]*f[ec][0][0]);
                    } else {
                      if (stride_any_direction[R]==1) {
                        for (int iR=0; iR<num_any_direction[R]; iR += 1) {
                          np->s[ec][iR] = np->energy[ec][iR]*fac;
                          op->P[ec][0][iR] = funinv*((2-om*om)*np->P[ec][0][iR] + (0.5*g-1)*op->P[ec][0][iR] + np->s[ec][iR]*f[ec][0][iR]);
                        }
                      } else {
                        for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                          np->s[ec][iR*stride_any_direction[R]] = np->energy[ec][iR*stride_any_direction[R]]*fac;
                          op->P[ec][0][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][0][iR*stride_any_direction[R]]
                             + (0.5*g-1)*op->P[ec][0][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][0][iR*stride_any_direction[R]]);
                        }
                      }
                    }
                  } else {
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              np->s[ec][0] = np->energy[ec][0]*fac;
                              op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                                 + np->s[ec][0]*f[ec][0][0]);
                            } else {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                np->s[ec][iY*stride_any_direction[Y]] = np->energy[ec][iY*stride_any_direction[Y]]*fac;
                                op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                                   + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                              }
                            }
                          } else {
                            np->s[ec][0] = np->energy[ec][0]*fac;
                            op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                               + np->s[ec][0]*f[ec][0][0]);
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              np->s[ec][0] = np->energy[ec][0]*fac;
                              op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                                 + np->s[ec][0]*f[ec][0][0]);
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                  np->s[ec][iY] = np->energy[ec][iY]*fac;
                                  op->P[ec][0][iY] = funinv*((2-om*om)*np->P[ec][0][iY] + (0.5*g-1)*op->P[ec][0][iY] + np->s[ec][iY]*f[ec][0][iY]);
                                }
                              } else {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  np->s[ec][iY*stride_any_direction[Y]] = np->energy[ec][iY*stride_any_direction[Y]]*fac;
                                  op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                                     + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                                }
                              }
                            }
                          } else {
                            np->s[ec][0] = np->energy[ec][0]*fac;
                            op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                               + np->s[ec][0]*f[ec][0][0]);
                          }
                        }
                      } else {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                np->s[ec][iX] = np->energy[ec][iX]*fac;
                                op->P[ec][0][iX] = funinv*((2-om*om)*np->P[ec][0][iX] + (0.5*g-1)*op->P[ec][0][iX] + np->s[ec][iX]*f[ec][0][iX]);
                              }
                            } else {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  np->s[ec][iX + iY*stride_any_direction[Y]] = np->energy[ec][iX
                                     + iY*stride_any_direction[Y]]*fac;
                                  op->P[ec][0][iX + iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iX
                                     + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX
                                        + iY*stride_any_direction[Y]] + np->s[ec][iX + iY*stride_any_direction[Y]]*f[ec][0][iX
                                           + iY*stride_any_direction[Y]]);
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                              np->s[ec][iX] = np->energy[ec][iX]*fac;
                              op->P[ec][0][iX] = funinv*((2-om*om)*np->P[ec][0][iX] + (0.5*g-1)*op->P[ec][0][iX] + np->s[ec][iX]*f[ec][0][iX]);
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                np->s[ec][iX*stride_any_direction[X]] = np->energy[ec][iX*stride_any_direction[X]]*fac;
                                op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                                   + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                    np->s[ec][iX*stride_any_direction[X] + iY] = np->energy[ec][iX*stride_any_direction[X]
                                       + iY]*fac;
                                    op->P[ec][0][iX*stride_any_direction[X] + iY] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                       + iY] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                          + iY] + np->s[ec][iX*stride_any_direction[X] + iY]*f[ec][0][iX*stride_any_direction[X]
                                             + iY]);
                                  }
                                }
                              } else {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                    np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                       = np->energy[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]*fac;
                                    op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                       = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                          + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                             + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                                + iY*stride_any_direction[Y]]*f[ec][0][iX*stride_any_direction[X]
                                                   + iY*stride_any_direction[Y]]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              np->s[ec][iX*stride_any_direction[X]] = np->energy[ec][iX*stride_any_direction[X]]*fac;
                              op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                                 + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                            }
                          }
                        }
                      }
                    } else {
                      np->s[ec][0] = np->energy[ec][0]*fac;
                      op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                         + np->s[ec][0]*f[ec][0][0]);
                    }
                  }
                }
              } else {
                if (stride_any_direction[Z]==1) {
                  if (stride_any_direction[R]) {
                    if (num_any_direction[R]==1) {
                      for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                        np->s[ec][iZ] = np->energy[ec][iZ]*fac;
                        op->P[ec][0][iZ] = funinv*((2-om*om)*np->P[ec][0][iZ] + (0.5*g-1)*op->P[ec][0][iZ] + np->s[ec][iZ]*f[ec][0][iZ]);
                      }
                    } else {
                      for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                        for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                          np->s[ec][iZ + iR*stride_any_direction[R]] = np->energy[ec][iZ
                             + iR*stride_any_direction[R]]*fac;
                          op->P[ec][0][iZ + iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][0][iZ
                             + iR*stride_any_direction[R]] + (0.5*g-1)*op->P[ec][0][iZ
                                + iR*stride_any_direction[R]] + np->s[ec][iZ + iR*stride_any_direction[R]]*f[ec][0][iZ
                                   + iR*stride_any_direction[R]]);
                        }
                      }
                    }
                  } else {
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[Y]) {
                          if (num_any_direction[Y]==1) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              np->s[ec][iZ] = np->energy[ec][iZ]*fac;
                              op->P[ec][0][iZ] = funinv*((2-om*om)*np->P[ec][0][iZ] + (0.5*g-1)*op->P[ec][0][iZ] + np->s[ec][iZ]*f[ec][0][iZ]);
                            }
                          } else {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                np->s[ec][iY*stride_any_direction[Y] + iZ] = np->energy[ec][iY*stride_any_direction[Y]
                                   + iZ]*fac;
                                op->P[ec][0][iY*stride_any_direction[Y] + iZ] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]
                                   + iZ] + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]
                                      + iZ] + np->s[ec][iY*stride_any_direction[Y] + iZ]*f[ec][0][iY*stride_any_direction[Y]
                                         + iZ]);
                              }
                            }
                          }
                        } else {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                            np->s[ec][iZ] = np->energy[ec][iZ]*fac;
                            op->P[ec][0][iZ] = funinv*((2-om*om)*np->P[ec][0][iZ] + (0.5*g-1)*op->P[ec][0][iZ] + np->s[ec][iZ]*f[ec][0][iZ]);
                          }
                        }
                      } else {
                        if (stride_any_direction[Y]) {
                          if (num_any_direction[Y]==1) {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                np->s[ec][iX*stride_any_direction[X] + iZ] = np->energy[ec][iX*stride_any_direction[X]
                                   + iZ]*fac;
                                op->P[ec][0][iX*stride_any_direction[X] + iZ] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                   + iZ] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                      + iZ] + np->s[ec][iX*stride_any_direction[X] + iZ]*f[ec][0][iX*stride_any_direction[X]
                                         + iZ]);
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                  np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                     + iZ] = np->energy[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                        + iZ]*fac;
                                  op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                     + iZ] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                        + iY*stride_any_direction[Y] + iZ] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                           + iY*stride_any_direction[Y] + iZ] + np->s[ec][iX*stride_any_direction[X]
                                              + iY*stride_any_direction[Y] + iZ]*f[ec][0][iX*stride_any_direction[X]
                                                 + iY*stride_any_direction[Y] + iZ]);
                                }
                              }
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              np->s[ec][iX*stride_any_direction[X] + iZ] = np->energy[ec][iX*stride_any_direction[X]
                                 + iZ]*fac;
                              op->P[ec][0][iX*stride_any_direction[X] + iZ] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                 + iZ] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                    + iZ] + np->s[ec][iX*stride_any_direction[X] + iZ]*f[ec][0][iX*stride_any_direction[X]
                                       + iZ]);
                            }
                          }
                        }
                      }
                    } else {
                      for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                        np->s[ec][iZ] = np->energy[ec][iZ]*fac;
                        op->P[ec][0][iZ] = funinv*((2-om*om)*np->P[ec][0][iZ] + (0.5*g-1)*op->P[ec][0][iZ] + np->s[ec][iZ]*f[ec][0][iZ]);
                      }
                    }
                  }
                } else {
                  if (stride_any_direction[R]) {
                    if (num_any_direction[R]==1) {
                      for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                        np->s[ec][iZ*stride_any_direction[Z]] = np->energy[ec][iZ*stride_any_direction[Z]]*fac;
                        op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                           + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                      }
                    } else {
                      if (stride_any_direction[R]==1) {
                        for (int iR=0; iR<num_any_direction[R]; iR += 1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                            np->s[ec][iZ*stride_any_direction[Z] + iR] = np->energy[ec][iZ*stride_any_direction[Z]
                               + iR]*fac;
                            op->P[ec][0][iZ*stride_any_direction[Z] + iR] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]
                               + iR] + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]
                                  + iR] + np->s[ec][iZ*stride_any_direction[Z] + iR]*f[ec][0][iZ*stride_any_direction[Z]
                                     + iR]);
                          }
                        }
                      } else {
                        for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                            np->s[ec][iZ*stride_any_direction[Z] + iR*stride_any_direction[R]]
                               = np->energy[ec][iZ*stride_any_direction[Z] + iR*stride_any_direction[R]]*fac;
                            op->P[ec][0][iZ*stride_any_direction[Z] + iR*stride_any_direction[R]]
                               = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]
                                  + iR*stride_any_direction[R]] + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]
                                     + iR*stride_any_direction[R]] + np->s[ec][iZ*stride_any_direction[Z]
                                        + iR*stride_any_direction[R]]*f[ec][0][iZ*stride_any_direction[Z]
                                           + iR*stride_any_direction[R]]);
                          }
                        }
                      }
                    }
                  } else {
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                np->s[ec][iZ*stride_any_direction[Z]] = np->energy[ec][iZ*stride_any_direction[Z]]*fac;
                                op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                                   + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                              }
                            } else {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                  np->s[ec][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                     = np->energy[ec][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*fac;
                                  op->P[ec][0][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                     = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]
                                        + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]
                                           + iZ*stride_any_direction[Z]] + np->s[ec][iY*stride_any_direction[Y]
                                              + iZ*stride_any_direction[Z]]*f[ec][0][iY*stride_any_direction[Y]
                                                 + iZ*stride_any_direction[Z]]);
                                }
                              }
                            }
                          } else {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                              np->s[ec][iZ*stride_any_direction[Z]] = np->energy[ec][iZ*stride_any_direction[Z]]*fac;
                              op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                                 + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                np->s[ec][iZ*stride_any_direction[Z]] = np->energy[ec][iZ*stride_any_direction[Z]]*fac;
                                op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                                   + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                    np->s[ec][iY + iZ*stride_any_direction[Z]] = np->energy[ec][iY
                                       + iZ*stride_any_direction[Z]]*fac;
                                    op->P[ec][0][iY + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iY
                                       + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iY
                                          + iZ*stride_any_direction[Z]] + np->s[ec][iY + iZ*stride_any_direction[Z]]*f[ec][0][iY
                                             + iZ*stride_any_direction[Z]]);
                                  }
                                }
                              } else {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                    np->s[ec][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = np->energy[ec][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*fac;
                                    op->P[ec][0][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]
                                          + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]
                                             + iZ*stride_any_direction[Z]] + np->s[ec][iY*stride_any_direction[Y]
                                                + iZ*stride_any_direction[Z]]*f[ec][0][iY*stride_any_direction[Y]
                                                   + iZ*stride_any_direction[Z]]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                              np->s[ec][iZ*stride_any_direction[Z]] = np->energy[ec][iZ*stride_any_direction[Z]]*fac;
                              op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                                 + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                            }
                          }
                        }
                      } else {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                  np->s[ec][iX + iZ*stride_any_direction[Z]] = np->energy[ec][iX
                                     + iZ*stride_any_direction[Z]]*fac;
                                  op->P[ec][0][iX + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iX
                                     + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX
                                        + iZ*stride_any_direction[Z]] + np->s[ec][iX + iZ*stride_any_direction[Z]]*f[ec][0][iX
                                           + iZ*stride_any_direction[Z]]);
                                }
                              }
                            } else {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                    np->s[ec][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = np->energy[ec][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*fac;
                                    op->P[ec][0][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = funinv*((2-om*om)*np->P[ec][0][iX + iY*stride_any_direction[Y]
                                          + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX
                                             + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                                + np->s[ec][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*f[ec][0][iX
                                                   + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                np->s[ec][iX + iZ*stride_any_direction[Z]] = np->energy[ec][iX
                                   + iZ*stride_any_direction[Z]]*fac;
                                op->P[ec][0][iX + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iX
                                   + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX
                                      + iZ*stride_any_direction[Z]] + np->s[ec][iX + iZ*stride_any_direction[Z]]*f[ec][0][iX
                                         + iZ*stride_any_direction[Z]]);
                              }
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                  np->s[ec][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                     = np->energy[ec][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]*fac;
                                  op->P[ec][0][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                     = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                        + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                           + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                              + iZ*stride_any_direction[Z]]*f[ec][0][iX*stride_any_direction[X]
                                                 + iZ*stride_any_direction[Z]]);
                                }
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                      np->s[ec][iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z]]
                                         = np->energy[ec][iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z]]*fac;
                                      op->P[ec][0][iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z]]
                                         = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                            + iY + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                               + iY + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                                  + iY + iZ*stride_any_direction[Z]]*f[ec][0][iX*stride_any_direction[X]
                                                     + iY + iZ*stride_any_direction[Z]]);
                                    }
                                  }
                                }
                              } else {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                      np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                         + iZ*stride_any_direction[Z]] = np->energy[ec][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*fac;
                                      op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                         + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                               + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                                  + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                                     + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*f[ec][0][iX*stride_any_direction[X]
                                                        + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]);
                                    }
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                np->s[ec][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                   = np->energy[ec][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]*fac;
                                op->P[ec][0][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                   = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                      + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                         + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                            + iZ*stride_any_direction[Z]]*f[ec][0][iX*stride_any_direction[X]
                                               + iZ*stride_any_direction[Z]]);
                              }
                            }
                          }
                        }
                      }
                    } else {
                      for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                        np->s[ec][iZ*stride_any_direction[Z]] = np->energy[ec][iZ*stride_any_direction[Z]]*fac;
                        op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                           + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                      }
                    }
                  }
                }
              }
            } else {
              if (stride_any_direction[R]) {
                if (num_any_direction[R]==1) {
                  np->s[ec][0] = np->energy[ec][0]*fac;
                  op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                     + np->s[ec][0]*f[ec][0][0]);
                } else {
                  if (stride_any_direction[R]==1) {
                    for (int iR=0; iR<num_any_direction[R]; iR += 1) {
                      np->s[ec][iR] = np->energy[ec][iR]*fac;
                      op->P[ec][0][iR] = funinv*((2-om*om)*np->P[ec][0][iR] + (0.5*g-1)*op->P[ec][0][iR] + np->s[ec][iR]*f[ec][0][iR]);
                    }
                  } else {
                    for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                      np->s[ec][iR*stride_any_direction[R]] = np->energy[ec][iR*stride_any_direction[R]]*fac;
                      op->P[ec][0][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][0][iR*stride_any_direction[R]]
                         + (0.5*g-1)*op->P[ec][0][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][0][iR*stride_any_direction[R]]);
                    }
                  }
                }
              } else {
                if (stride_any_direction[X]) {
                  if (num_any_direction[X]==1) {
                    if (stride_any_direction[X]==1) {
                      if (stride_any_direction[Y]) {
                        if (num_any_direction[Y]==1) {
                          np->s[ec][0] = np->energy[ec][0]*fac;
                          op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                             + np->s[ec][0]*f[ec][0][0]);
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            np->s[ec][iY*stride_any_direction[Y]] = np->energy[ec][iY*stride_any_direction[Y]]*fac;
                            op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                               + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                          }
                        }
                      } else {
                        np->s[ec][0] = np->energy[ec][0]*fac;
                        op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                           + np->s[ec][0]*f[ec][0][0]);
                      }
                    } else {
                      if (stride_any_direction[Y]) {
                        if (num_any_direction[Y]==1) {
                          np->s[ec][0] = np->energy[ec][0]*fac;
                          op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                             + np->s[ec][0]*f[ec][0][0]);
                        } else {
                          if (stride_any_direction[Y]==1) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              np->s[ec][iY] = np->energy[ec][iY]*fac;
                              op->P[ec][0][iY] = funinv*((2-om*om)*np->P[ec][0][iY] + (0.5*g-1)*op->P[ec][0][iY] + np->s[ec][iY]*f[ec][0][iY]);
                            }
                          } else {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              np->s[ec][iY*stride_any_direction[Y]] = np->energy[ec][iY*stride_any_direction[Y]]*fac;
                              op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                                 + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                            }
                          }
                        }
                      } else {
                        np->s[ec][0] = np->energy[ec][0]*fac;
                        op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                           + np->s[ec][0]*f[ec][0][0]);
                      }
                    }
                  } else {
                    if (stride_any_direction[X]==1) {
                      if (stride_any_direction[Y]) {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                            np->s[ec][iX] = np->energy[ec][iX]*fac;
                            op->P[ec][0][iX] = funinv*((2-om*om)*np->P[ec][0][iX] + (0.5*g-1)*op->P[ec][0][iX] + np->s[ec][iX]*f[ec][0][iX]);
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              np->s[ec][iX + iY*stride_any_direction[Y]] = np->energy[ec][iX
                                 + iY*stride_any_direction[Y]]*fac;
                              op->P[ec][0][iX + iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iX
                                 + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX
                                    + iY*stride_any_direction[Y]] + np->s[ec][iX + iY*stride_any_direction[Y]]*f[ec][0][iX
                                       + iY*stride_any_direction[Y]]);
                            }
                          }
                        }
                      } else {
                        for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                          np->s[ec][iX] = np->energy[ec][iX]*fac;
                          op->P[ec][0][iX] = funinv*((2-om*om)*np->P[ec][0][iX] + (0.5*g-1)*op->P[ec][0][iX] + np->s[ec][iX]*f[ec][0][iX]);
                        }
                      }
                    } else {
                      if (stride_any_direction[Y]) {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            np->s[ec][iX*stride_any_direction[X]] = np->energy[ec][iX*stride_any_direction[X]]*fac;
                            op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                               + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                          }
                        } else {
                          if (stride_any_direction[Y]==1) {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                np->s[ec][iX*stride_any_direction[X] + iY] = np->energy[ec][iX*stride_any_direction[X]
                                   + iY]*fac;
                                op->P[ec][0][iX*stride_any_direction[X] + iY] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                   + iY] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                      + iY] + np->s[ec][iX*stride_any_direction[X] + iY]*f[ec][0][iX*stride_any_direction[X]
                                         + iY]);
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = np->energy[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]*fac;
                                op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                      + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                         + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y]]*f[ec][0][iX*stride_any_direction[X]
                                               + iY*stride_any_direction[Y]]);
                              }
                            }
                          }
                        }
                      } else {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          np->s[ec][iX*stride_any_direction[X]] = np->energy[ec][iX*stride_any_direction[X]]*fac;
                          op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                             + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                        }
                      }
                    }
                  }
                } else {
                  np->s[ec][0] = np->energy[ec][0]*fac;
                  op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                     + np->s[ec][0]*f[ec][0][0]);
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
            if (stride_any_direction[Z]) {
              if (num_any_direction[Z]==1) {
                if (stride_any_direction[Z]==1) {
                  if (stride_any_direction[R]) {
                    if (num_any_direction[R]==1) {
                      np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                      op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                         + np->s[ec][0]*f[ec][0][0]);
                      op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                         + np->s[ec][0]*f[ec][1][0]);
                    } else {
                      for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                        np->s[ec][iR*stride_any_direction[R]] = max(-np->energy[ec][iR*stride_any_direction[R]]*fac,
                           0.0);
                        op->P[ec][0][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][0][iR*stride_any_direction[R]]
                           + (0.5*g-1)*op->P[ec][0][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][0][iR*stride_any_direction[R]]);
                        op->P[ec][1][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][1][iR*stride_any_direction[R]]
                           + (0.5*g-1)*op->P[ec][1][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][1][iR*stride_any_direction[R]]);
                      }
                    }
                  } else {
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[Y]) {
                          if (num_any_direction[Y]==1) {
                            np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                            op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                               + np->s[ec][0]*f[ec][0][0]);
                            op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                               + np->s[ec][0]*f[ec][1][0]);
                          } else {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              np->s[ec][iY*stride_any_direction[Y]] = max(-np->energy[ec][iY*stride_any_direction[Y]]*fac,
                                 0.0);
                              op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                                 + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                              op->P[ec][1][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][1][iY*stride_any_direction[Y]]
                                 + (0.5*g-1)*op->P[ec][1][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][1][iY*stride_any_direction[Y]]);
                            }
                          }
                        } else {
                          np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                          op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                             + np->s[ec][0]*f[ec][0][0]);
                          op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                             + np->s[ec][0]*f[ec][1][0]);
                        }
                      } else {
                        if (stride_any_direction[Y]) {
                          if (num_any_direction[Y]==1) {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              np->s[ec][iX*stride_any_direction[X]] = max(-np->energy[ec][iX*stride_any_direction[X]]*fac,
                                 0.0);
                              op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                                 + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                              op->P[ec][1][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]]
                                 + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][1][iX*stride_any_direction[X]]);
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = max(-np->energy[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]*fac,
                                      0.0);
                                op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                      + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                         + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y]]*f[ec][0][iX*stride_any_direction[X]
                                               + iY*stride_any_direction[Y]]);
                                op->P[ec][1][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                      + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                         + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y]]*f[ec][1][iX*stride_any_direction[X]
                                               + iY*stride_any_direction[Y]]);
                              }
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            np->s[ec][iX*stride_any_direction[X]] = max(-np->energy[ec][iX*stride_any_direction[X]]*fac,
                               0.0);
                            op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                               + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                            op->P[ec][1][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]]
                               + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][1][iX*stride_any_direction[X]]);
                          }
                        }
                      }
                    } else {
                      np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                      op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                         + np->s[ec][0]*f[ec][0][0]);
                      op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                         + np->s[ec][0]*f[ec][1][0]);
                    }
                  }
                } else {
                  if (stride_any_direction[R]) {
                    if (num_any_direction[R]==1) {
                      np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                      op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                         + np->s[ec][0]*f[ec][0][0]);
                      op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                         + np->s[ec][0]*f[ec][1][0]);
                    } else {
                      if (stride_any_direction[R]==1) {
                        for (int iR=0; iR<num_any_direction[R]; iR += 1) {
                          np->s[ec][iR] = max(-np->energy[ec][iR]*fac, 0.0);
                          op->P[ec][0][iR] = funinv*((2-om*om)*np->P[ec][0][iR] + (0.5*g-1)*op->P[ec][0][iR] + np->s[ec][iR]*f[ec][0][iR]);
                          op->P[ec][1][iR] = funinv*((2-om*om)*np->P[ec][1][iR] + (0.5*g-1)*op->P[ec][1][iR] + np->s[ec][iR]*f[ec][1][iR]);
                        }
                      } else {
                        for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                          np->s[ec][iR*stride_any_direction[R]] = max(-np->energy[ec][iR*stride_any_direction[R]]*fac,
                             0.0);
                          op->P[ec][0][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][0][iR*stride_any_direction[R]]
                             + (0.5*g-1)*op->P[ec][0][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][0][iR*stride_any_direction[R]]);
                          op->P[ec][1][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][1][iR*stride_any_direction[R]]
                             + (0.5*g-1)*op->P[ec][1][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][1][iR*stride_any_direction[R]]);
                        }
                      }
                    }
                  } else {
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                              op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                                 + np->s[ec][0]*f[ec][0][0]);
                              op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                                 + np->s[ec][0]*f[ec][1][0]);
                            } else {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                np->s[ec][iY*stride_any_direction[Y]] = max(-np->energy[ec][iY*stride_any_direction[Y]]*fac,
                                   0.0);
                                op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                                   + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                                op->P[ec][1][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][1][iY*stride_any_direction[Y]]
                                   + (0.5*g-1)*op->P[ec][1][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][1][iY*stride_any_direction[Y]]);
                              }
                            }
                          } else {
                            np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                            op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                               + np->s[ec][0]*f[ec][0][0]);
                            op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                               + np->s[ec][0]*f[ec][1][0]);
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                              op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                                 + np->s[ec][0]*f[ec][0][0]);
                              op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                                 + np->s[ec][0]*f[ec][1][0]);
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                  np->s[ec][iY] = max(-np->energy[ec][iY]*fac, 0.0);
                                  op->P[ec][0][iY] = funinv*((2-om*om)*np->P[ec][0][iY] + (0.5*g-1)*op->P[ec][0][iY] + np->s[ec][iY]*f[ec][0][iY]);
                                  op->P[ec][1][iY] = funinv*((2-om*om)*np->P[ec][1][iY] + (0.5*g-1)*op->P[ec][1][iY] + np->s[ec][iY]*f[ec][1][iY]);
                                }
                              } else {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  np->s[ec][iY*stride_any_direction[Y]] = max(-np->energy[ec][iY*stride_any_direction[Y]]*fac,
                                     0.0);
                                  op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                                     + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                                  op->P[ec][1][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][1][iY*stride_any_direction[Y]]
                                     + (0.5*g-1)*op->P[ec][1][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][1][iY*stride_any_direction[Y]]);
                                }
                              }
                            }
                          } else {
                            np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                            op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                               + np->s[ec][0]*f[ec][0][0]);
                            op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                               + np->s[ec][0]*f[ec][1][0]);
                          }
                        }
                      } else {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                np->s[ec][iX] = max(-np->energy[ec][iX]*fac, 0.0);
                                op->P[ec][0][iX] = funinv*((2-om*om)*np->P[ec][0][iX] + (0.5*g-1)*op->P[ec][0][iX] + np->s[ec][iX]*f[ec][0][iX]);
                                op->P[ec][1][iX] = funinv*((2-om*om)*np->P[ec][1][iX] + (0.5*g-1)*op->P[ec][1][iX] + np->s[ec][iX]*f[ec][1][iX]);
                              }
                            } else {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  np->s[ec][iX + iY*stride_any_direction[Y]] = max(-np->energy[ec][iX
                                     + iY*stride_any_direction[Y]]*fac, 0.0);
                                  op->P[ec][0][iX + iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iX
                                     + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX
                                        + iY*stride_any_direction[Y]] + np->s[ec][iX + iY*stride_any_direction[Y]]*f[ec][0][iX
                                           + iY*stride_any_direction[Y]]);
                                  op->P[ec][1][iX + iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][1][iX
                                     + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][1][iX
                                        + iY*stride_any_direction[Y]] + np->s[ec][iX + iY*stride_any_direction[Y]]*f[ec][1][iX
                                           + iY*stride_any_direction[Y]]);
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                              np->s[ec][iX] = max(-np->energy[ec][iX]*fac, 0.0);
                              op->P[ec][0][iX] = funinv*((2-om*om)*np->P[ec][0][iX] + (0.5*g-1)*op->P[ec][0][iX] + np->s[ec][iX]*f[ec][0][iX]);
                              op->P[ec][1][iX] = funinv*((2-om*om)*np->P[ec][1][iX] + (0.5*g-1)*op->P[ec][1][iX] + np->s[ec][iX]*f[ec][1][iX]);
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                np->s[ec][iX*stride_any_direction[X]] = max(-np->energy[ec][iX*stride_any_direction[X]]*fac,
                                   0.0);
                                op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                                   + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                                op->P[ec][1][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]]
                                   + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][1][iX*stride_any_direction[X]]);
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                    np->s[ec][iX*stride_any_direction[X] + iY] = max(-np->energy[ec][iX*stride_any_direction[X]
                                       + iY]*fac, 0.0);
                                    op->P[ec][0][iX*stride_any_direction[X] + iY] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                       + iY] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                          + iY] + np->s[ec][iX*stride_any_direction[X] + iY]*f[ec][0][iX*stride_any_direction[X]
                                             + iY]);
                                    op->P[ec][1][iX*stride_any_direction[X] + iY] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                       + iY] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                          + iY] + np->s[ec][iX*stride_any_direction[X] + iY]*f[ec][1][iX*stride_any_direction[X]
                                             + iY]);
                                  }
                                }
                              } else {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                    np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                       = max(-np->energy[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]*fac,
                                          0.0);
                                    op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                       = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                          + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                             + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                                + iY*stride_any_direction[Y]]*f[ec][0][iX*stride_any_direction[X]
                                                   + iY*stride_any_direction[Y]]);
                                    op->P[ec][1][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                       = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                          + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                             + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                                + iY*stride_any_direction[Y]]*f[ec][1][iX*stride_any_direction[X]
                                                   + iY*stride_any_direction[Y]]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              np->s[ec][iX*stride_any_direction[X]] = max(-np->energy[ec][iX*stride_any_direction[X]]*fac,
                                 0.0);
                              op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                                 + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                              op->P[ec][1][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]]
                                 + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][1][iX*stride_any_direction[X]]);
                            }
                          }
                        }
                      }
                    } else {
                      np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                      op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                         + np->s[ec][0]*f[ec][0][0]);
                      op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                         + np->s[ec][0]*f[ec][1][0]);
                    }
                  }
                }
              } else {
                if (stride_any_direction[Z]==1) {
                  if (stride_any_direction[R]) {
                    if (num_any_direction[R]==1) {
                      for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                        np->s[ec][iZ] = max(-np->energy[ec][iZ]*fac, 0.0);
                        op->P[ec][0][iZ] = funinv*((2-om*om)*np->P[ec][0][iZ] + (0.5*g-1)*op->P[ec][0][iZ] + np->s[ec][iZ]*f[ec][0][iZ]);
                        op->P[ec][1][iZ] = funinv*((2-om*om)*np->P[ec][1][iZ] + (0.5*g-1)*op->P[ec][1][iZ] + np->s[ec][iZ]*f[ec][1][iZ]);
                      }
                    } else {
                      for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                        for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                          np->s[ec][iZ + iR*stride_any_direction[R]] = max(-np->energy[ec][iZ
                             + iR*stride_any_direction[R]]*fac, 0.0);
                          op->P[ec][0][iZ + iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][0][iZ
                             + iR*stride_any_direction[R]] + (0.5*g-1)*op->P[ec][0][iZ
                                + iR*stride_any_direction[R]] + np->s[ec][iZ + iR*stride_any_direction[R]]*f[ec][0][iZ
                                   + iR*stride_any_direction[R]]);
                          op->P[ec][1][iZ + iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][1][iZ
                             + iR*stride_any_direction[R]] + (0.5*g-1)*op->P[ec][1][iZ
                                + iR*stride_any_direction[R]] + np->s[ec][iZ + iR*stride_any_direction[R]]*f[ec][1][iZ
                                   + iR*stride_any_direction[R]]);
                        }
                      }
                    }
                  } else {
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[Y]) {
                          if (num_any_direction[Y]==1) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              np->s[ec][iZ] = max(-np->energy[ec][iZ]*fac, 0.0);
                              op->P[ec][0][iZ] = funinv*((2-om*om)*np->P[ec][0][iZ] + (0.5*g-1)*op->P[ec][0][iZ] + np->s[ec][iZ]*f[ec][0][iZ]);
                              op->P[ec][1][iZ] = funinv*((2-om*om)*np->P[ec][1][iZ] + (0.5*g-1)*op->P[ec][1][iZ] + np->s[ec][iZ]*f[ec][1][iZ]);
                            }
                          } else {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                np->s[ec][iY*stride_any_direction[Y] + iZ] = max(-np->energy[ec][iY*stride_any_direction[Y]
                                   + iZ]*fac, 0.0);
                                op->P[ec][0][iY*stride_any_direction[Y] + iZ] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]
                                   + iZ] + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]
                                      + iZ] + np->s[ec][iY*stride_any_direction[Y] + iZ]*f[ec][0][iY*stride_any_direction[Y]
                                         + iZ]);
                                op->P[ec][1][iY*stride_any_direction[Y] + iZ] = funinv*((2-om*om)*np->P[ec][1][iY*stride_any_direction[Y]
                                   + iZ] + (0.5*g-1)*op->P[ec][1][iY*stride_any_direction[Y]
                                      + iZ] + np->s[ec][iY*stride_any_direction[Y] + iZ]*f[ec][1][iY*stride_any_direction[Y]
                                         + iZ]);
                              }
                            }
                          }
                        } else {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                            np->s[ec][iZ] = max(-np->energy[ec][iZ]*fac, 0.0);
                            op->P[ec][0][iZ] = funinv*((2-om*om)*np->P[ec][0][iZ] + (0.5*g-1)*op->P[ec][0][iZ] + np->s[ec][iZ]*f[ec][0][iZ]);
                            op->P[ec][1][iZ] = funinv*((2-om*om)*np->P[ec][1][iZ] + (0.5*g-1)*op->P[ec][1][iZ] + np->s[ec][iZ]*f[ec][1][iZ]);
                          }
                        }
                      } else {
                        if (stride_any_direction[Y]) {
                          if (num_any_direction[Y]==1) {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                np->s[ec][iX*stride_any_direction[X] + iZ] = max(-np->energy[ec][iX*stride_any_direction[X]
                                   + iZ]*fac, 0.0);
                                op->P[ec][0][iX*stride_any_direction[X] + iZ] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                   + iZ] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                      + iZ] + np->s[ec][iX*stride_any_direction[X] + iZ]*f[ec][0][iX*stride_any_direction[X]
                                         + iZ]);
                                op->P[ec][1][iX*stride_any_direction[X] + iZ] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                   + iZ] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                      + iZ] + np->s[ec][iX*stride_any_direction[X] + iZ]*f[ec][1][iX*stride_any_direction[X]
                                         + iZ]);
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                  np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                     + iZ] = max(-np->energy[ec][iX*stride_any_direction[X] +
                                        iY*stride_any_direction[Y] + iZ]*fac, 0.0);
                                  op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                     + iZ] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                        + iY*stride_any_direction[Y] + iZ] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                           + iY*stride_any_direction[Y] + iZ] + np->s[ec][iX*stride_any_direction[X]
                                              + iY*stride_any_direction[Y] + iZ]*f[ec][0][iX*stride_any_direction[X]
                                                 + iY*stride_any_direction[Y] + iZ]);
                                  op->P[ec][1][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                     + iZ] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                        + iY*stride_any_direction[Y] + iZ] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                           + iY*stride_any_direction[Y] + iZ] + np->s[ec][iX*stride_any_direction[X]
                                              + iY*stride_any_direction[Y] + iZ]*f[ec][1][iX*stride_any_direction[X]
                                                 + iY*stride_any_direction[Y] + iZ]);
                                }
                              }
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              np->s[ec][iX*stride_any_direction[X] + iZ] = max(-np->energy[ec][iX*stride_any_direction[X]
                                 + iZ]*fac, 0.0);
                              op->P[ec][0][iX*stride_any_direction[X] + iZ] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                 + iZ] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                    + iZ] + np->s[ec][iX*stride_any_direction[X] + iZ]*f[ec][0][iX*stride_any_direction[X]
                                       + iZ]);
                              op->P[ec][1][iX*stride_any_direction[X] + iZ] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                 + iZ] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                    + iZ] + np->s[ec][iX*stride_any_direction[X] + iZ]*f[ec][1][iX*stride_any_direction[X]
                                       + iZ]);
                            }
                          }
                        }
                      }
                    } else {
                      for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                        np->s[ec][iZ] = max(-np->energy[ec][iZ]*fac, 0.0);
                        op->P[ec][0][iZ] = funinv*((2-om*om)*np->P[ec][0][iZ] + (0.5*g-1)*op->P[ec][0][iZ] + np->s[ec][iZ]*f[ec][0][iZ]);
                        op->P[ec][1][iZ] = funinv*((2-om*om)*np->P[ec][1][iZ] + (0.5*g-1)*op->P[ec][1][iZ] + np->s[ec][iZ]*f[ec][1][iZ]);
                      }
                    }
                  }
                } else {
                  if (stride_any_direction[R]) {
                    if (num_any_direction[R]==1) {
                      for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                        np->s[ec][iZ*stride_any_direction[Z]] = max(-np->energy[ec][iZ*stride_any_direction[Z]]*fac,
                           0.0);
                        op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                           + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                        op->P[ec][1][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iZ*stride_any_direction[Z]]
                           + (0.5*g-1)*op->P[ec][1][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][1][iZ*stride_any_direction[Z]]);
                      }
                    } else {
                      if (stride_any_direction[R]==1) {
                        for (int iR=0; iR<num_any_direction[R]; iR += 1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                            np->s[ec][iZ*stride_any_direction[Z] + iR] = max(-np->energy[ec][iZ*stride_any_direction[Z]
                               + iR]*fac, 0.0);
                            op->P[ec][0][iZ*stride_any_direction[Z] + iR] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]
                               + iR] + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]
                                  + iR] + np->s[ec][iZ*stride_any_direction[Z] + iR]*f[ec][0][iZ*stride_any_direction[Z]
                                     + iR]);
                            op->P[ec][1][iZ*stride_any_direction[Z] + iR] = funinv*((2-om*om)*np->P[ec][1][iZ*stride_any_direction[Z]
                               + iR] + (0.5*g-1)*op->P[ec][1][iZ*stride_any_direction[Z]
                                  + iR] + np->s[ec][iZ*stride_any_direction[Z] + iR]*f[ec][1][iZ*stride_any_direction[Z]
                                     + iR]);
                          }
                        }
                      } else {
                        for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                            np->s[ec][iZ*stride_any_direction[Z] + iR*stride_any_direction[R]]
                               = max(-np->energy[ec][iZ*stride_any_direction[Z] + iR*stride_any_direction[R]]*fac,
                                  0.0);
                            op->P[ec][0][iZ*stride_any_direction[Z] + iR*stride_any_direction[R]]
                               = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]
                                  + iR*stride_any_direction[R]] + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]
                                     + iR*stride_any_direction[R]] + np->s[ec][iZ*stride_any_direction[Z]
                                        + iR*stride_any_direction[R]]*f[ec][0][iZ*stride_any_direction[Z]
                                           + iR*stride_any_direction[R]]);
                            op->P[ec][1][iZ*stride_any_direction[Z] + iR*stride_any_direction[R]]
                               = funinv*((2-om*om)*np->P[ec][1][iZ*stride_any_direction[Z]
                                  + iR*stride_any_direction[R]] + (0.5*g-1)*op->P[ec][1][iZ*stride_any_direction[Z]
                                     + iR*stride_any_direction[R]] + np->s[ec][iZ*stride_any_direction[Z]
                                        + iR*stride_any_direction[R]]*f[ec][1][iZ*stride_any_direction[Z]
                                           + iR*stride_any_direction[R]]);
                          }
                        }
                      }
                    }
                  } else {
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                np->s[ec][iZ*stride_any_direction[Z]] = max(-np->energy[ec][iZ*stride_any_direction[Z]]*fac,
                                   0.0);
                                op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                                   + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                                op->P[ec][1][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iZ*stride_any_direction[Z]]
                                   + (0.5*g-1)*op->P[ec][1][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][1][iZ*stride_any_direction[Z]]);
                              }
                            } else {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                  np->s[ec][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                     = max(-np->energy[ec][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*fac,
                                        0.0);
                                  op->P[ec][0][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                     = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]
                                        + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]
                                           + iZ*stride_any_direction[Z]] + np->s[ec][iY*stride_any_direction[Y]
                                              + iZ*stride_any_direction[Z]]*f[ec][0][iY*stride_any_direction[Y]
                                                 + iZ*stride_any_direction[Z]]);
                                  op->P[ec][1][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                     = funinv*((2-om*om)*np->P[ec][1][iY*stride_any_direction[Y]
                                        + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iY*stride_any_direction[Y]
                                           + iZ*stride_any_direction[Z]] + np->s[ec][iY*stride_any_direction[Y]
                                              + iZ*stride_any_direction[Z]]*f[ec][1][iY*stride_any_direction[Y]
                                                 + iZ*stride_any_direction[Z]]);
                                }
                              }
                            }
                          } else {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                              np->s[ec][iZ*stride_any_direction[Z]] = max(-np->energy[ec][iZ*stride_any_direction[Z]]*fac,
                                 0.0);
                              op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                                 + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                              op->P[ec][1][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iZ*stride_any_direction[Z]]
                                 + (0.5*g-1)*op->P[ec][1][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][1][iZ*stride_any_direction[Z]]);
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                np->s[ec][iZ*stride_any_direction[Z]] = max(-np->energy[ec][iZ*stride_any_direction[Z]]*fac,
                                   0.0);
                                op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                                   + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                                op->P[ec][1][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iZ*stride_any_direction[Z]]
                                   + (0.5*g-1)*op->P[ec][1][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][1][iZ*stride_any_direction[Z]]);
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                    np->s[ec][iY + iZ*stride_any_direction[Z]] = max(-np->energy[ec][iY
                                       + iZ*stride_any_direction[Z]]*fac, 0.0);
                                    op->P[ec][0][iY + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iY
                                       + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iY
                                          + iZ*stride_any_direction[Z]] + np->s[ec][iY + iZ*stride_any_direction[Z]]*f[ec][0][iY
                                             + iZ*stride_any_direction[Z]]);
                                    op->P[ec][1][iY + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iY
                                       + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iY
                                          + iZ*stride_any_direction[Z]] + np->s[ec][iY + iZ*stride_any_direction[Z]]*f[ec][1][iY
                                             + iZ*stride_any_direction[Z]]);
                                  }
                                }
                              } else {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                    np->s[ec][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = max(-np->energy[ec][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*fac,
                                          0.0);
                                    op->P[ec][0][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]
                                          + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]
                                             + iZ*stride_any_direction[Z]] + np->s[ec][iY*stride_any_direction[Y]
                                                + iZ*stride_any_direction[Z]]*f[ec][0][iY*stride_any_direction[Y]
                                                   + iZ*stride_any_direction[Z]]);
                                    op->P[ec][1][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = funinv*((2-om*om)*np->P[ec][1][iY*stride_any_direction[Y]
                                          + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iY*stride_any_direction[Y]
                                             + iZ*stride_any_direction[Z]] + np->s[ec][iY*stride_any_direction[Y]
                                                + iZ*stride_any_direction[Z]]*f[ec][1][iY*stride_any_direction[Y]
                                                   + iZ*stride_any_direction[Z]]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                              np->s[ec][iZ*stride_any_direction[Z]] = max(-np->energy[ec][iZ*stride_any_direction[Z]]*fac,
                                 0.0);
                              op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                                 + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                              op->P[ec][1][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iZ*stride_any_direction[Z]]
                                 + (0.5*g-1)*op->P[ec][1][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][1][iZ*stride_any_direction[Z]]);
                            }
                          }
                        }
                      } else {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                  np->s[ec][iX + iZ*stride_any_direction[Z]] = max(-np->energy[ec][iX
                                     + iZ*stride_any_direction[Z]]*fac, 0.0);
                                  op->P[ec][0][iX + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iX
                                     + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX
                                        + iZ*stride_any_direction[Z]] + np->s[ec][iX + iZ*stride_any_direction[Z]]*f[ec][0][iX
                                           + iZ*stride_any_direction[Z]]);
                                  op->P[ec][1][iX + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iX
                                     + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iX
                                        + iZ*stride_any_direction[Z]] + np->s[ec][iX + iZ*stride_any_direction[Z]]*f[ec][1][iX
                                           + iZ*stride_any_direction[Z]]);
                                }
                              }
                            } else {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                    np->s[ec][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = max(-np->energy[ec][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*fac, 0.0);
                                    op->P[ec][0][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = funinv*((2-om*om)*np->P[ec][0][iX + iY*stride_any_direction[Y]
                                          + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX
                                             + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                                + np->s[ec][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*f[ec][0][iX
                                                   + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]);
                                    op->P[ec][1][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = funinv*((2-om*om)*np->P[ec][1][iX + iY*stride_any_direction[Y]
                                          + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iX
                                             + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                                + np->s[ec][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*f[ec][1][iX
                                                   + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                np->s[ec][iX + iZ*stride_any_direction[Z]] = max(-np->energy[ec][iX
                                   + iZ*stride_any_direction[Z]]*fac, 0.0);
                                op->P[ec][0][iX + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iX
                                   + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX
                                      + iZ*stride_any_direction[Z]] + np->s[ec][iX + iZ*stride_any_direction[Z]]*f[ec][0][iX
                                         + iZ*stride_any_direction[Z]]);
                                op->P[ec][1][iX + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iX
                                   + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iX
                                      + iZ*stride_any_direction[Z]] + np->s[ec][iX + iZ*stride_any_direction[Z]]*f[ec][1][iX
                                         + iZ*stride_any_direction[Z]]);
                              }
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                  np->s[ec][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                     = max(-np->energy[ec][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]*fac,
                                        0.0);
                                  op->P[ec][0][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                     = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                        + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                           + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                              + iZ*stride_any_direction[Z]]*f[ec][0][iX*stride_any_direction[X]
                                                 + iZ*stride_any_direction[Z]]);
                                  op->P[ec][1][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                     = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                        + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                           + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                              + iZ*stride_any_direction[Z]]*f[ec][1][iX*stride_any_direction[X]
                                                 + iZ*stride_any_direction[Z]]);
                                }
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                      np->s[ec][iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z]]
                                         = max(-np->energy[ec][iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z]]*fac, 0.0);
                                      op->P[ec][0][iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z]]
                                         = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                            + iY + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                               + iY + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                                  + iY + iZ*stride_any_direction[Z]]*f[ec][0][iX*stride_any_direction[X]
                                                     + iY + iZ*stride_any_direction[Z]]);
                                      op->P[ec][1][iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z]]
                                         = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                            + iY + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                               + iY + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                                  + iY + iZ*stride_any_direction[Z]]*f[ec][1][iX*stride_any_direction[X]
                                                     + iY + iZ*stride_any_direction[Z]]);
                                    }
                                  }
                                }
                              } else {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                      np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                         + iZ*stride_any_direction[Z]] = max(-np->energy[ec][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*fac,
                                               0.0);
                                      op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                         + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                               + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                                  + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                                     + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*f[ec][0][iX*stride_any_direction[X]
                                                        + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]);
                                      op->P[ec][1][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                         + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                               + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                                  + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                                     + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*f[ec][1][iX*stride_any_direction[X]
                                                        + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]);
                                    }
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                np->s[ec][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                   = max(-np->energy[ec][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]*fac,
                                      0.0);
                                op->P[ec][0][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                   = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                      + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                         + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                            + iZ*stride_any_direction[Z]]*f[ec][0][iX*stride_any_direction[X]
                                               + iZ*stride_any_direction[Z]]);
                                op->P[ec][1][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                   = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                      + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                         + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                            + iZ*stride_any_direction[Z]]*f[ec][1][iX*stride_any_direction[X]
                                               + iZ*stride_any_direction[Z]]);
                              }
                            }
                          }
                        }
                      }
                    } else {
                      for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                        np->s[ec][iZ*stride_any_direction[Z]] = max(-np->energy[ec][iZ*stride_any_direction[Z]]*fac,
                           0.0);
                        op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                           + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                        op->P[ec][1][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iZ*stride_any_direction[Z]]
                           + (0.5*g-1)*op->P[ec][1][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][1][iZ*stride_any_direction[Z]]);
                      }
                    }
                  }
                }
              }
            } else {
              if (stride_any_direction[R]) {
                if (num_any_direction[R]==1) {
                  np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                  op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                     + np->s[ec][0]*f[ec][0][0]);
                  op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                     + np->s[ec][0]*f[ec][1][0]);
                } else {
                  if (stride_any_direction[R]==1) {
                    for (int iR=0; iR<num_any_direction[R]; iR += 1) {
                      np->s[ec][iR] = max(-np->energy[ec][iR]*fac, 0.0);
                      op->P[ec][0][iR] = funinv*((2-om*om)*np->P[ec][0][iR] + (0.5*g-1)*op->P[ec][0][iR] + np->s[ec][iR]*f[ec][0][iR]);
                      op->P[ec][1][iR] = funinv*((2-om*om)*np->P[ec][1][iR] + (0.5*g-1)*op->P[ec][1][iR] + np->s[ec][iR]*f[ec][1][iR]);
                    }
                  } else {
                    for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                      np->s[ec][iR*stride_any_direction[R]] = max(-np->energy[ec][iR*stride_any_direction[R]]*fac,
                         0.0);
                      op->P[ec][0][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][0][iR*stride_any_direction[R]]
                         + (0.5*g-1)*op->P[ec][0][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][0][iR*stride_any_direction[R]]);
                      op->P[ec][1][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][1][iR*stride_any_direction[R]]
                         + (0.5*g-1)*op->P[ec][1][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][1][iR*stride_any_direction[R]]);
                    }
                  }
                }
              } else {
                if (stride_any_direction[X]) {
                  if (num_any_direction[X]==1) {
                    if (stride_any_direction[X]==1) {
                      if (stride_any_direction[Y]) {
                        if (num_any_direction[Y]==1) {
                          np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                          op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                             + np->s[ec][0]*f[ec][0][0]);
                          op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                             + np->s[ec][0]*f[ec][1][0]);
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            np->s[ec][iY*stride_any_direction[Y]] = max(-np->energy[ec][iY*stride_any_direction[Y]]*fac,
                               0.0);
                            op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                               + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                            op->P[ec][1][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][1][iY*stride_any_direction[Y]]
                               + (0.5*g-1)*op->P[ec][1][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][1][iY*stride_any_direction[Y]]);
                          }
                        }
                      } else {
                        np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                        op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                           + np->s[ec][0]*f[ec][0][0]);
                        op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                           + np->s[ec][0]*f[ec][1][0]);
                      }
                    } else {
                      if (stride_any_direction[Y]) {
                        if (num_any_direction[Y]==1) {
                          np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                          op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                             + np->s[ec][0]*f[ec][0][0]);
                          op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                             + np->s[ec][0]*f[ec][1][0]);
                        } else {
                          if (stride_any_direction[Y]==1) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              np->s[ec][iY] = max(-np->energy[ec][iY]*fac, 0.0);
                              op->P[ec][0][iY] = funinv*((2-om*om)*np->P[ec][0][iY] + (0.5*g-1)*op->P[ec][0][iY] + np->s[ec][iY]*f[ec][0][iY]);
                              op->P[ec][1][iY] = funinv*((2-om*om)*np->P[ec][1][iY] + (0.5*g-1)*op->P[ec][1][iY] + np->s[ec][iY]*f[ec][1][iY]);
                            }
                          } else {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              np->s[ec][iY*stride_any_direction[Y]] = max(-np->energy[ec][iY*stride_any_direction[Y]]*fac,
                                 0.0);
                              op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                                 + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                              op->P[ec][1][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][1][iY*stride_any_direction[Y]]
                                 + (0.5*g-1)*op->P[ec][1][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][1][iY*stride_any_direction[Y]]);
                            }
                          }
                        }
                      } else {
                        np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                        op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                           + np->s[ec][0]*f[ec][0][0]);
                        op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                           + np->s[ec][0]*f[ec][1][0]);
                      }
                    }
                  } else {
                    if (stride_any_direction[X]==1) {
                      if (stride_any_direction[Y]) {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                            np->s[ec][iX] = max(-np->energy[ec][iX]*fac, 0.0);
                            op->P[ec][0][iX] = funinv*((2-om*om)*np->P[ec][0][iX] + (0.5*g-1)*op->P[ec][0][iX] + np->s[ec][iX]*f[ec][0][iX]);
                            op->P[ec][1][iX] = funinv*((2-om*om)*np->P[ec][1][iX] + (0.5*g-1)*op->P[ec][1][iX] + np->s[ec][iX]*f[ec][1][iX]);
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              np->s[ec][iX + iY*stride_any_direction[Y]] = max(-np->energy[ec][iX
                                 + iY*stride_any_direction[Y]]*fac, 0.0);
                              op->P[ec][0][iX + iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iX
                                 + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX
                                    + iY*stride_any_direction[Y]] + np->s[ec][iX + iY*stride_any_direction[Y]]*f[ec][0][iX
                                       + iY*stride_any_direction[Y]]);
                              op->P[ec][1][iX + iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][1][iX
                                 + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][1][iX
                                    + iY*stride_any_direction[Y]] + np->s[ec][iX + iY*stride_any_direction[Y]]*f[ec][1][iX
                                       + iY*stride_any_direction[Y]]);
                            }
                          }
                        }
                      } else {
                        for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                          np->s[ec][iX] = max(-np->energy[ec][iX]*fac, 0.0);
                          op->P[ec][0][iX] = funinv*((2-om*om)*np->P[ec][0][iX] + (0.5*g-1)*op->P[ec][0][iX] + np->s[ec][iX]*f[ec][0][iX]);
                          op->P[ec][1][iX] = funinv*((2-om*om)*np->P[ec][1][iX] + (0.5*g-1)*op->P[ec][1][iX] + np->s[ec][iX]*f[ec][1][iX]);
                        }
                      }
                    } else {
                      if (stride_any_direction[Y]) {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            np->s[ec][iX*stride_any_direction[X]] = max(-np->energy[ec][iX*stride_any_direction[X]]*fac,
                               0.0);
                            op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                               + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                            op->P[ec][1][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]]
                               + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][1][iX*stride_any_direction[X]]);
                          }
                        } else {
                          if (stride_any_direction[Y]==1) {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                np->s[ec][iX*stride_any_direction[X] + iY] = max(-np->energy[ec][iX*stride_any_direction[X]
                                   + iY]*fac, 0.0);
                                op->P[ec][0][iX*stride_any_direction[X] + iY] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                   + iY] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                      + iY] + np->s[ec][iX*stride_any_direction[X] + iY]*f[ec][0][iX*stride_any_direction[X]
                                         + iY]);
                                op->P[ec][1][iX*stride_any_direction[X] + iY] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                   + iY] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                      + iY] + np->s[ec][iX*stride_any_direction[X] + iY]*f[ec][1][iX*stride_any_direction[X]
                                         + iY]);
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = max(-np->energy[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]*fac,
                                      0.0);
                                op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                      + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                         + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y]]*f[ec][0][iX*stride_any_direction[X]
                                               + iY*stride_any_direction[Y]]);
                                op->P[ec][1][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                      + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                         + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y]]*f[ec][1][iX*stride_any_direction[X]
                                               + iY*stride_any_direction[Y]]);
                              }
                            }
                          }
                        }
                      } else {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          np->s[ec][iX*stride_any_direction[X]] = max(-np->energy[ec][iX*stride_any_direction[X]]*fac,
                             0.0);
                          op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                             + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                          op->P[ec][1][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]]
                             + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][1][iX*stride_any_direction[X]]);
                        }
                      }
                    }
                  }
                } else {
                  np->s[ec][0] = max(-np->energy[ec][0]*fac, 0.0);
                  op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                     + np->s[ec][0]*f[ec][0][0]);
                  op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                     + np->s[ec][0]*f[ec][1][0]);
                }
              }
            }
          } else {
            for (int i=0;i<ntot;i++) {
              np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
              np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
            }
            if (stride_any_direction[Z]) {
              if (num_any_direction[Z]==1) {
                if (stride_any_direction[Z]==1) {
                  if (stride_any_direction[R]) {
                    if (num_any_direction[R]==1) {
                      np->s[ec][0] = np->energy[ec][0]*fac;
                      op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                         + np->s[ec][0]*f[ec][0][0]);
                      op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                         + np->s[ec][0]*f[ec][1][0]);
                    } else {
                      for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                        np->s[ec][iR*stride_any_direction[R]] = np->energy[ec][iR*stride_any_direction[R]]*fac;
                        op->P[ec][0][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][0][iR*stride_any_direction[R]]
                           + (0.5*g-1)*op->P[ec][0][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][0][iR*stride_any_direction[R]]);
                        op->P[ec][1][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][1][iR*stride_any_direction[R]]
                           + (0.5*g-1)*op->P[ec][1][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][1][iR*stride_any_direction[R]]);
                      }
                    }
                  } else {
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[Y]) {
                          if (num_any_direction[Y]==1) {
                            np->s[ec][0] = np->energy[ec][0]*fac;
                            op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                               + np->s[ec][0]*f[ec][0][0]);
                            op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                               + np->s[ec][0]*f[ec][1][0]);
                          } else {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              np->s[ec][iY*stride_any_direction[Y]] = np->energy[ec][iY*stride_any_direction[Y]]*fac;
                              op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                                 + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                              op->P[ec][1][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][1][iY*stride_any_direction[Y]]
                                 + (0.5*g-1)*op->P[ec][1][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][1][iY*stride_any_direction[Y]]);
                            }
                          }
                        } else {
                          np->s[ec][0] = np->energy[ec][0]*fac;
                          op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                             + np->s[ec][0]*f[ec][0][0]);
                          op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                             + np->s[ec][0]*f[ec][1][0]);
                        }
                      } else {
                        if (stride_any_direction[Y]) {
                          if (num_any_direction[Y]==1) {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              np->s[ec][iX*stride_any_direction[X]] = np->energy[ec][iX*stride_any_direction[X]]*fac;
                              op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                                 + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                              op->P[ec][1][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]]
                                 + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][1][iX*stride_any_direction[X]]);
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = np->energy[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]*fac;
                                op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                      + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                         + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y]]*f[ec][0][iX*stride_any_direction[X]
                                               + iY*stride_any_direction[Y]]);
                                op->P[ec][1][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                      + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                         + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y]]*f[ec][1][iX*stride_any_direction[X]
                                               + iY*stride_any_direction[Y]]);
                              }
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            np->s[ec][iX*stride_any_direction[X]] = np->energy[ec][iX*stride_any_direction[X]]*fac;
                            op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                               + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                            op->P[ec][1][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]]
                               + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][1][iX*stride_any_direction[X]]);
                          }
                        }
                      }
                    } else {
                      np->s[ec][0] = np->energy[ec][0]*fac;
                      op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                         + np->s[ec][0]*f[ec][0][0]);
                      op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                         + np->s[ec][0]*f[ec][1][0]);
                    }
                  }
                } else {
                  if (stride_any_direction[R]) {
                    if (num_any_direction[R]==1) {
                      np->s[ec][0] = np->energy[ec][0]*fac;
                      op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                         + np->s[ec][0]*f[ec][0][0]);
                      op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                         + np->s[ec][0]*f[ec][1][0]);
                    } else {
                      if (stride_any_direction[R]==1) {
                        for (int iR=0; iR<num_any_direction[R]; iR += 1) {
                          np->s[ec][iR] = np->energy[ec][iR]*fac;
                          op->P[ec][0][iR] = funinv*((2-om*om)*np->P[ec][0][iR] + (0.5*g-1)*op->P[ec][0][iR] + np->s[ec][iR]*f[ec][0][iR]);
                          op->P[ec][1][iR] = funinv*((2-om*om)*np->P[ec][1][iR] + (0.5*g-1)*op->P[ec][1][iR] + np->s[ec][iR]*f[ec][1][iR]);
                        }
                      } else {
                        for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                          np->s[ec][iR*stride_any_direction[R]] = np->energy[ec][iR*stride_any_direction[R]]*fac;
                          op->P[ec][0][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][0][iR*stride_any_direction[R]]
                             + (0.5*g-1)*op->P[ec][0][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][0][iR*stride_any_direction[R]]);
                          op->P[ec][1][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][1][iR*stride_any_direction[R]]
                             + (0.5*g-1)*op->P[ec][1][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][1][iR*stride_any_direction[R]]);
                        }
                      }
                    }
                  } else {
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              np->s[ec][0] = np->energy[ec][0]*fac;
                              op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                                 + np->s[ec][0]*f[ec][0][0]);
                              op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                                 + np->s[ec][0]*f[ec][1][0]);
                            } else {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                np->s[ec][iY*stride_any_direction[Y]] = np->energy[ec][iY*stride_any_direction[Y]]*fac;
                                op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                                   + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                                op->P[ec][1][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][1][iY*stride_any_direction[Y]]
                                   + (0.5*g-1)*op->P[ec][1][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][1][iY*stride_any_direction[Y]]);
                              }
                            }
                          } else {
                            np->s[ec][0] = np->energy[ec][0]*fac;
                            op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                               + np->s[ec][0]*f[ec][0][0]);
                            op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                               + np->s[ec][0]*f[ec][1][0]);
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              np->s[ec][0] = np->energy[ec][0]*fac;
                              op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                                 + np->s[ec][0]*f[ec][0][0]);
                              op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                                 + np->s[ec][0]*f[ec][1][0]);
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                  np->s[ec][iY] = np->energy[ec][iY]*fac;
                                  op->P[ec][0][iY] = funinv*((2-om*om)*np->P[ec][0][iY] + (0.5*g-1)*op->P[ec][0][iY] + np->s[ec][iY]*f[ec][0][iY]);
                                  op->P[ec][1][iY] = funinv*((2-om*om)*np->P[ec][1][iY] + (0.5*g-1)*op->P[ec][1][iY] + np->s[ec][iY]*f[ec][1][iY]);
                                }
                              } else {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  np->s[ec][iY*stride_any_direction[Y]] = np->energy[ec][iY*stride_any_direction[Y]]*fac;
                                  op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                                     + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                                  op->P[ec][1][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][1][iY*stride_any_direction[Y]]
                                     + (0.5*g-1)*op->P[ec][1][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][1][iY*stride_any_direction[Y]]);
                                }
                              }
                            }
                          } else {
                            np->s[ec][0] = np->energy[ec][0]*fac;
                            op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                               + np->s[ec][0]*f[ec][0][0]);
                            op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                               + np->s[ec][0]*f[ec][1][0]);
                          }
                        }
                      } else {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                np->s[ec][iX] = np->energy[ec][iX]*fac;
                                op->P[ec][0][iX] = funinv*((2-om*om)*np->P[ec][0][iX] + (0.5*g-1)*op->P[ec][0][iX] + np->s[ec][iX]*f[ec][0][iX]);
                                op->P[ec][1][iX] = funinv*((2-om*om)*np->P[ec][1][iX] + (0.5*g-1)*op->P[ec][1][iX] + np->s[ec][iX]*f[ec][1][iX]);
                              }
                            } else {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  np->s[ec][iX + iY*stride_any_direction[Y]] = np->energy[ec][iX
                                     + iY*stride_any_direction[Y]]*fac;
                                  op->P[ec][0][iX + iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iX
                                     + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX
                                        + iY*stride_any_direction[Y]] + np->s[ec][iX + iY*stride_any_direction[Y]]*f[ec][0][iX
                                           + iY*stride_any_direction[Y]]);
                                  op->P[ec][1][iX + iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][1][iX
                                     + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][1][iX
                                        + iY*stride_any_direction[Y]] + np->s[ec][iX + iY*stride_any_direction[Y]]*f[ec][1][iX
                                           + iY*stride_any_direction[Y]]);
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                              np->s[ec][iX] = np->energy[ec][iX]*fac;
                              op->P[ec][0][iX] = funinv*((2-om*om)*np->P[ec][0][iX] + (0.5*g-1)*op->P[ec][0][iX] + np->s[ec][iX]*f[ec][0][iX]);
                              op->P[ec][1][iX] = funinv*((2-om*om)*np->P[ec][1][iX] + (0.5*g-1)*op->P[ec][1][iX] + np->s[ec][iX]*f[ec][1][iX]);
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                np->s[ec][iX*stride_any_direction[X]] = np->energy[ec][iX*stride_any_direction[X]]*fac;
                                op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                                   + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                                op->P[ec][1][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]]
                                   + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][1][iX*stride_any_direction[X]]);
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                    np->s[ec][iX*stride_any_direction[X] + iY] = np->energy[ec][iX*stride_any_direction[X]
                                       + iY]*fac;
                                    op->P[ec][0][iX*stride_any_direction[X] + iY] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                       + iY] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                          + iY] + np->s[ec][iX*stride_any_direction[X] + iY]*f[ec][0][iX*stride_any_direction[X]
                                             + iY]);
                                    op->P[ec][1][iX*stride_any_direction[X] + iY] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                       + iY] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                          + iY] + np->s[ec][iX*stride_any_direction[X] + iY]*f[ec][1][iX*stride_any_direction[X]
                                             + iY]);
                                  }
                                }
                              } else {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                    np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                       = np->energy[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]*fac;
                                    op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                       = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                          + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                             + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                                + iY*stride_any_direction[Y]]*f[ec][0][iX*stride_any_direction[X]
                                                   + iY*stride_any_direction[Y]]);
                                    op->P[ec][1][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                       = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                          + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                             + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                                + iY*stride_any_direction[Y]]*f[ec][1][iX*stride_any_direction[X]
                                                   + iY*stride_any_direction[Y]]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              np->s[ec][iX*stride_any_direction[X]] = np->energy[ec][iX*stride_any_direction[X]]*fac;
                              op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                                 + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                              op->P[ec][1][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]]
                                 + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][1][iX*stride_any_direction[X]]);
                            }
                          }
                        }
                      }
                    } else {
                      np->s[ec][0] = np->energy[ec][0]*fac;
                      op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                         + np->s[ec][0]*f[ec][0][0]);
                      op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                         + np->s[ec][0]*f[ec][1][0]);
                    }
                  }
                }
              } else {
                if (stride_any_direction[Z]==1) {
                  if (stride_any_direction[R]) {
                    if (num_any_direction[R]==1) {
                      for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                        np->s[ec][iZ] = np->energy[ec][iZ]*fac;
                        op->P[ec][0][iZ] = funinv*((2-om*om)*np->P[ec][0][iZ] + (0.5*g-1)*op->P[ec][0][iZ] + np->s[ec][iZ]*f[ec][0][iZ]);
                        op->P[ec][1][iZ] = funinv*((2-om*om)*np->P[ec][1][iZ] + (0.5*g-1)*op->P[ec][1][iZ] + np->s[ec][iZ]*f[ec][1][iZ]);
                      }
                    } else {
                      for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                        for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                          np->s[ec][iZ + iR*stride_any_direction[R]] = np->energy[ec][iZ
                             + iR*stride_any_direction[R]]*fac;
                          op->P[ec][0][iZ + iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][0][iZ
                             + iR*stride_any_direction[R]] + (0.5*g-1)*op->P[ec][0][iZ
                                + iR*stride_any_direction[R]] + np->s[ec][iZ + iR*stride_any_direction[R]]*f[ec][0][iZ
                                   + iR*stride_any_direction[R]]);
                          op->P[ec][1][iZ + iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][1][iZ
                             + iR*stride_any_direction[R]] + (0.5*g-1)*op->P[ec][1][iZ
                                + iR*stride_any_direction[R]] + np->s[ec][iZ + iR*stride_any_direction[R]]*f[ec][1][iZ
                                   + iR*stride_any_direction[R]]);
                        }
                      }
                    }
                  } else {
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[Y]) {
                          if (num_any_direction[Y]==1) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              np->s[ec][iZ] = np->energy[ec][iZ]*fac;
                              op->P[ec][0][iZ] = funinv*((2-om*om)*np->P[ec][0][iZ] + (0.5*g-1)*op->P[ec][0][iZ] + np->s[ec][iZ]*f[ec][0][iZ]);
                              op->P[ec][1][iZ] = funinv*((2-om*om)*np->P[ec][1][iZ] + (0.5*g-1)*op->P[ec][1][iZ] + np->s[ec][iZ]*f[ec][1][iZ]);
                            }
                          } else {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                np->s[ec][iY*stride_any_direction[Y] + iZ] = np->energy[ec][iY*stride_any_direction[Y]
                                   + iZ]*fac;
                                op->P[ec][0][iY*stride_any_direction[Y] + iZ] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]
                                   + iZ] + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]
                                      + iZ] + np->s[ec][iY*stride_any_direction[Y] + iZ]*f[ec][0][iY*stride_any_direction[Y]
                                         + iZ]);
                                op->P[ec][1][iY*stride_any_direction[Y] + iZ] = funinv*((2-om*om)*np->P[ec][1][iY*stride_any_direction[Y]
                                   + iZ] + (0.5*g-1)*op->P[ec][1][iY*stride_any_direction[Y]
                                      + iZ] + np->s[ec][iY*stride_any_direction[Y] + iZ]*f[ec][1][iY*stride_any_direction[Y]
                                         + iZ]);
                              }
                            }
                          }
                        } else {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                            np->s[ec][iZ] = np->energy[ec][iZ]*fac;
                            op->P[ec][0][iZ] = funinv*((2-om*om)*np->P[ec][0][iZ] + (0.5*g-1)*op->P[ec][0][iZ] + np->s[ec][iZ]*f[ec][0][iZ]);
                            op->P[ec][1][iZ] = funinv*((2-om*om)*np->P[ec][1][iZ] + (0.5*g-1)*op->P[ec][1][iZ] + np->s[ec][iZ]*f[ec][1][iZ]);
                          }
                        }
                      } else {
                        if (stride_any_direction[Y]) {
                          if (num_any_direction[Y]==1) {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                np->s[ec][iX*stride_any_direction[X] + iZ] = np->energy[ec][iX*stride_any_direction[X]
                                   + iZ]*fac;
                                op->P[ec][0][iX*stride_any_direction[X] + iZ] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                   + iZ] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                      + iZ] + np->s[ec][iX*stride_any_direction[X] + iZ]*f[ec][0][iX*stride_any_direction[X]
                                         + iZ]);
                                op->P[ec][1][iX*stride_any_direction[X] + iZ] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                   + iZ] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                      + iZ] + np->s[ec][iX*stride_any_direction[X] + iZ]*f[ec][1][iX*stride_any_direction[X]
                                         + iZ]);
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                  np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                     + iZ] = np->energy[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                        + iZ]*fac;
                                  op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                     + iZ] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                        + iY*stride_any_direction[Y] + iZ] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                           + iY*stride_any_direction[Y] + iZ] + np->s[ec][iX*stride_any_direction[X]
                                              + iY*stride_any_direction[Y] + iZ]*f[ec][0][iX*stride_any_direction[X]
                                                 + iY*stride_any_direction[Y] + iZ]);
                                  op->P[ec][1][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                     + iZ] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                        + iY*stride_any_direction[Y] + iZ] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                           + iY*stride_any_direction[Y] + iZ] + np->s[ec][iX*stride_any_direction[X]
                                              + iY*stride_any_direction[Y] + iZ]*f[ec][1][iX*stride_any_direction[X]
                                                 + iY*stride_any_direction[Y] + iZ]);
                                }
                              }
                            }
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              np->s[ec][iX*stride_any_direction[X] + iZ] = np->energy[ec][iX*stride_any_direction[X]
                                 + iZ]*fac;
                              op->P[ec][0][iX*stride_any_direction[X] + iZ] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                 + iZ] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                    + iZ] + np->s[ec][iX*stride_any_direction[X] + iZ]*f[ec][0][iX*stride_any_direction[X]
                                       + iZ]);
                              op->P[ec][1][iX*stride_any_direction[X] + iZ] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                 + iZ] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                    + iZ] + np->s[ec][iX*stride_any_direction[X] + iZ]*f[ec][1][iX*stride_any_direction[X]
                                       + iZ]);
                            }
                          }
                        }
                      }
                    } else {
                      for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                        np->s[ec][iZ] = np->energy[ec][iZ]*fac;
                        op->P[ec][0][iZ] = funinv*((2-om*om)*np->P[ec][0][iZ] + (0.5*g-1)*op->P[ec][0][iZ] + np->s[ec][iZ]*f[ec][0][iZ]);
                        op->P[ec][1][iZ] = funinv*((2-om*om)*np->P[ec][1][iZ] + (0.5*g-1)*op->P[ec][1][iZ] + np->s[ec][iZ]*f[ec][1][iZ]);
                      }
                    }
                  }
                } else {
                  if (stride_any_direction[R]) {
                    if (num_any_direction[R]==1) {
                      for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                        np->s[ec][iZ*stride_any_direction[Z]] = np->energy[ec][iZ*stride_any_direction[Z]]*fac;
                        op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                           + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                        op->P[ec][1][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iZ*stride_any_direction[Z]]
                           + (0.5*g-1)*op->P[ec][1][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][1][iZ*stride_any_direction[Z]]);
                      }
                    } else {
                      if (stride_any_direction[R]==1) {
                        for (int iR=0; iR<num_any_direction[R]; iR += 1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                            np->s[ec][iZ*stride_any_direction[Z] + iR] = np->energy[ec][iZ*stride_any_direction[Z]
                               + iR]*fac;
                            op->P[ec][0][iZ*stride_any_direction[Z] + iR] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]
                               + iR] + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]
                                  + iR] + np->s[ec][iZ*stride_any_direction[Z] + iR]*f[ec][0][iZ*stride_any_direction[Z]
                                     + iR]);
                            op->P[ec][1][iZ*stride_any_direction[Z] + iR] = funinv*((2-om*om)*np->P[ec][1][iZ*stride_any_direction[Z]
                               + iR] + (0.5*g-1)*op->P[ec][1][iZ*stride_any_direction[Z]
                                  + iR] + np->s[ec][iZ*stride_any_direction[Z] + iR]*f[ec][1][iZ*stride_any_direction[Z]
                                     + iR]);
                          }
                        }
                      } else {
                        for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                            np->s[ec][iZ*stride_any_direction[Z] + iR*stride_any_direction[R]]
                               = np->energy[ec][iZ*stride_any_direction[Z] + iR*stride_any_direction[R]]*fac;
                            op->P[ec][0][iZ*stride_any_direction[Z] + iR*stride_any_direction[R]]
                               = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]
                                  + iR*stride_any_direction[R]] + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]
                                     + iR*stride_any_direction[R]] + np->s[ec][iZ*stride_any_direction[Z]
                                        + iR*stride_any_direction[R]]*f[ec][0][iZ*stride_any_direction[Z]
                                           + iR*stride_any_direction[R]]);
                            op->P[ec][1][iZ*stride_any_direction[Z] + iR*stride_any_direction[R]]
                               = funinv*((2-om*om)*np->P[ec][1][iZ*stride_any_direction[Z]
                                  + iR*stride_any_direction[R]] + (0.5*g-1)*op->P[ec][1][iZ*stride_any_direction[Z]
                                     + iR*stride_any_direction[R]] + np->s[ec][iZ*stride_any_direction[Z]
                                        + iR*stride_any_direction[R]]*f[ec][1][iZ*stride_any_direction[Z]
                                           + iR*stride_any_direction[R]]);
                          }
                        }
                      }
                    }
                  } else {
                    if (stride_any_direction[X]) {
                      if (num_any_direction[X]==1) {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                np->s[ec][iZ*stride_any_direction[Z]] = np->energy[ec][iZ*stride_any_direction[Z]]*fac;
                                op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                                   + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                                op->P[ec][1][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iZ*stride_any_direction[Z]]
                                   + (0.5*g-1)*op->P[ec][1][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][1][iZ*stride_any_direction[Z]]);
                              }
                            } else {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                  np->s[ec][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                     = np->energy[ec][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*fac;
                                  op->P[ec][0][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                     = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]
                                        + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]
                                           + iZ*stride_any_direction[Z]] + np->s[ec][iY*stride_any_direction[Y]
                                              + iZ*stride_any_direction[Z]]*f[ec][0][iY*stride_any_direction[Y]
                                                 + iZ*stride_any_direction[Z]]);
                                  op->P[ec][1][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                     = funinv*((2-om*om)*np->P[ec][1][iY*stride_any_direction[Y]
                                        + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iY*stride_any_direction[Y]
                                           + iZ*stride_any_direction[Z]] + np->s[ec][iY*stride_any_direction[Y]
                                              + iZ*stride_any_direction[Z]]*f[ec][1][iY*stride_any_direction[Y]
                                                 + iZ*stride_any_direction[Z]]);
                                }
                              }
                            }
                          } else {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                              np->s[ec][iZ*stride_any_direction[Z]] = np->energy[ec][iZ*stride_any_direction[Z]]*fac;
                              op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                                 + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                              op->P[ec][1][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iZ*stride_any_direction[Z]]
                                 + (0.5*g-1)*op->P[ec][1][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][1][iZ*stride_any_direction[Z]]);
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                np->s[ec][iZ*stride_any_direction[Z]] = np->energy[ec][iZ*stride_any_direction[Z]]*fac;
                                op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                                   + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                                op->P[ec][1][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iZ*stride_any_direction[Z]]
                                   + (0.5*g-1)*op->P[ec][1][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][1][iZ*stride_any_direction[Z]]);
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                    np->s[ec][iY + iZ*stride_any_direction[Z]] = np->energy[ec][iY
                                       + iZ*stride_any_direction[Z]]*fac;
                                    op->P[ec][0][iY + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iY
                                       + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iY
                                          + iZ*stride_any_direction[Z]] + np->s[ec][iY + iZ*stride_any_direction[Z]]*f[ec][0][iY
                                             + iZ*stride_any_direction[Z]]);
                                    op->P[ec][1][iY + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iY
                                       + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iY
                                          + iZ*stride_any_direction[Z]] + np->s[ec][iY + iZ*stride_any_direction[Z]]*f[ec][1][iY
                                             + iZ*stride_any_direction[Z]]);
                                  }
                                }
                              } else {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                    np->s[ec][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = np->energy[ec][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*fac;
                                    op->P[ec][0][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]
                                          + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]
                                             + iZ*stride_any_direction[Z]] + np->s[ec][iY*stride_any_direction[Y]
                                                + iZ*stride_any_direction[Z]]*f[ec][0][iY*stride_any_direction[Y]
                                                   + iZ*stride_any_direction[Z]]);
                                    op->P[ec][1][iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = funinv*((2-om*om)*np->P[ec][1][iY*stride_any_direction[Y]
                                          + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iY*stride_any_direction[Y]
                                             + iZ*stride_any_direction[Z]] + np->s[ec][iY*stride_any_direction[Y]
                                                + iZ*stride_any_direction[Z]]*f[ec][1][iY*stride_any_direction[Y]
                                                   + iZ*stride_any_direction[Z]]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                              np->s[ec][iZ*stride_any_direction[Z]] = np->energy[ec][iZ*stride_any_direction[Z]]*fac;
                              op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                                 + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                              op->P[ec][1][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iZ*stride_any_direction[Z]]
                                 + (0.5*g-1)*op->P[ec][1][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][1][iZ*stride_any_direction[Z]]);
                            }
                          }
                        }
                      } else {
                        if (stride_any_direction[X]==1) {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                  np->s[ec][iX + iZ*stride_any_direction[Z]] = np->energy[ec][iX
                                     + iZ*stride_any_direction[Z]]*fac;
                                  op->P[ec][0][iX + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iX
                                     + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX
                                        + iZ*stride_any_direction[Z]] + np->s[ec][iX + iZ*stride_any_direction[Z]]*f[ec][0][iX
                                           + iZ*stride_any_direction[Z]]);
                                  op->P[ec][1][iX + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iX
                                     + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iX
                                        + iZ*stride_any_direction[Z]] + np->s[ec][iX + iZ*stride_any_direction[Z]]*f[ec][1][iX
                                           + iZ*stride_any_direction[Z]]);
                                }
                              }
                            } else {
                              for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                                for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                    np->s[ec][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = np->energy[ec][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*fac;
                                    op->P[ec][0][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = funinv*((2-om*om)*np->P[ec][0][iX + iY*stride_any_direction[Y]
                                          + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX
                                             + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                                + np->s[ec][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*f[ec][0][iX
                                                   + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]);
                                    op->P[ec][1][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                       = funinv*((2-om*om)*np->P[ec][1][iX + iY*stride_any_direction[Y]
                                          + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iX
                                             + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                                + np->s[ec][iX + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*f[ec][1][iX
                                                   + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]);
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                np->s[ec][iX + iZ*stride_any_direction[Z]] = np->energy[ec][iX
                                   + iZ*stride_any_direction[Z]]*fac;
                                op->P[ec][0][iX + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iX
                                   + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX
                                      + iZ*stride_any_direction[Z]] + np->s[ec][iX + iZ*stride_any_direction[Z]]*f[ec][0][iX
                                         + iZ*stride_any_direction[Z]]);
                                op->P[ec][1][iX + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iX
                                   + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iX
                                      + iZ*stride_any_direction[Z]] + np->s[ec][iX + iZ*stride_any_direction[Z]]*f[ec][1][iX
                                         + iZ*stride_any_direction[Z]]);
                              }
                            }
                          }
                        } else {
                          if (stride_any_direction[Y]) {
                            if (num_any_direction[Y]==1) {
                              for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                  np->s[ec][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                     = np->energy[ec][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]*fac;
                                  op->P[ec][0][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                     = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                        + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                           + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                              + iZ*stride_any_direction[Z]]*f[ec][0][iX*stride_any_direction[X]
                                                 + iZ*stride_any_direction[Z]]);
                                  op->P[ec][1][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                     = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                        + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                           + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                              + iZ*stride_any_direction[Z]]*f[ec][1][iX*stride_any_direction[X]
                                                 + iZ*stride_any_direction[Z]]);
                                }
                              }
                            } else {
                              if (stride_any_direction[Y]==1) {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                      np->s[ec][iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z]]
                                         = np->energy[ec][iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z]]*fac;
                                      op->P[ec][0][iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z]]
                                         = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                            + iY + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                               + iY + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                                  + iY + iZ*stride_any_direction[Z]]*f[ec][0][iX*stride_any_direction[X]
                                                     + iY + iZ*stride_any_direction[Z]]);
                                      op->P[ec][1][iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z]]
                                         = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                            + iY + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                               + iY + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                                  + iY + iZ*stride_any_direction[Z]]*f[ec][1][iX*stride_any_direction[X]
                                                     + iY + iZ*stride_any_direction[Z]]);
                                    }
                                  }
                                }
                              } else {
                                for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                                  for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                      np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                         + iZ*stride_any_direction[Z]] = np->energy[ec][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*fac;
                                      op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                         + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                               + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                                  + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                                     + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*f[ec][0][iX*stride_any_direction[X]
                                                        + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]);
                                      op->P[ec][1][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                         + iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]
                                               + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                                  + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                                     + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]*f[ec][1][iX*stride_any_direction[X]
                                                        + iY*stride_any_direction[Y] + iZ*stride_any_direction[Z]]);
                                    }
                                  }
                                }
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                np->s[ec][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                   = np->energy[ec][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]*fac;
                                op->P[ec][0][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                   = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                      + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                         + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                            + iZ*stride_any_direction[Z]]*f[ec][0][iX*stride_any_direction[X]
                                               + iZ*stride_any_direction[Z]]);
                                op->P[ec][1][iX*stride_any_direction[X] + iZ*stride_any_direction[Z]]
                                   = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                      + iZ*stride_any_direction[Z]] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                         + iZ*stride_any_direction[Z]] + np->s[ec][iX*stride_any_direction[X]
                                            + iZ*stride_any_direction[Z]]*f[ec][1][iX*stride_any_direction[X]
                                               + iZ*stride_any_direction[Z]]);
                              }
                            }
                          }
                        }
                      }
                    } else {
                      for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                        np->s[ec][iZ*stride_any_direction[Z]] = np->energy[ec][iZ*stride_any_direction[Z]]*fac;
                        op->P[ec][0][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][0][iZ*stride_any_direction[Z]]
                           + (0.5*g-1)*op->P[ec][0][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][0][iZ*stride_any_direction[Z]]);
                        op->P[ec][1][iZ*stride_any_direction[Z]] = funinv*((2-om*om)*np->P[ec][1][iZ*stride_any_direction[Z]]
                           + (0.5*g-1)*op->P[ec][1][iZ*stride_any_direction[Z]] + np->s[ec][iZ*stride_any_direction[Z]]*f[ec][1][iZ*stride_any_direction[Z]]);
                      }
                    }
                  }
                }
              }
            } else {
              if (stride_any_direction[R]) {
                if (num_any_direction[R]==1) {
                  np->s[ec][0] = np->energy[ec][0]*fac;
                  op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                     + np->s[ec][0]*f[ec][0][0]);
                  op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                     + np->s[ec][0]*f[ec][1][0]);
                } else {
                  if (stride_any_direction[R]==1) {
                    for (int iR=0; iR<num_any_direction[R]; iR += 1) {
                      np->s[ec][iR] = np->energy[ec][iR]*fac;
                      op->P[ec][0][iR] = funinv*((2-om*om)*np->P[ec][0][iR] + (0.5*g-1)*op->P[ec][0][iR] + np->s[ec][iR]*f[ec][0][iR]);
                      op->P[ec][1][iR] = funinv*((2-om*om)*np->P[ec][1][iR] + (0.5*g-1)*op->P[ec][1][iR] + np->s[ec][iR]*f[ec][1][iR]);
                    }
                  } else {
                    for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                      np->s[ec][iR*stride_any_direction[R]] = np->energy[ec][iR*stride_any_direction[R]]*fac;
                      op->P[ec][0][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][0][iR*stride_any_direction[R]]
                         + (0.5*g-1)*op->P[ec][0][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][0][iR*stride_any_direction[R]]);
                      op->P[ec][1][iR*stride_any_direction[R]] = funinv*((2-om*om)*np->P[ec][1][iR*stride_any_direction[R]]
                         + (0.5*g-1)*op->P[ec][1][iR*stride_any_direction[R]] + np->s[ec][iR*stride_any_direction[R]]*f[ec][1][iR*stride_any_direction[R]]);
                    }
                  }
                }
              } else {
                if (stride_any_direction[X]) {
                  if (num_any_direction[X]==1) {
                    if (stride_any_direction[X]==1) {
                      if (stride_any_direction[Y]) {
                        if (num_any_direction[Y]==1) {
                          np->s[ec][0] = np->energy[ec][0]*fac;
                          op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                             + np->s[ec][0]*f[ec][0][0]);
                          op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                             + np->s[ec][0]*f[ec][1][0]);
                        } else {
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            np->s[ec][iY*stride_any_direction[Y]] = np->energy[ec][iY*stride_any_direction[Y]]*fac;
                            op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                               + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                            op->P[ec][1][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][1][iY*stride_any_direction[Y]]
                               + (0.5*g-1)*op->P[ec][1][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][1][iY*stride_any_direction[Y]]);
                          }
                        }
                      } else {
                        np->s[ec][0] = np->energy[ec][0]*fac;
                        op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                           + np->s[ec][0]*f[ec][0][0]);
                        op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                           + np->s[ec][0]*f[ec][1][0]);
                      }
                    } else {
                      if (stride_any_direction[Y]) {
                        if (num_any_direction[Y]==1) {
                          np->s[ec][0] = np->energy[ec][0]*fac;
                          op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                             + np->s[ec][0]*f[ec][0][0]);
                          op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                             + np->s[ec][0]*f[ec][1][0]);
                        } else {
                          if (stride_any_direction[Y]==1) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              np->s[ec][iY] = np->energy[ec][iY]*fac;
                              op->P[ec][0][iY] = funinv*((2-om*om)*np->P[ec][0][iY] + (0.5*g-1)*op->P[ec][0][iY] + np->s[ec][iY]*f[ec][0][iY]);
                              op->P[ec][1][iY] = funinv*((2-om*om)*np->P[ec][1][iY] + (0.5*g-1)*op->P[ec][1][iY] + np->s[ec][iY]*f[ec][1][iY]);
                            }
                          } else {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              np->s[ec][iY*stride_any_direction[Y]] = np->energy[ec][iY*stride_any_direction[Y]]*fac;
                              op->P[ec][0][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iY*stride_any_direction[Y]]
                                 + (0.5*g-1)*op->P[ec][0][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][0][iY*stride_any_direction[Y]]);
                              op->P[ec][1][iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][1][iY*stride_any_direction[Y]]
                                 + (0.5*g-1)*op->P[ec][1][iY*stride_any_direction[Y]] + np->s[ec][iY*stride_any_direction[Y]]*f[ec][1][iY*stride_any_direction[Y]]);
                            }
                          }
                        }
                      } else {
                        np->s[ec][0] = np->energy[ec][0]*fac;
                        op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                           + np->s[ec][0]*f[ec][0][0]);
                        op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                           + np->s[ec][0]*f[ec][1][0]);
                      }
                    }
                  } else {
                    if (stride_any_direction[X]==1) {
                      if (stride_any_direction[Y]) {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                            np->s[ec][iX] = np->energy[ec][iX]*fac;
                            op->P[ec][0][iX] = funinv*((2-om*om)*np->P[ec][0][iX] + (0.5*g-1)*op->P[ec][0][iX] + np->s[ec][iX]*f[ec][0][iX]);
                            op->P[ec][1][iX] = funinv*((2-om*om)*np->P[ec][1][iX] + (0.5*g-1)*op->P[ec][1][iX] + np->s[ec][iX]*f[ec][1][iX]);
                          }
                        } else {
                          for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              np->s[ec][iX + iY*stride_any_direction[Y]] = np->energy[ec][iX
                                 + iY*stride_any_direction[Y]]*fac;
                              op->P[ec][0][iX + iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][0][iX
                                 + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX
                                    + iY*stride_any_direction[Y]] + np->s[ec][iX + iY*stride_any_direction[Y]]*f[ec][0][iX
                                       + iY*stride_any_direction[Y]]);
                              op->P[ec][1][iX + iY*stride_any_direction[Y]] = funinv*((2-om*om)*np->P[ec][1][iX
                                 + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][1][iX
                                    + iY*stride_any_direction[Y]] + np->s[ec][iX + iY*stride_any_direction[Y]]*f[ec][1][iX
                                       + iY*stride_any_direction[Y]]);
                            }
                          }
                        }
                      } else {
                        for (int iX=0; iX<num_any_direction[X]; iX += 1) {
                          np->s[ec][iX] = np->energy[ec][iX]*fac;
                          op->P[ec][0][iX] = funinv*((2-om*om)*np->P[ec][0][iX] + (0.5*g-1)*op->P[ec][0][iX] + np->s[ec][iX]*f[ec][0][iX]);
                          op->P[ec][1][iX] = funinv*((2-om*om)*np->P[ec][1][iX] + (0.5*g-1)*op->P[ec][1][iX] + np->s[ec][iX]*f[ec][1][iX]);
                        }
                      }
                    } else {
                      if (stride_any_direction[Y]) {
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            np->s[ec][iX*stride_any_direction[X]] = np->energy[ec][iX*stride_any_direction[X]]*fac;
                            op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                               + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                            op->P[ec][1][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]]
                               + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][1][iX*stride_any_direction[X]]);
                          }
                        } else {
                          if (stride_any_direction[Y]==1) {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                                np->s[ec][iX*stride_any_direction[X] + iY] = np->energy[ec][iX*stride_any_direction[X]
                                   + iY]*fac;
                                op->P[ec][0][iX*stride_any_direction[X] + iY] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                   + iY] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                      + iY] + np->s[ec][iX*stride_any_direction[X] + iY]*f[ec][0][iX*stride_any_direction[X]
                                         + iY]);
                                op->P[ec][1][iX*stride_any_direction[X] + iY] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                   + iY] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                      + iY] + np->s[ec][iX*stride_any_direction[X] + iY]*f[ec][1][iX*stride_any_direction[X]
                                         + iY]);
                              }
                            }
                          } else {
                            for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                              for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                                np->s[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = np->energy[ec][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]*fac;
                                op->P[ec][0][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]
                                      + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]
                                         + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y]]*f[ec][0][iX*stride_any_direction[X]
                                               + iY*stride_any_direction[Y]]);
                                op->P[ec][1][iX*stride_any_direction[X] + iY*stride_any_direction[Y]]
                                   = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]
                                      + iY*stride_any_direction[Y]] + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]
                                         + iY*stride_any_direction[Y]] + np->s[ec][iX*stride_any_direction[X]
                                            + iY*stride_any_direction[Y]]*f[ec][1][iX*stride_any_direction[X]
                                               + iY*stride_any_direction[Y]]);
                              }
                            }
                          }
                        }
                      } else {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          np->s[ec][iX*stride_any_direction[X]] = np->energy[ec][iX*stride_any_direction[X]]*fac;
                          op->P[ec][0][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][0][iX*stride_any_direction[X]]
                             + (0.5*g-1)*op->P[ec][0][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][0][iX*stride_any_direction[X]]);
                          op->P[ec][1][iX*stride_any_direction[X]] = funinv*((2-om*om)*np->P[ec][1][iX*stride_any_direction[X]]
                             + (0.5*g-1)*op->P[ec][1][iX*stride_any_direction[X]] + np->s[ec][iX*stride_any_direction[X]]*f[ec][1][iX*stride_any_direction[X]]);
                        }
                      }
                    }
                  }
                } else {
                  np->s[ec][0] = np->energy[ec][0]*fac;
                  op->P[ec][0][0] = funinv*((2-om*om)*np->P[ec][0][0] + (0.5*g-1)*op->P[ec][0][0]
                     + np->s[ec][0]*f[ec][0][0]);
                  op->P[ec][1][0] = funinv*((2-om*om)*np->P[ec][1][0] + (0.5*g-1)*op->P[ec][1][0]
                     + np->s[ec][0]*f[ec][1][0]);
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
