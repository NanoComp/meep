if (f[Er][0]) {
  // Am in cylindrical coordinates.
  FOR_ELECTRIC_COMPONENTS(ec) if (f[ec][0]) {
    const int yee_idx = v.yee_index(ec);
    const int s_ec = stride_any_direction[component_direction(ec)];
    const direction d_1 = (direction)(((ec-2)+1)%3+2);
    const component ec_1 = direction_component(ec,d_1);
    const int s_1 = stride_any_direction[d_1];
    const direction d_2 = (direction)(((ec-2)+2)%3+2);
    const component ec_2 = direction_component(ec,d_2);
    const int s_2 = stride_any_direction[d_2];
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
                  const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                     + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                        + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                           + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                              + np->energy[ec_2][i-s_2+s_ec]);
                  np->s[ec][i] = max(-energy_here*fac, 0.0);
                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                    const int i = yee_idx + iR*stride_any_direction[R];
                    const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                       + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                          + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                             + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                + np->energy[ec_2][i-s_2+s_ec]);
                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                    const int i = yee_idx + iZ;
                    const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                       + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                          + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                             + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                + np->energy[ec_2][i-s_2+s_ec]);
                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                         + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                            + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                               + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                  + np->energy[ec_2][i-s_2+s_ec]);
                      np->s[ec][i] = max(-energy_here*fac, 0.0);
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    }
                  }
                }
              }
            } else { // not fac > 0
              for (int i=0;i<ntot;i++) {
                np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
              }
              if (num_any_direction[Z]==1) {
                if (num_any_direction[R]==1) {
                  const int i = yee_idx;
                  const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                     + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                        + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                           + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                              + np->energy[ec_2][i-s_2+s_ec]);
                  np->s[ec][i] = energy_here*fac;
                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                    const int i = yee_idx + iR*stride_any_direction[R];
                    const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                       + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                          + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                             + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                + np->energy[ec_2][i-s_2+s_ec]);
                    np->s[ec][i] = energy_here*fac;
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                    const int i = yee_idx + iZ;
                    const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                       + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                          + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                             + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                + np->energy[ec_2][i-s_2+s_ec]);
                    np->s[ec][i] = energy_here*fac;
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                         + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                            + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                               + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                  + np->energy[ec_2][i-s_2+s_ec]);
                      np->s[ec][i] = energy_here*fac;
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    }
                  }
                }
              }
            }
          } else { // not fac
            for (int i=0;i<ntot;i++) {
              np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                 + np->s[ec][i]*f[ec][0][i]);
            }
          }
        }
      }
    } else { // not is_real
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
                  const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                     + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                        + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                           + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                              + np->energy[ec_2][i-s_2+s_ec]);
                  np->s[ec][i] = max(-energy_here*fac, 0.0);
                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                     + np->s[ec][i]*f[ec][1][i]);
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                    const int i = yee_idx + iR*stride_any_direction[R];
                    const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                       + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                          + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                             + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                + np->energy[ec_2][i-s_2+s_ec]);
                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  }
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                    const int i = yee_idx + iZ;
                    const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                       + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                          + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                             + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                + np->energy[ec_2][i-s_2+s_ec]);
                    np->s[ec][i] = max(-energy_here*fac, 0.0);
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                         + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                            + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                               + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                  + np->energy[ec_2][i-s_2+s_ec]);
                      np->s[ec][i] = max(-energy_here*fac, 0.0);
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                      op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                         + np->s[ec][i]*f[ec][1][i]);
                    }
                  }
                }
              }
            } else { // not fac > 0
              for (int i=0;i<ntot;i++) {
                np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
              }
              if (num_any_direction[Z]==1) {
                if (num_any_direction[R]==1) {
                  const int i = yee_idx;
                  const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                     + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                        + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                           + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                              + np->energy[ec_2][i-s_2+s_ec]);
                  np->s[ec][i] = energy_here*fac;
                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                  op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                     + np->s[ec][i]*f[ec][1][i]);
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                    const int i = yee_idx + iR*stride_any_direction[R];
                    const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                       + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                          + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                             + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                + np->energy[ec_2][i-s_2+s_ec]);
                    np->s[ec][i] = energy_here*fac;
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  }
                }
              } else { // not num_any_direction[Z]==1
                if (num_any_direction[R]==1) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                    const int i = yee_idx + iZ;
                    const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                       + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                          + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                             + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                + np->energy[ec_2][i-s_2+s_ec]);
                    np->s[ec][i] = energy_here*fac;
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  }
                } else { // not num_any_direction[R]==1
                  for (int iR=0; iR<num_any_direction[R]; iR += stride_any_direction[R]) {
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                      const int i = yee_idx + iZ + iR*stride_any_direction[R];
                      const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                         + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                            + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                               + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                  + np->energy[ec_2][i-s_2+s_ec]);
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
          } else { // not fac
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
} else { // not f[Er][0]
  if (f[Ey][0]) {
    if (f[Ez][0]) {
      if (stride_any_direction[Z]) {
        // Am in 3D
        FOR_ELECTRIC_COMPONENTS(ec) if (f[ec][0]) {
          const int yee_idx = v.yee_index(ec);
          const int s_ec = stride_any_direction[component_direction(ec)];
          const direction d_1 = (direction)((ec+1)%3);
          const component ec_1 = direction_component(ec,d_1);
          const int s_1 = stride_any_direction[d_1];
          const direction d_2 = (direction)((ec+2)%3);
          const component ec_2 = direction_component(ec,d_2);
          const int s_2 = stride_any_direction[d_2];
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
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        } else { // not num_any_direction[Y]==1
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            const int i = yee_idx + iX*stride_any_direction[X];
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                       + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                          + np->energy[ec_2][i-s_2+s_ec]);
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                            const int i = yee_idx + iZ;
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                              const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                       + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                          + np->energy[ec_2][i-s_2+s_ec]);
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                       + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                          + np->energy[ec_2][i-s_2+s_ec]);
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                   + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                      + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                         + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                            + np->energy[ec_2][i-s_2+s_ec]);
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            }
                          }
                        }
                      }
                    }
                  } else { // not fac > 0
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    }
                    if (num_any_direction[Z]==1) {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          const int i = yee_idx;
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        } else { // not num_any_direction[Y]==1
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            const int i = yee_idx + iX*stride_any_direction[X];
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                       + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                          + np->energy[ec_2][i-s_2+s_ec]);
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                            const int i = yee_idx + iZ;
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                              const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                       + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                          + np->energy[ec_2][i-s_2+s_ec]);
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                       + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                          + np->energy[ec_2][i-s_2+s_ec]);
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                   + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                      + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                         + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                            + np->energy[ec_2][i-s_2+s_ec]);
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
                } else { // not fac
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                }
              }
            }
          } else { // not is_real
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
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        } else { // not num_any_direction[Y]==1
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            const int i = yee_idx + iX*stride_any_direction[X];
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                       + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                          + np->energy[ec_2][i-s_2+s_ec]);
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                            const int i = yee_idx + iZ;
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                              const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                       + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                          + np->energy[ec_2][i-s_2+s_ec]);
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                       + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                          + np->energy[ec_2][i-s_2+s_ec]);
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                   + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                      + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                         + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                            + np->energy[ec_2][i-s_2+s_ec]);
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
                  } else { // not fac > 0
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                      np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    }
                    if (num_any_direction[Z]==1) {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          const int i = yee_idx;
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        } else { // not num_any_direction[Y]==1
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            const int i = yee_idx + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            const int i = yee_idx + iX*stride_any_direction[X];
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                              const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                       + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                          + np->energy[ec_2][i-s_2+s_ec]);
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                            const int i = yee_idx + iZ;
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                              const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                       + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                          + np->energy[ec_2][i-s_2+s_ec]);
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                              const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                       + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                          + np->energy[ec_2][i-s_2+s_ec]);
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += stride_any_direction[Y]) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                   + iZ;
                                const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                   + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                      + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                         + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                            + np->energy[ec_2][i-s_2+s_ec]);
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
                } else { // not fac
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
      } else { // not stride_any_direction[Z]
        // Am in 2D
        FOR_ELECTRIC_COMPONENTS(ec) if (f[ec][0]) {
          const int yee_idx = v.yee_index(ec);
          const int s_ec = stride_any_direction[component_direction(ec)];
          const direction d_1 = (direction)((ec+1)%3);
          const component ec_1 = direction_component(ec,d_1);
          const int s_1 = stride_any_direction[d_1];
          const direction d_2 = (direction)((ec+2)%3);
          const component ec_2 = direction_component(ec,d_2);
          const int s_2 = stride_any_direction[d_2];
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
                        const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                           + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                              + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                 + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                    + np->energy[ec_2][i-s_2+s_ec]);
                        np->s[ec][i] = max(-energy_here*fac, 0.0);
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      } else { // not num_any_direction[Y]==1
                        for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                          const int i = yee_idx + iY;
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      }
                    }
                  } else { // not fac > 0
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    }
                    if (num_any_direction[X]==1) {
                      if (num_any_direction[Y]==1) {
                        const int i = yee_idx;
                        const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                           + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                              + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                 + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                    + np->energy[ec_2][i-s_2+s_ec]);
                        np->s[ec][i] = energy_here*fac;
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      } else { // not num_any_direction[Y]==1
                        for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                          const int i = yee_idx + iY;
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      }
                    }
                  }
                } else { // not fac
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                }
              }
            }
          } else { // not is_real
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
                        const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                           + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                              + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                 + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                    + np->energy[ec_2][i-s_2+s_ec]);
                        np->s[ec][i] = max(-energy_here*fac, 0.0);
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                        op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                           + np->s[ec][i]*f[ec][1][i]);
                      } else { // not num_any_direction[Y]==1
                        for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                          const int i = yee_idx + iY;
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        }
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      }
                    }
                  } else { // not fac > 0
                    for (int i=0;i<ntot;i++) {
                      np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                      np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    }
                    if (num_any_direction[X]==1) {
                      if (num_any_direction[Y]==1) {
                        const int i = yee_idx;
                        const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                           + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                              + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                 + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                    + np->energy[ec_2][i-s_2+s_ec]);
                        np->s[ec][i] = energy_here*fac;
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                        op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                           + np->s[ec][i]*f[ec][1][i]);
                      } else { // not num_any_direction[Y]==1
                        for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                          const int i = yee_idx + iY;
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        }
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
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
                } else { // not fac
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
    } else { // not f[Ez][0]
      // Am in 2D TE
      FOR_ELECTRIC_COMPONENTS(ec) if (f[ec][0]) {
        const int yee_idx = v.yee_index(ec);
        const int s_ec = stride_any_direction[component_direction(ec)];
        const direction d_1 = (direction)((ec+1)%2);
        const component ec_1 = direction_component(ec,d_1);
        const int s_1 = stride_any_direction[d_1];
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
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          const int i = yee_idx;
                          const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec]);
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        } else { // not num_any_direction[Y]==1
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iY;
                            const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec]);
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            const int i = yee_idx + iX*stride_any_direction[X];
                            const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec]);
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY;
                              const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec]);
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                            const int i = yee_idx + iZ*stride_any_direction[Z];
                            const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec]);
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                              const int i = yee_idx + iY + iZ*stride_any_direction[Z];
                              const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec]);
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ*stride_any_direction[Z];
                              const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec]);
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z];
                                const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                                   + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                      + np->energy[ec_1][i-s_1+s_ec]);
                                np->s[ec][i] = max(-energy_here*fac, 0.0);
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            }
                          }
                        }
                      }
                    }
                  } else { // not stride_any_direction[Z]
                    if (num_any_direction[X]==1) {
                      if (num_any_direction[Y]==1) {
                        const int i = yee_idx;
                        const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                           + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                              + np->energy[ec_1][i-s_1+s_ec]);
                        np->s[ec][i] = max(-energy_here*fac, 0.0);
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      } else { // not num_any_direction[Y]==1
                        for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                          const int i = yee_idx + iY;
                          const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec]);
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec]);
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec]);
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      }
                    }
                  }
                } else { // not fac > 0
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  }
                  if (stride_any_direction[Z]) {
                    if (num_any_direction[Z]==1) {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          const int i = yee_idx;
                          const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec]);
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        } else { // not num_any_direction[Y]==1
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iY;
                            const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec]);
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            const int i = yee_idx + iX*stride_any_direction[X];
                            const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec]);
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY;
                              const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec]);
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                            const int i = yee_idx + iZ*stride_any_direction[Z];
                            const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec]);
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                              const int i = yee_idx + iY + iZ*stride_any_direction[Z];
                              const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec]);
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ*stride_any_direction[Z];
                              const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec]);
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z];
                                const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                                   + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                      + np->energy[ec_1][i-s_1+s_ec]);
                                np->s[ec][i] = energy_here*fac;
                                op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                   + np->s[ec][i]*f[ec][0][i]);
                              }
                            }
                          }
                        }
                      }
                    }
                  } else { // not stride_any_direction[Z]
                    if (num_any_direction[X]==1) {
                      if (num_any_direction[Y]==1) {
                        const int i = yee_idx;
                        const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                           + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                              + np->energy[ec_1][i-s_1+s_ec]);
                        np->s[ec][i] = energy_here*fac;
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      } else { // not num_any_direction[Y]==1
                        for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                          const int i = yee_idx + iY;
                          const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec]);
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec]);
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec]);
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      }
                    }
                  }
                }
              } else { // not fac
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                }
              }
            }
          }
        } else { // not is_real
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
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          const int i = yee_idx;
                          const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec]);
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        } else { // not num_any_direction[Y]==1
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iY;
                            const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec]);
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            const int i = yee_idx + iX*stride_any_direction[X];
                            const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec]);
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY;
                              const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec]);
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                            const int i = yee_idx + iZ*stride_any_direction[Z];
                            const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec]);
                            np->s[ec][i] = max(-energy_here*fac, 0.0);
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                              const int i = yee_idx + iY + iZ*stride_any_direction[Z];
                              const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec]);
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ*stride_any_direction[Z];
                              const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec]);
                              np->s[ec][i] = max(-energy_here*fac, 0.0);
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z];
                                const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                                   + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                      + np->energy[ec_1][i-s_1+s_ec]);
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
                  } else { // not stride_any_direction[Z]
                    if (num_any_direction[X]==1) {
                      if (num_any_direction[Y]==1) {
                        const int i = yee_idx;
                        const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                           + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                              + np->energy[ec_1][i-s_1+s_ec]);
                        np->s[ec][i] = max(-energy_here*fac, 0.0);
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                        op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                           + np->s[ec][i]*f[ec][1][i]);
                      } else { // not num_any_direction[Y]==1
                        for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                          const int i = yee_idx + iY;
                          const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec]);
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        }
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec]);
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec]);
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
                } else { // not fac > 0
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                  }
                  if (stride_any_direction[Z]) {
                    if (num_any_direction[Z]==1) {
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          const int i = yee_idx;
                          const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec]);
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        } else { // not num_any_direction[Y]==1
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iY;
                            const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec]);
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            const int i = yee_idx + iX*stride_any_direction[X];
                            const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec]);
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY;
                              const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec]);
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    } else { // not num_any_direction[Z]==1
                      if (num_any_direction[X]==1) {
                        if (num_any_direction[Y]==1) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                            const int i = yee_idx + iZ*stride_any_direction[Z];
                            const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec]);
                            np->s[ec][i] = energy_here*fac;
                            op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                              const int i = yee_idx + iY + iZ*stride_any_direction[Z];
                              const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec]);
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      } else { // not num_any_direction[X]==1
                        if (num_any_direction[Y]==1) {
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iZ*stride_any_direction[Z];
                              const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec]);
                              np->s[ec][i] = energy_here*fac;
                              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        } else { // not num_any_direction[Y]==1
                          for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                            for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                              for (int iZ=0; iZ<num_any_direction[Z]; iZ += stride_any_direction[Z]) {
                                const int i = yee_idx + iX*stride_any_direction[X] + iY + iZ*stride_any_direction[Z];
                                const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                                   + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                      + np->energy[ec_1][i-s_1+s_ec]);
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
                  } else { // not stride_any_direction[Z]
                    if (num_any_direction[X]==1) {
                      if (num_any_direction[Y]==1) {
                        const int i = yee_idx;
                        const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                           + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                              + np->energy[ec_1][i-s_1+s_ec]);
                        np->s[ec][i] = energy_here*fac;
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                        op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                           + np->s[ec][i]*f[ec][1][i]);
                      } else { // not num_any_direction[Y]==1
                        for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                          const int i = yee_idx + iY;
                          const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec]);
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        }
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec]);
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                          for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY;
                            const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec]);
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
              } else { // not fac
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
  } else { // not f[Ey][0]
    if (f[Ex][0]) {
      // Am in 1D
      FOR_ELECTRIC_COMPONENTS(ec) if (f[ec][0]) {
        const int yee_idx = v.yee_index(ec);
        const int s_ec = stride_any_direction[component_direction(ec)];
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
                  } else { // not num_any_direction[Z]==1
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                      const int i = yee_idx + iZ;
                      const double energy_here = np->energy[Ex][i];
                      np->s[ec][i] = max(-energy_here*fac, 0.0);
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    }
                  }
                } else { // not fac > 0
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  }
                  if (num_any_direction[Z]==1) {
                    const int i = yee_idx;
                    const double energy_here = np->energy[Ex][i];
                    np->s[ec][i] = energy_here*fac;
                    op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  } else { // not num_any_direction[Z]==1
                    for (int iZ=0; iZ<num_any_direction[Z]; iZ += 1) {
                      const int i = yee_idx + iZ;
                      const double energy_here = np->energy[Ex][i];
                      np->s[ec][i] = energy_here*fac;
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    }
                  }
                }
              } else { // not fac
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                }
              }
            }
          }
        } else { // not is_real
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
                  } else { // not num_any_direction[Z]==1
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
                } else { // not fac > 0
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
                  } else { // not num_any_direction[Z]==1
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
              } else { // not fac
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
    } else { // not f[Ex][0]
      // Am in 2D TM
      FOR_ELECTRIC_COMPONENTS(ec) if (f[ec][0]) {
        const int yee_idx = v.yee_index(ec);
        const int s_ec = stride_any_direction[component_direction(ec)];
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
                      const double energy_here = np->energy[Ez][i];
                      np->s[ec][i] = max(-energy_here*fac, 0.0);
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    } else { // not num_any_direction[Y]==1
                      for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                        const int i = yee_idx + iY;
                        const double energy_here = np->energy[Ez][i];
                        np->s[ec][i] = max(-energy_here*fac, 0.0);
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      }
                    }
                  } else { // not num_any_direction[X]==1
                    if (num_any_direction[Y]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                        const int i = yee_idx + iX*stride_any_direction[X];
                        const double energy_here = np->energy[Ez][i];
                        np->s[ec][i] = max(-energy_here*fac, 0.0);
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      }
                    } else { // not num_any_direction[Y]==1
                      for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                        for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY;
                          const double energy_here = np->energy[Ez][i];
                          np->s[ec][i] = max(-energy_here*fac, 0.0);
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      }
                    }
                  }
                } else { // not fac > 0
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  }
                  if (num_any_direction[X]==1) {
                    if (num_any_direction[Y]==1) {
                      const int i = yee_idx;
                      const double energy_here = np->energy[Ez][i];
                      np->s[ec][i] = energy_here*fac;
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    } else { // not num_any_direction[Y]==1
                      for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                        const int i = yee_idx + iY;
                        const double energy_here = np->energy[Ez][i];
                        np->s[ec][i] = energy_here*fac;
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      }
                    }
                  } else { // not num_any_direction[X]==1
                    if (num_any_direction[Y]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                        const int i = yee_idx + iX*stride_any_direction[X];
                        const double energy_here = np->energy[Ez][i];
                        np->s[ec][i] = energy_here*fac;
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      }
                    } else { // not num_any_direction[Y]==1
                      for (int iX=0; iX<num_any_direction[X]; iX += stride_any_direction[X]) {
                        for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY;
                          const double energy_here = np->energy[Ez][i];
                          np->s[ec][i] = energy_here*fac;
                          op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      }
                    }
                  }
                }
              } else { // not fac
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                }
              }
            }
          }
        } else { // not is_real
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
                      const double energy_here = np->energy[Ez][i];
                      np->s[ec][i] = max(-energy_here*fac, 0.0);
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                      op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                         + np->s[ec][i]*f[ec][1][i]);
                    } else { // not num_any_direction[Y]==1
                      for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                        const int i = yee_idx + iY;
                        const double energy_here = np->energy[Ez][i];
                        np->s[ec][i] = max(-energy_here*fac, 0.0);
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                        op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                           + np->s[ec][i]*f[ec][1][i]);
                      }
                    }
                  } else { // not num_any_direction[X]==1
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
                    } else { // not num_any_direction[Y]==1
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
                    }
                  }
                } else { // not fac > 0
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                  }
                  if (num_any_direction[X]==1) {
                    if (num_any_direction[Y]==1) {
                      const int i = yee_idx;
                      const double energy_here = np->energy[Ez][i];
                      np->s[ec][i] = energy_here*fac;
                      op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                      op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                         + np->s[ec][i]*f[ec][1][i]);
                    } else { // not num_any_direction[Y]==1
                      for (int iY=0; iY<num_any_direction[Y]; iY += 1) {
                        const int i = yee_idx + iY;
                        const double energy_here = np->energy[Ez][i];
                        np->s[ec][i] = energy_here*fac;
                        op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                        op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                           + np->s[ec][i]*f[ec][1][i]);
                      }
                    }
                  } else { // not num_any_direction[X]==1
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
                    } else { // not num_any_direction[Y]==1
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
                    }
                  }
                }
              } else { // not fac
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
