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
      if (pol) {
        for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
          const double om_psqr = fabs(op->pb->saturated_sigma);
          const double oo_ep_om_psqr = 1.0/(8*pi*om_psqr);
          const double oo_e_sat = 1.0/op->pb->energy_saturation;
          const double om_sqr = op->pb->omeganot*op->pb->omeganot;
          const double g = op->pb->gamma;
          const double funinv = 1.0/(1+0.5*g);
          if (np->saturation_factor) {
            for (int i=0;i<ntot;i++) {
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
            }
            if (num_any_direction[Z]==1) {
              if (num_any_direction[R]==1) {
                const int i = yee_idx;
                const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                   + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                      + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                         + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                            + np->energy[ec_2][i-s_2+s_ec]);
                const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                   + op->P[ec][0][i])) + (0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                      + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                         + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                            op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                               + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                  + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*((((np->P[ec_2][0][i])
                                     + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i]))
                                        + (((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2])
                                           + op->P[ec_2][0][i-s_2])) + (((np->P[ec_2][0][i+s_ec]) +
                                              op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec]))
                                                 + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec])
                                                    + op->P[ec_2][0][i-s_2+s_ec]));
                const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                   - (op->P[ec][0][i]))) + (0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                      - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                         - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                            - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*((((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + (((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + (((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                   + np->s[ec][i]*f[ec][0][i]);
              } else { // not num_any_direction[R]==1
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                     + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                        + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                           + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                              + np->energy[ec_2][i-s_2+s_ec]);
                  const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                     + op->P[ec][0][i])) + (0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                        + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                           + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                              op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                 + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                    + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*((((np->P[ec_2][0][i])
                                       + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i]))
                                          + (((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2])
                                             + op->P[ec_2][0][i-s_2])) + (((np->P[ec_2][0][i+s_ec]) +
                                                op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec]))
                                                   + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec])
                                                      + op->P[ec_2][0][i-s_2+s_ec]));
                  const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                     - (op->P[ec][0][i]))) + (0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                        - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                           - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                              - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*((((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + (((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + (((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                  const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                  np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                  op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                }
              }
            } else { // not num_any_direction[Z]==1
              if (num_any_direction[R]==1) {
                for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                  const int i = yee_idx + iZ;
                  const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                     + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                        + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                           + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                              + np->energy[ec_2][i-s_2+s_ec]);
                  const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                     + op->P[ec][0][i])) + (0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                        + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                           + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                              op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                 + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                    + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*((((np->P[ec_2][0][i])
                                       + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i]))
                                          + (((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2])
                                             + op->P[ec_2][0][i-s_2])) + (((np->P[ec_2][0][i+s_ec]) +
                                                op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec]))
                                                   + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec])
                                                      + op->P[ec_2][0][i-s_2+s_ec]));
                  const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                     - (op->P[ec][0][i]))) + (0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                        - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                           - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                              - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*((((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + (((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + (((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                  const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                  np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                  op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                }
              } else { // not num_any_direction[R]==1
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ + iR*stride_any_direction[R];
                    const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                       + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                          + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                             + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                + np->energy[ec_2][i-s_2+s_ec]);
                    const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                       + op->P[ec][0][i])) + (0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                          + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                             + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                   + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                      + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*((((np->P[ec_2][0][i])
                                         + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i]))
                                            + (((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2])
                                               + op->P[ec_2][0][i-s_2])) + (((np->P[ec_2][0][i+s_ec]) +
                                                  op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec]))
                                                     + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec])
                                                        + op->P[ec_2][0][i-s_2+s_ec]));
                    const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                       - (op->P[ec][0][i]))) + (0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                          - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                             - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*((((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + (((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + (((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                    const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                    np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                    op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                }
              }
            }
          } else { // not np->saturation_factor
            for (int i=0;i<ntot;i++) {
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
              op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                 + np->s[ec][i]*f[ec][0][i]);
            }
          }
        }
      }
    } else { // not is_real
      if (pol) {
        for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
          const double om_psqr = fabs(op->pb->saturated_sigma);
          const double oo_ep_om_psqr = 1.0/(8*pi*om_psqr);
          const double oo_e_sat = 1.0/op->pb->energy_saturation;
          const double om_sqr = op->pb->omeganot*op->pb->omeganot;
          const double g = op->pb->gamma;
          const double funinv = 1.0/(1+0.5*g);
          if (np->saturation_factor) {
            for (int i=0;i<ntot;i++) {
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
            }
            if (num_any_direction[Z]==1) {
              if (num_any_direction[R]==1) {
                const int i = yee_idx;
                const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                   + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                      + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                         + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                            + np->energy[ec_2][i-s_2+s_ec]);
                const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                   + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                      + op->P[ec][0][i])) + (0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                         + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                            + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                               + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                  + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*(((((np->P[ec_2][1][i]) + op->P[ec_2][1][i])*((np->P[ec_2][1][i]) + op->P[ec_2][1][i])) + ((np->P[ec_2][0][i]) + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i])) + ((((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])*((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])) + ((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])) + ((((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])*((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])) + ((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])) + (((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])*((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])) + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec]));
                const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                   - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                      - (op->P[ec][0][i]))) + (0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))) + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*(((((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))*((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))) + ((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + ((((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))*((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))) + ((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + ((((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))*((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))) + ((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + (((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))*((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                   + np->s[ec][i]*f[ec][0][i]);
                op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                   + np->s[ec][i]*f[ec][1][i]);
              } else { // not num_any_direction[R]==1
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  const int i = yee_idx + iR*stride_any_direction[R];
                  const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                     + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                        + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                           + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                              + np->energy[ec_2][i-s_2+s_ec]);
                  const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                     + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                        + op->P[ec][0][i])) + (0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                           + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                              + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                 + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                    + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*(((((np->P[ec_2][1][i]) + op->P[ec_2][1][i])*((np->P[ec_2][1][i]) + op->P[ec_2][1][i])) + ((np->P[ec_2][0][i]) + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i])) + ((((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])*((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])) + ((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])) + ((((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])*((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])) + ((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])) + (((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])*((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])) + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec]));
                  const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                     - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                        - (op->P[ec][0][i]))) + (0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))) + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*(((((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))*((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))) + ((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + ((((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))*((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))) + ((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + ((((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))*((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))) + ((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + (((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))*((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                  const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                  np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                  op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                  op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                     + np->s[ec][i]*f[ec][1][i]);
                }
              }
            } else { // not num_any_direction[Z]==1
              if (num_any_direction[R]==1) {
                for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                  const int i = yee_idx + iZ;
                  const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                     + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                        + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                           + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                              + np->energy[ec_2][i-s_2+s_ec]);
                  const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                     + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                        + op->P[ec][0][i])) + (0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                           + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                              + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                 + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                    + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*(((((np->P[ec_2][1][i]) + op->P[ec_2][1][i])*((np->P[ec_2][1][i]) + op->P[ec_2][1][i])) + ((np->P[ec_2][0][i]) + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i])) + ((((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])*((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])) + ((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])) + ((((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])*((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])) + ((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])) + (((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])*((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])) + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec]));
                  const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                     - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                        - (op->P[ec][0][i]))) + (0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))) + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*(((((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))*((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))) + ((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + ((((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))*((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))) + ((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + ((((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))*((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))) + ((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + (((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))*((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                  const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                  np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                  op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                  op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                     + np->s[ec][i]*f[ec][1][i]);
                }
              } else { // not num_any_direction[R]==1
                for (int iR=0; iR<num_any_direction[R]; iR++) {
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ + iR*stride_any_direction[R];
                    const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                       + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                          + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                             + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                + np->energy[ec_2][i-s_2+s_ec]);
                    const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                       + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                          + op->P[ec][0][i])) + (0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                             + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                   + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                      + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*(((((np->P[ec_2][1][i]) + op->P[ec_2][1][i])*((np->P[ec_2][1][i]) + op->P[ec_2][1][i])) + ((np->P[ec_2][0][i]) + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i])) + ((((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])*((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])) + ((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])) + ((((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])*((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])) + ((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])) + (((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])*((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])) + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec]));
                    const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                       - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                          - (op->P[ec][0][i]))) + (0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))) + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*(((((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))*((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))) + ((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + ((((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))*((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))) + ((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + ((((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))*((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))) + ((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + (((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))*((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                    const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                    np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                    op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  }
                }
              }
            }
          } else { // not np->saturation_factor
            for (int i=0;i<ntot;i++) {
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
              np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
              op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                 + np->s[ec][i]*f[ec][0][i]);
              op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
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
            if (pol) {
              for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                const double om_psqr = fabs(op->pb->saturated_sigma);
                const double oo_ep_om_psqr = 1.0/(8*pi*om_psqr);
                const double oo_e_sat = 1.0/op->pb->energy_saturation;
                const double om_sqr = op->pb->omeganot*op->pb->omeganot;
                const double g = op->pb->gamma;
                const double funinv = 1.0/(1+0.5*g);
                if (np->saturation_factor) {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
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
                        const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                           + op->P[ec][0][i])) + (0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                              + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                 + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                    op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                       + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                          + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*((((np->P[ec_2][0][i])
                                             + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i]))
                                                + (((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2])
                                                   + op->P[ec_2][0][i-s_2])) + (((np->P[ec_2][0][i+s_ec]) +
                                                      op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec]))
                                                         + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec])
                                                            + op->P[ec_2][0][i-s_2+s_ec]));
                        const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                           - (op->P[ec][0][i]))) + (0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                              - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                                 - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                    - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*((((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + (((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + (((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                        const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                        np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                        op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      } else { // not num_any_direction[Y]==1
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iY*stride_any_direction[Y];
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                             + op->P[ec][0][i])) + (0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                   + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                      op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                         + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                            + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*((((np->P[ec_2][0][i])
                                               + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i]))
                                                  + (((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2])
                                                     + op->P[ec_2][0][i-s_2])) + (((np->P[ec_2][0][i+s_ec]) +
                                                        op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec]))
                                                           + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec])
                                                              + op->P[ec_2][0][i-s_2+s_ec]));
                          const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                             - (op->P[ec][0][i]))) + (0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                                - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                                   - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                      - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*((((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + (((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + (((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                          const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                          np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                          op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                             + op->P[ec][0][i])) + (0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                   + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                      op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                         + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                            + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*((((np->P[ec_2][0][i])
                                               + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i]))
                                                  + (((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2])
                                                     + op->P[ec_2][0][i-s_2])) + (((np->P[ec_2][0][i+s_ec]) +
                                                        op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec]))
                                                           + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec])
                                                              + op->P[ec_2][0][i-s_2+s_ec]));
                          const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                             - (op->P[ec][0][i]))) + (0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                                - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                                   - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                      - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*((((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + (((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + (((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                          const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                          np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                          op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                               + op->P[ec][0][i])) + (0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                  + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                     + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                        op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                           + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                              + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*((((np->P[ec_2][0][i])
                                                 + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i]))
                                                    + (((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2])
                                                       + op->P[ec_2][0][i-s_2])) + (((np->P[ec_2][0][i+s_ec]) +
                                                          op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec]))
                                                             + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec])
                                                                + op->P[ec_2][0][i-s_2+s_ec]));
                            const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                               - (op->P[ec][0][i]))) + (0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                                  - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                                     - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                        - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*((((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + (((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + (((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                            const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                            np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                            op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      }
                    }
                  } else { // not num_any_direction[Z]==1
                    if (num_any_direction[X]==1) {
                      if (num_any_direction[Y]==1) {
                        for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                          const int i = yee_idx + iZ;
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                             + op->P[ec][0][i])) + (0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                   + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                      op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                         + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                            + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*((((np->P[ec_2][0][i])
                                               + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i]))
                                                  + (((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2])
                                                     + op->P[ec_2][0][i-s_2])) + (((np->P[ec_2][0][i+s_ec]) +
                                                        op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec]))
                                                           + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec])
                                                              + op->P[ec_2][0][i-s_2+s_ec]));
                          const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                             - (op->P[ec][0][i]))) + (0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                                - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                                   - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                      - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*((((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + (((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + (((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                          const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                          np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                          op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                               + op->P[ec][0][i])) + (0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                  + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                     + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                        op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                           + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                              + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*((((np->P[ec_2][0][i])
                                                 + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i]))
                                                    + (((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2])
                                                       + op->P[ec_2][0][i-s_2])) + (((np->P[ec_2][0][i+s_ec]) +
                                                          op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec]))
                                                             + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec])
                                                                + op->P[ec_2][0][i-s_2+s_ec]));
                            const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                               - (op->P[ec][0][i]))) + (0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                                  - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                                     - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                        - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*((((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + (((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + (((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                            const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                            np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                            op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                               + op->P[ec][0][i])) + (0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                  + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                     + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                        op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                           + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                              + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*((((np->P[ec_2][0][i])
                                                 + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i]))
                                                    + (((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2])
                                                       + op->P[ec_2][0][i-s_2])) + (((np->P[ec_2][0][i+s_ec]) +
                                                          op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec]))
                                                             + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec])
                                                                + op->P[ec_2][0][i-s_2+s_ec]));
                            const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                               - (op->P[ec][0][i]))) + (0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                                  - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                                     - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                        - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*((((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + (((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + (((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                            const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                            np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                            op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                          }
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                 + iZ;
                              const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                       + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                          + np->energy[ec_2][i-s_2+s_ec]);
                              const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                                 + op->P[ec][0][i])) + (0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                    + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                       + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                          op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                             + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                                + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*((((np->P[ec_2][0][i])
                                                   + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i]))
                                                      + (((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2])
                                                         + op->P[ec_2][0][i-s_2])) + (((np->P[ec_2][0][i+s_ec]) +
                                                            op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec]))
                                                               + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec])
                                                                  + op->P[ec_2][0][i-s_2+s_ec]));
                              const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                                 - (op->P[ec][0][i]))) + (0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                                    - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                                       - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                          - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*((((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + (((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + (((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                              const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                              np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                              op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                            }
                          }
                        }
                      }
                    }
                  }
                } else { // not np->saturation_factor
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                }
              }
            }
          } else { // not is_real
            if (pol) {
              for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                const double om_psqr = fabs(op->pb->saturated_sigma);
                const double oo_ep_om_psqr = 1.0/(8*pi*om_psqr);
                const double oo_e_sat = 1.0/op->pb->energy_saturation;
                const double om_sqr = op->pb->omeganot*op->pb->omeganot;
                const double g = op->pb->gamma;
                const double funinv = 1.0/(1+0.5*g);
                if (np->saturation_factor) {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
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
                        const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                           + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                              + op->P[ec][0][i])) + (0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                                 + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                    + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                       + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                          + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*(((((np->P[ec_2][1][i]) + op->P[ec_2][1][i])*((np->P[ec_2][1][i]) + op->P[ec_2][1][i])) + ((np->P[ec_2][0][i]) + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i])) + ((((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])*((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])) + ((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])) + ((((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])*((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])) + ((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])) + (((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])*((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])) + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec]));
                        const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                           - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                              - (op->P[ec][0][i]))) + (0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))) + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*(((((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))*((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))) + ((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + ((((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))*((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))) + ((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + ((((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))*((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))) + ((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + (((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))*((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                        const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                        np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                        op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                        op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                           + np->s[ec][i]*f[ec][1][i]);
                      } else { // not num_any_direction[Y]==1
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iY*stride_any_direction[Y];
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                             + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                                + op->P[ec][0][i])) + (0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                                   + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                      + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                         + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                            + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*(((((np->P[ec_2][1][i]) + op->P[ec_2][1][i])*((np->P[ec_2][1][i]) + op->P[ec_2][1][i])) + ((np->P[ec_2][0][i]) + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i])) + ((((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])*((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])) + ((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])) + ((((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])*((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])) + ((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])) + (((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])*((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])) + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec]));
                          const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                             - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                                - (op->P[ec][0][i]))) + (0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))) + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*(((((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))*((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))) + ((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + ((((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))*((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))) + ((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + ((((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))*((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))) + ((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + (((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))*((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                          const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                          np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                          op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        }
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          const int i = yee_idx + iX*stride_any_direction[X];
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                             + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                                + op->P[ec][0][i])) + (0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                                   + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                      + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                         + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                            + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*(((((np->P[ec_2][1][i]) + op->P[ec_2][1][i])*((np->P[ec_2][1][i]) + op->P[ec_2][1][i])) + ((np->P[ec_2][0][i]) + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i])) + ((((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])*((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])) + ((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])) + ((((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])*((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])) + ((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])) + (((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])*((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])) + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec]));
                          const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                             - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                                - (op->P[ec][0][i]))) + (0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))) + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*(((((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))*((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))) + ((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + ((((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))*((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))) + ((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + ((((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))*((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))) + ((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + (((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))*((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                          const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                          np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                          op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y];
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                               + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                                  + op->P[ec][0][i])) + (0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                                     + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                        + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                           + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                              + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*(((((np->P[ec_2][1][i]) + op->P[ec_2][1][i])*((np->P[ec_2][1][i]) + op->P[ec_2][1][i])) + ((np->P[ec_2][0][i]) + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i])) + ((((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])*((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])) + ((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])) + ((((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])*((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])) + ((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])) + (((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])*((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])) + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec]));
                            const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                               - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                                  - (op->P[ec][0][i]))) + (0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))) + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*(((((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))*((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))) + ((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + ((((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))*((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))) + ((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + ((((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))*((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))) + ((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + (((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))*((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                            const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                            np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                            op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      }
                    }
                  } else { // not num_any_direction[Z]==1
                    if (num_any_direction[X]==1) {
                      if (num_any_direction[Y]==1) {
                        for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                          const int i = yee_idx + iZ;
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                             + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                                + op->P[ec][0][i])) + (0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                                   + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                      + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                         + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                            + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*(((((np->P[ec_2][1][i]) + op->P[ec_2][1][i])*((np->P[ec_2][1][i]) + op->P[ec_2][1][i])) + ((np->P[ec_2][0][i]) + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i])) + ((((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])*((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])) + ((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])) + ((((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])*((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])) + ((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])) + (((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])*((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])) + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec]));
                          const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                             - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                                - (op->P[ec][0][i]))) + (0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))) + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*(((((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))*((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))) + ((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + ((((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))*((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))) + ((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + ((((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))*((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))) + ((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + (((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))*((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                          const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                          np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                          op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iY*stride_any_direction[Y] + iZ;
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                               + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                                  + op->P[ec][0][i])) + (0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                                     + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                        + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                           + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                              + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*(((((np->P[ec_2][1][i]) + op->P[ec_2][1][i])*((np->P[ec_2][1][i]) + op->P[ec_2][1][i])) + ((np->P[ec_2][0][i]) + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i])) + ((((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])*((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])) + ((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])) + ((((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])*((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])) + ((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])) + (((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])*((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])) + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec]));
                            const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                               - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                                  - (op->P[ec][0][i]))) + (0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))) + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*(((((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))*((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))) + ((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + ((((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))*((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))) + ((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + ((((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))*((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))) + ((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + (((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))*((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                            const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                            np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                            op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      }
                    } else { // not num_any_direction[X]==1
                      if (num_any_direction[Y]==1) {
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                            const int i = yee_idx + iX*stride_any_direction[X] + iZ;
                            const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                               + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                  + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                     + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                        + np->energy[ec_2][i-s_2+s_ec]);
                            const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                               + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                                  + op->P[ec][0][i])) + (0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                                     + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                        + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                           + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                              + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*(((((np->P[ec_2][1][i]) + op->P[ec_2][1][i])*((np->P[ec_2][1][i]) + op->P[ec_2][1][i])) + ((np->P[ec_2][0][i]) + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i])) + ((((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])*((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])) + ((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])) + ((((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])*((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])) + ((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])) + (((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])*((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])) + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec]));
                            const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                               - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                                  - (op->P[ec][0][i]))) + (0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))) + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*(((((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))*((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))) + ((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + ((((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))*((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))) + ((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + ((((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))*((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))) + ((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + (((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))*((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                            const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                            np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                            op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                               + np->s[ec][i]*f[ec][0][i]);
                            op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                               + np->s[ec][i]*f[ec][1][i]);
                          }
                        }
                      } else { // not num_any_direction[Y]==1
                        for (int iX=0; iX<num_any_direction[X]; iX++) {
                          for (int iY=0; iY<num_any_direction[Y]; iY++) {
                            for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                              const int i = yee_idx + iX*stride_any_direction[X] + iY*stride_any_direction[Y]
                                 + iZ;
                              const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                                 + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                    + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                       + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                          + np->energy[ec_2][i-s_2+s_ec]);
                              const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                                 + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                                    + op->P[ec][0][i])) + (0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                                       + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                          + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                             + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                                + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*(((((np->P[ec_2][1][i]) + op->P[ec_2][1][i])*((np->P[ec_2][1][i]) + op->P[ec_2][1][i])) + ((np->P[ec_2][0][i]) + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i])) + ((((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])*((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])) + ((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])) + ((((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])*((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])) + ((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])) + (((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])*((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])) + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec]));
                              const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                                 - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                                    - (op->P[ec][0][i]))) + (0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))) + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*(((((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))*((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))) + ((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + ((((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))*((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))) + ((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + ((((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))*((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))) + ((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + (((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))*((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                              const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                              np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                              op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                                 + np->s[ec][i]*f[ec][0][i]);
                              op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                                 + np->s[ec][i]*f[ec][1][i]);
                            }
                          }
                        }
                      }
                    }
                  }
                } else { // not np->saturation_factor
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
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
            if (pol) {
              for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                const double om_psqr = fabs(op->pb->saturated_sigma);
                const double oo_ep_om_psqr = 1.0/(8*pi*om_psqr);
                const double oo_e_sat = 1.0/op->pb->energy_saturation;
                const double om_sqr = op->pb->omeganot*op->pb->omeganot;
                const double g = op->pb->gamma;
                const double funinv = 1.0/(1+0.5*g);
                if (np->saturation_factor) {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  }
                  if (num_any_direction[X]==1) {
                    if (num_any_direction[Y]==1) {
                      const int i = yee_idx;
                      const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                         + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                            + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                               + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                  + np->energy[ec_2][i-s_2+s_ec]);
                      const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                         + op->P[ec][0][i])) + (0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                            + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                               + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                  op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                     + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                        + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*((((np->P[ec_2][0][i])
                                           + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i]))
                                              + (((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2])
                                                 + op->P[ec_2][0][i-s_2])) + (((np->P[ec_2][0][i+s_ec]) +
                                                    op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec]))
                                                       + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec])
                                                          + op->P[ec_2][0][i-s_2+s_ec]));
                      const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                         - (op->P[ec][0][i]))) + (0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                            - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                               - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                  - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*((((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + (((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + (((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                      const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                      np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                      op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    } else { // not num_any_direction[Y]==1
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                           + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                              + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                 + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                    + np->energy[ec_2][i-s_2+s_ec]);
                        const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                           + op->P[ec][0][i])) + (0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                              + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                 + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                    op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                       + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                          + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*((((np->P[ec_2][0][i])
                                             + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i]))
                                                + (((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2])
                                                   + op->P[ec_2][0][i-s_2])) + (((np->P[ec_2][0][i+s_ec]) +
                                                      op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec]))
                                                         + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec])
                                                            + op->P[ec_2][0][i-s_2+s_ec]));
                        const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                           - (op->P[ec][0][i]))) + (0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                              - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                                 - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                    - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*((((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + (((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + (((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                        const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                        np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                        op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      }
                    }
                  } else { // not num_any_direction[X]==1
                    if (num_any_direction[Y]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        const int i = yee_idx + iX*stride_any_direction[X];
                        const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                           + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                              + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                 + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                    + np->energy[ec_2][i-s_2+s_ec]);
                        const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                           + op->P[ec][0][i])) + (0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                              + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                 + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                    op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                       + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                          + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*((((np->P[ec_2][0][i])
                                             + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i]))
                                                + (((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2])
                                                   + op->P[ec_2][0][i-s_2])) + (((np->P[ec_2][0][i+s_ec]) +
                                                      op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec]))
                                                         + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec])
                                                            + op->P[ec_2][0][i-s_2+s_ec]));
                        const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                           - (op->P[ec][0][i]))) + (0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                              - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                                 - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                    - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*((((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + (((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + (((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                        const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                        np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                        op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      }
                    } else { // not num_any_direction[Y]==1
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY;
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                             + op->P[ec][0][i])) + (0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                   + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                      op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                         + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                            + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*((((np->P[ec_2][0][i])
                                               + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i]))
                                                  + (((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2])
                                                     + op->P[ec_2][0][i-s_2])) + (((np->P[ec_2][0][i+s_ec]) +
                                                        op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec]))
                                                           + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec])
                                                              + op->P[ec_2][0][i-s_2+s_ec]));
                          const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                             - (op->P[ec][0][i]))) + (0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                                - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                                   - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                      - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*((((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + (((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + (((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                          const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                          np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                          op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                        }
                      }
                    }
                  }
                } else { // not np->saturation_factor
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                }
              }
            }
          } else { // not is_real
            if (pol) {
              for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
                const double om_psqr = fabs(op->pb->saturated_sigma);
                const double oo_ep_om_psqr = 1.0/(8*pi*om_psqr);
                const double oo_e_sat = 1.0/op->pb->energy_saturation;
                const double om_sqr = op->pb->omeganot*op->pb->omeganot;
                const double g = op->pb->gamma;
                const double funinv = 1.0/(1+0.5*g);
                if (np->saturation_factor) {
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                  }
                  if (num_any_direction[X]==1) {
                    if (num_any_direction[Y]==1) {
                      const int i = yee_idx;
                      const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                         + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                            + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                               + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                  + np->energy[ec_2][i-s_2+s_ec]);
                      const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                         + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                            + op->P[ec][0][i])) + (0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                               + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                  + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                     + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                        + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*(((((np->P[ec_2][1][i]) + op->P[ec_2][1][i])*((np->P[ec_2][1][i]) + op->P[ec_2][1][i])) + ((np->P[ec_2][0][i]) + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i])) + ((((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])*((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])) + ((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])) + ((((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])*((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])) + ((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])) + (((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])*((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])) + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec]));
                      const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                         - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                            - (op->P[ec][0][i]))) + (0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))) + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*(((((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))*((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))) + ((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + ((((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))*((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))) + ((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + ((((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))*((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))) + ((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + (((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))*((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                      const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                      np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                      op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                      op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                         + np->s[ec][i]*f[ec][1][i]);
                    } else { // not num_any_direction[Y]==1
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iY;
                        const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                           + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                              + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                 + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                    + np->energy[ec_2][i-s_2+s_ec]);
                        const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                           + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                              + op->P[ec][0][i])) + (0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                                 + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                    + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                       + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                          + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*(((((np->P[ec_2][1][i]) + op->P[ec_2][1][i])*((np->P[ec_2][1][i]) + op->P[ec_2][1][i])) + ((np->P[ec_2][0][i]) + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i])) + ((((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])*((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])) + ((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])) + ((((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])*((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])) + ((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])) + (((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])*((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])) + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec]));
                        const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                           - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                              - (op->P[ec][0][i]))) + (0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))) + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*(((((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))*((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))) + ((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + ((((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))*((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))) + ((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + ((((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))*((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))) + ((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + (((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))*((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                        const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                        np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                        op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                        op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                           + np->s[ec][i]*f[ec][1][i]);
                      }
                    }
                  } else { // not num_any_direction[X]==1
                    if (num_any_direction[Y]==1) {
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        const int i = yee_idx + iX*stride_any_direction[X];
                        const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                           + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                              + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                 + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                    + np->energy[ec_2][i-s_2+s_ec]);
                        const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                           + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                              + op->P[ec][0][i])) + (0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                                 + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                    + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                       + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                          + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*(((((np->P[ec_2][1][i]) + op->P[ec_2][1][i])*((np->P[ec_2][1][i]) + op->P[ec_2][1][i])) + ((np->P[ec_2][0][i]) + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i])) + ((((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])*((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])) + ((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])) + ((((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])*((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])) + ((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])) + (((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])*((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])) + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec]));
                        const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                           - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                              - (op->P[ec][0][i]))) + (0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))) + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*(((((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))*((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))) + ((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + ((((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))*((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))) + ((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + ((((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))*((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))) + ((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + (((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))*((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                        const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                        np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                        op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                        op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                           + np->s[ec][i]*f[ec][1][i]);
                      }
                    } else { // not num_any_direction[Y]==1
                      for (int iX=0; iX<num_any_direction[X]; iX++) {
                        for (int iY=0; iY<num_any_direction[Y]; iY++) {
                          const int i = yee_idx + iX*stride_any_direction[X] + iY;
                          const double energy_here = (np->energy[ec][i]) + (0.25*((np->energy[ec_1][i])
                             + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                                + np->energy[ec_1][i-s_1+s_ec])) + 0.25*((np->energy[ec_2][i])
                                   + (np->energy[ec_2][i-s_2]) + (np->energy[ec_2][i+s_ec])
                                      + np->energy[ec_2][i-s_2+s_ec]);
                          const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                             + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                                + op->P[ec][0][i])) + (0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                                   + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                      + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                         + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                            + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]))) + 0.25*(((((np->P[ec_2][1][i]) + op->P[ec_2][1][i])*((np->P[ec_2][1][i]) + op->P[ec_2][1][i])) + ((np->P[ec_2][0][i]) + op->P[ec_2][0][i])*((np->P[ec_2][0][i]) + op->P[ec_2][0][i])) + ((((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])*((np->P[ec_2][1][i-s_2]) + op->P[ec_2][1][i-s_2])) + ((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])*((np->P[ec_2][0][i-s_2]) + op->P[ec_2][0][i-s_2])) + ((((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])*((np->P[ec_2][1][i+s_ec]) + op->P[ec_2][1][i+s_ec])) + ((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])*((np->P[ec_2][0][i+s_ec]) + op->P[ec_2][0][i+s_ec])) + (((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])*((np->P[ec_2][1][i-s_2+s_ec]) + op->P[ec_2][1][i-s_2+s_ec])) + ((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec])*((np->P[ec_2][0][i-s_2+s_ec]) + op->P[ec_2][0][i-s_2+s_ec]));
                          const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                             - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                                - (op->P[ec][0][i]))) + (0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))) + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])))) + 0.25*(((((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))*((np->P[ec_2][1][i]) - (op->P[ec_2][1][i]))) + ((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))*((np->P[ec_2][0][i]) - (op->P[ec_2][0][i]))) + ((((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))*((np->P[ec_2][1][i-s_2]) - (op->P[ec_2][1][i-s_2]))) + ((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))*((np->P[ec_2][0][i-s_2]) - (op->P[ec_2][0][i-s_2]))) + ((((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))*((np->P[ec_2][1][i+s_ec]) - (op->P[ec_2][1][i+s_ec]))) + ((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))*((np->P[ec_2][0][i+s_ec]) - (op->P[ec_2][0][i+s_ec]))) + (((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))*((np->P[ec_2][1][i-s_2+s_ec]) - (op->P[ec_2][1][i-s_2+s_ec]))) + ((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec]))*((np->P[ec_2][0][i-s_2+s_ec]) - (op->P[ec_2][0][i-s_2+s_ec])));
                          const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                          np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                          op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                             + np->s[ec][i]*f[ec][0][i]);
                          op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                             + np->s[ec][i]*f[ec][1][i]);
                        }
                      }
                    }
                  }
                } else { // not np->saturation_factor
                  for (int i=0;i<ntot;i++) {
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                    np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                    op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
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
      FOR_E_AND_D(ec,dc) if (f[ec][0]) {
        const int yee_idx = v.yee_index(ec);
        const int d_ec = component_direction(ec);
        const int s_ec = stride_any_direction[d_ec];
        const direction d_1 = (direction)((d_ec+1)%2);
        const component ec_1 = direction_component(ec,d_1);
        const int s_1 = stride_any_direction[d_1];
        if (is_real) {
          if (pol) {
            for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
              const double om_psqr = fabs(op->pb->saturated_sigma);
              const double oo_ep_om_psqr = 1.0/(8*pi*om_psqr);
              const double oo_e_sat = 1.0/op->pb->energy_saturation;
              const double om_sqr = op->pb->omeganot*op->pb->omeganot;
              const double g = op->pb->gamma;
              const double funinv = 1.0/(1+0.5*g);
              if (np->saturation_factor) {
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                }
                if (num_any_direction[X]==1) {
                  if (num_any_direction[Y]==1) {
                    const int i = yee_idx;
                    const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                       + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                          + np->energy[ec_1][i-s_1+s_ec]);
                    const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                       + op->P[ec][0][i])) + 0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                          + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                             + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                   + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                      + op->P[ec_1][0][i-s_1+s_ec]));
                    const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                       - (op->P[ec][0][i]))) + 0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                          - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                             - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])));
                    const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                    np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                    op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  } else { // not num_any_direction[Y]==1
                    for (int iY=0; iY<num_any_direction[Y]; iY++) {
                      const int i = yee_idx + iY;
                      const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                         + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                            + np->energy[ec_1][i-s_1+s_ec]);
                      const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                         + op->P[ec][0][i])) + 0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                            + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                               + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                  op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                     + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                        + op->P[ec_1][0][i-s_1+s_ec]));
                      const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                         - (op->P[ec][0][i]))) + 0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                            - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                               - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                  - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])));
                      const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                      np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                      op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    }
                  }
                } else { // not num_any_direction[X]==1
                  if (num_any_direction[Y]==1) {
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      const int i = yee_idx + iX*stride_any_direction[X];
                      const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                         + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                            + np->energy[ec_1][i-s_1+s_ec]);
                      const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                         + op->P[ec][0][i])) + 0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                            + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                               + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                  op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                     + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                        + op->P[ec_1][0][i-s_1+s_ec]));
                      const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                         - (op->P[ec][0][i]))) + 0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                            - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                               - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                  - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])));
                      const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                      np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                      op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    }
                  } else { // not num_any_direction[Y]==1
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iX*stride_any_direction[X] + iY;
                        const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                           + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                              + np->energy[ec_1][i-s_1+s_ec]);
                        const double P_sqr = (((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                           + op->P[ec][0][i])) + 0.25*((((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                              + op->P[ec_1][0][i])) + (((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                 + op->P[ec_1][0][i-s_1])) + (((np->P[ec_1][0][i+s_ec]) +
                                    op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec]))
                                       + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec])
                                          + op->P[ec_1][0][i-s_1+s_ec]));
                        const double Pdot_sqr = (((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                           - (op->P[ec][0][i]))) + 0.25*((((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                              - (op->P[ec_1][0][i]))) + (((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1])
                                 - (op->P[ec_1][0][i-s_1]))) + (((np->P[ec_1][0][i+s_ec])
                                    - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])));
                        const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                        np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                        op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      }
                    }
                  }
                }
              } else { // not np->saturation_factor
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                }
              }
            }
          }
        } else { // not is_real
          if (pol) {
            for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
              const double om_psqr = fabs(op->pb->saturated_sigma);
              const double oo_ep_om_psqr = 1.0/(8*pi*om_psqr);
              const double oo_e_sat = 1.0/op->pb->energy_saturation;
              const double om_sqr = op->pb->omeganot*op->pb->omeganot;
              const double g = op->pb->gamma;
              const double funinv = 1.0/(1+0.5*g);
              if (np->saturation_factor) {
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                }
                if (num_any_direction[X]==1) {
                  if (num_any_direction[Y]==1) {
                    const int i = yee_idx;
                    const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                       + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                          + np->energy[ec_1][i-s_1+s_ec]);
                    const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                       + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                          + op->P[ec][0][i])) + 0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                             + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                   + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                      + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]));
                    const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                       - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                          - (op->P[ec][0][i]))) + 0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i])
                             - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                                - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1])
                                   - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) -
                                      (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1])))
                                         + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec])
                                            - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec])
                                               - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])));
                    const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                    np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                    op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  } else { // not num_any_direction[Y]==1
                    for (int iY=0; iY<num_any_direction[Y]; iY++) {
                      const int i = yee_idx + iY;
                      const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                         + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                            + np->energy[ec_1][i-s_1+s_ec]);
                      const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                         + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                            + op->P[ec][0][i])) + 0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                               + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                  + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                     + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                        + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]));
                      const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                         - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                            - (op->P[ec][0][i]))) + 0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i])
                               - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                                  - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1])
                                     - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) -
                                        (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1])))
                                           + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec])
                                              - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec])
                                                 - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])));
                      const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                      np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                      op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                      op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                         + np->s[ec][i]*f[ec][1][i]);
                    }
                  }
                } else { // not num_any_direction[X]==1
                  if (num_any_direction[Y]==1) {
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      const int i = yee_idx + iX*stride_any_direction[X];
                      const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                         + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                            + np->energy[ec_1][i-s_1+s_ec]);
                      const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                         + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                            + op->P[ec][0][i])) + 0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                               + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                  + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                     + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                        + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]));
                      const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                         - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                            - (op->P[ec][0][i]))) + 0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i])
                               - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                                  - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1])
                                     - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) -
                                        (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1])))
                                           + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec])
                                              - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec])
                                                 - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])));
                      const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                      np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                      op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                      op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                         + np->s[ec][i]*f[ec][1][i]);
                    }
                  } else { // not num_any_direction[Y]==1
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iX*stride_any_direction[X] + iY;
                        const double energy_here = (np->energy[ec][i]) + 0.25*((np->energy[ec_1][i])
                           + (np->energy[ec_1][i-s_1]) + (np->energy[ec_1][i+s_ec])
                              + np->energy[ec_1][i-s_1+s_ec]);
                        const double P_sqr = ((((np->P[ec][1][i]) + op->P[ec][1][i])*((np->P[ec][1][i])
                           + op->P[ec][1][i])) + ((np->P[ec][0][i]) + op->P[ec][0][i])*((np->P[ec][0][i])
                              + op->P[ec][0][i])) + 0.25*(((((np->P[ec_1][1][i]) + op->P[ec_1][1][i])*((np->P[ec_1][1][i])
                                 + op->P[ec_1][1][i])) + ((np->P[ec_1][0][i]) + op->P[ec_1][0][i])*((np->P[ec_1][0][i])
                                    + op->P[ec_1][0][i])) + ((((np->P[ec_1][1][i-s_1]) + op->P[ec_1][1][i-s_1])*((np->P[ec_1][1][i-s_1])
                                       + op->P[ec_1][1][i-s_1])) + ((np->P[ec_1][0][i-s_1]) + op->P[ec_1][0][i-s_1])*((np->P[ec_1][0][i-s_1])
                                          + op->P[ec_1][0][i-s_1])) + ((((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])*((np->P[ec_1][1][i+s_ec]) + op->P[ec_1][1][i+s_ec])) + ((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])*((np->P[ec_1][0][i+s_ec]) + op->P[ec_1][0][i+s_ec])) + (((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])*((np->P[ec_1][1][i-s_1+s_ec]) + op->P[ec_1][1][i-s_1+s_ec])) + ((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec])*((np->P[ec_1][0][i-s_1+s_ec]) + op->P[ec_1][0][i-s_1+s_ec]));
                        const double Pdot_sqr = ((((np->P[ec][1][i]) - (op->P[ec][1][i]))*((np->P[ec][1][i])
                           - (op->P[ec][1][i]))) + ((np->P[ec][0][i]) - (op->P[ec][0][i]))*((np->P[ec][0][i])
                              - (op->P[ec][0][i]))) + 0.25*(((((np->P[ec_1][1][i]) - (op->P[ec_1][1][i]))*((np->P[ec_1][1][i])
                                 - (op->P[ec_1][1][i]))) + ((np->P[ec_1][0][i]) - (op->P[ec_1][0][i]))*((np->P[ec_1][0][i])
                                    - (op->P[ec_1][0][i]))) + ((((np->P[ec_1][1][i-s_1]) - (op->P[ec_1][1][i-s_1]))*((np->P[ec_1][1][i-s_1])
                                       - (op->P[ec_1][1][i-s_1]))) + ((np->P[ec_1][0][i-s_1]) -
                                          (op->P[ec_1][0][i-s_1]))*((np->P[ec_1][0][i-s_1]) - (op->P[ec_1][0][i-s_1])))
                                             + ((((np->P[ec_1][1][i+s_ec]) - (op->P[ec_1][1][i+s_ec]))*((np->P[ec_1][1][i+s_ec])
                                                - (op->P[ec_1][1][i+s_ec]))) + ((np->P[ec_1][0][i+s_ec])
                                                   - (op->P[ec_1][0][i+s_ec]))*((np->P[ec_1][0][i+s_ec]) - (op->P[ec_1][0][i+s_ec]))) + (((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))*((np->P[ec_1][1][i-s_1+s_ec]) - (op->P[ec_1][1][i-s_1+s_ec]))) + ((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec]))*((np->P[ec_1][0][i-s_1+s_ec]) - (op->P[ec_1][0][i-s_1+s_ec])));
                        const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                        np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                        op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                        op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                           + np->s[ec][i]*f[ec][1][i]);
                      }
                    }
                  }
                }
              } else { // not np->saturation_factor
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                  op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                  op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
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
      FOR_E_AND_D(ec,dc) if (f[ec][0]) {
        const int yee_idx = v.yee_index(ec);
        const int d_ec = component_direction(ec);
        const int s_ec = stride_any_direction[d_ec];
        if (is_real) {
          if (pol) {
            for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
              const double om_psqr = fabs(op->pb->saturated_sigma);
              const double oo_ep_om_psqr = 1.0/(8*pi*om_psqr);
              const double oo_e_sat = 1.0/op->pb->energy_saturation;
              const double om_sqr = op->pb->omeganot*op->pb->omeganot;
              const double g = op->pb->gamma;
              const double funinv = 1.0/(1+0.5*g);
              if (np->saturation_factor) {
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                }
                if (num_any_direction[Z]==1) {
                  const int i = yee_idx;
                  const double energy_here = np->energy[Ex][i];
                  const double P_sqr = ((np->P[Ex][0][i]) + op->P[Ex][0][i])*((np->P[Ex][0][i])
                     + op->P[Ex][0][i]);
                  const double Pdot_sqr = ((np->P[Ex][0][i]) - (op->P[Ex][0][i]))*((np->P[Ex][0][i])
                     - (op->P[Ex][0][i]));
                  const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                  np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                  op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                } else { // not num_any_direction[Z]==1
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    const double energy_here = np->energy[Ex][i];
                    const double P_sqr = ((np->P[Ex][0][i]) + op->P[Ex][0][i])*((np->P[Ex][0][i])
                       + op->P[Ex][0][i]);
                    const double Pdot_sqr = ((np->P[Ex][0][i]) - (op->P[Ex][0][i]))*((np->P[Ex][0][i])
                       - (op->P[Ex][0][i]));
                    const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                    np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                    op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  }
                }
              } else { // not np->saturation_factor
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                }
              }
            }
          }
        } else { // not is_real
          if (pol) {
            for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
              const double om_psqr = fabs(op->pb->saturated_sigma);
              const double oo_ep_om_psqr = 1.0/(8*pi*om_psqr);
              const double oo_e_sat = 1.0/op->pb->energy_saturation;
              const double om_sqr = op->pb->omeganot*op->pb->omeganot;
              const double g = op->pb->gamma;
              const double funinv = 1.0/(1+0.5*g);
              if (np->saturation_factor) {
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                }
                if (num_any_direction[Z]==1) {
                  const int i = yee_idx;
                  const double energy_here = np->energy[Ex][i];
                  const double P_sqr = (((np->P[Ex][1][i]) + op->P[Ex][1][i])*((np->P[Ex][1][i])
                     + op->P[Ex][1][i])) + ((np->P[Ex][0][i]) + op->P[Ex][0][i])*((np->P[Ex][0][i])
                        + op->P[Ex][0][i]);
                  const double Pdot_sqr = (((np->P[Ex][1][i]) - (op->P[Ex][1][i]))*((np->P[Ex][1][i])
                     - (op->P[Ex][1][i]))) + ((np->P[Ex][0][i]) - (op->P[Ex][0][i]))*((np->P[Ex][0][i])
                        - (op->P[Ex][0][i]));
                  const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                  np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                  op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                  op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                     + np->s[ec][i]*f[ec][1][i]);
                } else { // not num_any_direction[Z]==1
                  for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
                    const int i = yee_idx + iZ;
                    const double energy_here = np->energy[Ex][i];
                    const double P_sqr = (((np->P[Ex][1][i]) + op->P[Ex][1][i])*((np->P[Ex][1][i])
                       + op->P[Ex][1][i])) + ((np->P[Ex][0][i]) + op->P[Ex][0][i])*((np->P[Ex][0][i])
                          + op->P[Ex][0][i]);
                    const double Pdot_sqr = (((np->P[Ex][1][i]) - (op->P[Ex][1][i]))*((np->P[Ex][1][i])
                       - (op->P[Ex][1][i]))) + ((np->P[Ex][0][i]) - (op->P[Ex][0][i]))*((np->P[Ex][0][i])
                          - (op->P[Ex][0][i]));
                    const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                    np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                    op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  }
                }
              } else { // not np->saturation_factor
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                  op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                  op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
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
      FOR_E_AND_D(ec,dc) if (f[ec][0]) {
        const int yee_idx = v.yee_index(ec);
        const int d_ec = component_direction(ec);
        const int s_ec = stride_any_direction[d_ec];
        if (is_real) {
          if (pol) {
            for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
              const double om_psqr = fabs(op->pb->saturated_sigma);
              const double oo_ep_om_psqr = 1.0/(8*pi*om_psqr);
              const double oo_e_sat = 1.0/op->pb->energy_saturation;
              const double om_sqr = op->pb->omeganot*op->pb->omeganot;
              const double g = op->pb->gamma;
              const double funinv = 1.0/(1+0.5*g);
              if (np->saturation_factor) {
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                }
                if (num_any_direction[X]==1) {
                  if (num_any_direction[Y]==1) {
                    const int i = yee_idx;
                    const double energy_here = np->energy[Ez][i];
                    const double P_sqr = ((np->P[Ez][0][i]) + op->P[Ez][0][i])*((np->P[Ez][0][i])
                       + op->P[Ez][0][i]);
                    const double Pdot_sqr = ((np->P[Ez][0][i]) - (op->P[Ez][0][i]))*((np->P[Ez][0][i])
                       - (op->P[Ez][0][i]));
                    const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                    np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                    op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                  } else { // not num_any_direction[Y]==1
                    for (int iY=0; iY<num_any_direction[Y]; iY++) {
                      const int i = yee_idx + iY;
                      const double energy_here = np->energy[Ez][i];
                      const double P_sqr = ((np->P[Ez][0][i]) + op->P[Ez][0][i])*((np->P[Ez][0][i])
                         + op->P[Ez][0][i]);
                      const double Pdot_sqr = ((np->P[Ez][0][i]) - (op->P[Ez][0][i]))*((np->P[Ez][0][i])
                         - (op->P[Ez][0][i]));
                      const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                      np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                      op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    }
                  }
                } else { // not num_any_direction[X]==1
                  if (num_any_direction[Y]==1) {
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      const int i = yee_idx + iX*stride_any_direction[X];
                      const double energy_here = np->energy[Ez][i];
                      const double P_sqr = ((np->P[Ez][0][i]) + op->P[Ez][0][i])*((np->P[Ez][0][i])
                         + op->P[Ez][0][i]);
                      const double Pdot_sqr = ((np->P[Ez][0][i]) - (op->P[Ez][0][i]))*((np->P[Ez][0][i])
                         - (op->P[Ez][0][i]));
                      const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                      np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                      op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                    }
                  } else { // not num_any_direction[Y]==1
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iX*stride_any_direction[X] + iY;
                        const double energy_here = np->energy[Ez][i];
                        const double P_sqr = ((np->P[Ez][0][i]) + op->P[Ez][0][i])*((np->P[Ez][0][i])
                           + op->P[Ez][0][i]);
                        const double Pdot_sqr = ((np->P[Ez][0][i]) - (op->P[Ez][0][i]))*((np->P[Ez][0][i])
                           - (op->P[Ez][0][i]));
                        const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                        np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                        op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                      }
                    }
                  }
                }
              } else { // not np->saturation_factor
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                }
              }
            }
          }
        } else { // not is_real
          if (pol) {
            for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
              const double om_psqr = fabs(op->pb->saturated_sigma);
              const double oo_ep_om_psqr = 1.0/(8*pi*om_psqr);
              const double oo_e_sat = 1.0/op->pb->energy_saturation;
              const double om_sqr = op->pb->omeganot*op->pb->omeganot;
              const double g = op->pb->gamma;
              const double funinv = 1.0/(1+0.5*g);
              if (np->saturation_factor) {
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                }
                if (num_any_direction[X]==1) {
                  if (num_any_direction[Y]==1) {
                    const int i = yee_idx;
                    const double energy_here = np->energy[Ez][i];
                    const double P_sqr = (((np->P[Ez][1][i]) + op->P[Ez][1][i])*((np->P[Ez][1][i])
                       + op->P[Ez][1][i])) + ((np->P[Ez][0][i]) + op->P[Ez][0][i])*((np->P[Ez][0][i])
                          + op->P[Ez][0][i]);
                    const double Pdot_sqr = (((np->P[Ez][1][i]) - (op->P[Ez][1][i]))*((np->P[Ez][1][i])
                       - (op->P[Ez][1][i]))) + ((np->P[Ez][0][i]) - (op->P[Ez][0][i]))*((np->P[Ez][0][i])
                          - (op->P[Ez][0][i]));
                    const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                    np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                    op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                       + np->s[ec][i]*f[ec][0][i]);
                    op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                       + np->s[ec][i]*f[ec][1][i]);
                  } else { // not num_any_direction[Y]==1
                    for (int iY=0; iY<num_any_direction[Y]; iY++) {
                      const int i = yee_idx + iY;
                      const double energy_here = np->energy[Ez][i];
                      const double P_sqr = (((np->P[Ez][1][i]) + op->P[Ez][1][i])*((np->P[Ez][1][i])
                         + op->P[Ez][1][i])) + ((np->P[Ez][0][i]) + op->P[Ez][0][i])*((np->P[Ez][0][i])
                            + op->P[Ez][0][i]);
                      const double Pdot_sqr = (((np->P[Ez][1][i]) - (op->P[Ez][1][i]))*((np->P[Ez][1][i])
                         - (op->P[Ez][1][i]))) + ((np->P[Ez][0][i]) - (op->P[Ez][0][i]))*((np->P[Ez][0][i])
                            - (op->P[Ez][0][i]));
                      const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                      np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                      op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                      op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                         + np->s[ec][i]*f[ec][1][i]);
                    }
                  }
                } else { // not num_any_direction[X]==1
                  if (num_any_direction[Y]==1) {
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      const int i = yee_idx + iX*stride_any_direction[X];
                      const double energy_here = np->energy[Ez][i];
                      const double P_sqr = (((np->P[Ez][1][i]) + op->P[Ez][1][i])*((np->P[Ez][1][i])
                         + op->P[Ez][1][i])) + ((np->P[Ez][0][i]) + op->P[Ez][0][i])*((np->P[Ez][0][i])
                            + op->P[Ez][0][i]);
                      const double Pdot_sqr = (((np->P[Ez][1][i]) - (op->P[Ez][1][i]))*((np->P[Ez][1][i])
                         - (op->P[Ez][1][i]))) + ((np->P[Ez][0][i]) - (op->P[Ez][0][i]))*((np->P[Ez][0][i])
                            - (op->P[Ez][0][i]));
                      const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                      np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                      op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                         + np->s[ec][i]*f[ec][0][i]);
                      op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                         + np->s[ec][i]*f[ec][1][i]);
                    }
                  } else { // not num_any_direction[Y]==1
                    for (int iX=0; iX<num_any_direction[X]; iX++) {
                      for (int iY=0; iY<num_any_direction[Y]; iY++) {
                        const int i = yee_idx + iX*stride_any_direction[X] + iY;
                        const double energy_here = np->energy[Ez][i];
                        const double P_sqr = (((np->P[Ez][1][i]) + op->P[Ez][1][i])*((np->P[Ez][1][i])
                           + op->P[Ez][1][i])) + ((np->P[Ez][0][i]) + op->P[Ez][0][i])*((np->P[Ez][0][i])
                              + op->P[Ez][0][i]);
                        const double Pdot_sqr = (((np->P[Ez][1][i]) - (op->P[Ez][1][i]))*((np->P[Ez][1][i])
                           - (op->P[Ez][1][i]))) + ((np->P[Ez][0][i]) - (op->P[Ez][0][i]))*((np->P[Ez][0][i])
                              - (op->P[Ez][0][i]));
                        const double energy_P_here = (om_sqr*P_sqr + (1.0/c/c)*Pdot_sqr)*oo_ep_om_psqr;
                        np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat);
                        op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                           + np->s[ec][i]*f[ec][0][i]);
                        op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                           + np->s[ec][i]*f[ec][1][i]);
                      }
                    }
                  }
                }
              } else { // not np->saturation_factor
                for (int i=0;i<ntot;i++) {
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
                  np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
                  op->P[ec][0][i] = funinv*((2-om_sqr)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                     + np->s[ec][i]*f[ec][0][i]);
                  op->P[ec][1][i] = funinv*((2-om_sqr)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
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
