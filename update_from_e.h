FOR_E_AND_D(ec,dc) if (f[ec][0]) {
  if (is_real) {
    if (pol) {
      for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
        const double fac = np->saturation_factor;
        const double g = op->pb->gamma;;
        const double om = op->pb->omeganot;;
        const double funinv = 1.0/(1+0.5*g);;
        if (fac) {
          if (fac > 0) {
            for (int i=0;i<ntot;i++) {
              np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
              np->s[ec][i] = max(-np->energy[ec][i]*fac, 0.0);
              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                 + np->s[ec][i]*f[ec][0][i]);
            }
          } else {
            for (int i=0;i<ntot;i++) {
              np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
              np->s[ec][i] = np->energy[ec][i]*fac;
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
        const double g = op->pb->gamma;;
        const double om = op->pb->omeganot;;
        const double funinv = 1.0/(1+0.5*g);;
        if (fac) {
          if (fac > 0) {
            for (int i=0;i<ntot;i++) {
              np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
              np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
              np->s[ec][i] = max(-np->energy[ec][i]*fac, 0.0);
              op->P[ec][0][i] = funinv*((2-om*om)*np->P[ec][0][i] + (0.5*g-1)*op->P[ec][0][i]
                 + np->s[ec][i]*f[ec][0][i]);
              op->P[ec][1][i] = funinv*((2-om*om)*np->P[ec][1][i] + (0.5*g-1)*op->P[ec][1][i]
                 + np->s[ec][i]*f[ec][1][i]);
            }
          } else {
            for (int i=0;i<ntot;i++) {
              np->energy[ec][i] += 0.5*(np->P[ec][0][i] - op->P[ec][0][i])*f[ec][0][i];
              np->energy[ec][i] += 0.5*(np->P[ec][1][i] - op->P[ec][1][i])*f[ec][1][i];
              np->s[ec][i] = np->energy[ec][i]*fac;
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
