FOR_E_AND_D(ec,dc) if (f[ec][0]) {
  if (pol) {
    for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
      const double fac = np->saturation_factor;
      const double g = op->pb->gamma;;
      const double om = op->pb->omeganot;;
      const double funinv = 1.0/(1+0.5*g);;
      if (fac) {
        if (fac > 0) {
          for (int i=0;i<ntot;i++) {
            DOCMP {
              np->energy[ec][i] += 0.5*(np->P[ec][cmp][i] - op->P[ec][cmp][i])*f[ec][cmp][i];
            }
            np->s[ec][i] = max(-np->energy[ec][i]*fac, 0.0);
            DOCMP {
              op->P[ec][cmp][i] = funinv*((2-om*om)*np->P[ec][cmp][i] + (0.5*g-1)*op->P[ec][cmp][i]) + np->s[ec][i]*f[ec][cmp][i];
            }
          }
        } else {
          for (int i=0;i<ntot;i++) {
            DOCMP {
              np->energy[ec][i] += 0.5*(np->P[ec][cmp][i] - op->P[ec][cmp][i])*f[ec][cmp][i];
            }
            np->s[ec][i] = np->energy[ec][i]*fac;
            DOCMP {
              op->P[ec][cmp][i] = funinv*((2-om*om)*np->P[ec][cmp][i] + (0.5*g-1)*op->P[ec][cmp][i]) + np->s[ec][i]*f[ec][cmp][i];
            }
          }
        }
      } else {
        for (int i=0;i<ntot;i++) {
          DOCMP {
            np->energy[ec][i] += 0.5*(np->P[ec][cmp][i] - op->P[ec][cmp][i])*f[ec][cmp][i];
          }
          DOCMP {
            op->P[ec][cmp][i] = funinv*((2-om*om)*np->P[ec][cmp][i] + (0.5*g-1)*op->P[ec][cmp][i]) + np->s[ec][i]*f[ec][cmp][i];
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
