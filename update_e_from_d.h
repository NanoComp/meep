if (pol) {
  FOR_E_AND_D(ec,dc) if (f[ec][0]) {
    const double *the_inveps = ma->inveps[ec][component_direction(ec)];
    DOCMP {
      for (int i=0;i<ntot;i++) {
        double the_polarization = 0.0;
        FOR_POLARIZATIONS(pol, p) the_polarization += p->P[ec][cmp][i];
        f[ec][cmp][i] = the_inveps[i]*(f[dc][cmp][i] - the_polarization);
      }
    }
  }
} else {
  FOR_E_AND_D(ec,dc) if (f[ec][0]) {
    const double *the_inveps = ma->inveps[ec][component_direction(ec)];
    DOCMP {
      for (int i=0;i<ntot;i++) {
        f[ec][cmp][i] = the_inveps[i]*f[dc][cmp][i];
      }
    }
  }
}
