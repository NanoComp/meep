if (have_m) {
  if (have_p) {
    if (have_m_pml) {
      if (have_p_pml) {
        if (n3==1) {
          if (s2==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
                const double em = the_f_m_pml[ind];
                const double ep = the_f_p_pml[ind];
                const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                const double deriv_p = f_p[ind+stride_p]-f_p[ind];
                const double m_change = decay_m[ind]*(c*m_deriv_m - C_m[ind]*em);
                const double p_change = decay_p[ind]*(c*deriv_p - C_p[ind]*ep);
                the_f_p_pml[ind] += p_change;
                the_f_m_pml[ind] += m_change;
                the_f[ind] += m_change + p_change;
              }
          } else {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
                const double em = the_f_m_pml[ind];
                const double ep = the_f_p_pml[ind];
                const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                const double deriv_p = f_p[ind+stride_p]-f_p[ind];
                const double m_change = decay_m[ind]*(c*m_deriv_m - C_m[ind]*em);
                const double p_change = decay_p[ind]*(c*deriv_p - C_p[ind]*ep);
                the_f_p_pml[ind] += p_change;
                the_f_m_pml[ind] += m_change;
                the_f[ind] += m_change + p_change;
              }
          }
        } else {
          if (s3==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                  const double em = the_f_m_pml[ind];
                  const double ep = the_f_p_pml[ind];
                  const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                  const double deriv_p = f_p[ind+stride_p]-f_p[ind];
                  const double m_change = decay_m[ind]*(c*m_deriv_m - C_m[ind]*em);
                  const double p_change = decay_p[ind]*(c*deriv_p - C_p[ind]*ep);
                  the_f_p_pml[ind] += p_change;
                  the_f_m_pml[ind] += m_change;
                  the_f[ind] += m_change + p_change;
                }
          } else {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                  const double em = the_f_m_pml[ind];
                  const double ep = the_f_p_pml[ind];
                  const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                  const double deriv_p = f_p[ind+stride_p]-f_p[ind];
                  const double m_change = decay_m[ind]*(c*m_deriv_m - C_m[ind]*em);
                  const double p_change = decay_p[ind]*(c*deriv_p - C_p[ind]*ep);
                  the_f_p_pml[ind] += p_change;
                  the_f_m_pml[ind] += m_change;
                  the_f[ind] += m_change + p_change;
                }
          }
        }
      } else {
        if (n3==1) {
          if (s2==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
                const double em = the_f_m_pml[ind];
                const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                const double deriv_p = f_p[ind+stride_p]-f_p[ind];
                const double m_change = decay_m[ind]*(c*m_deriv_m - C_m[ind]*em);
                const double p_change = inveps[ind]*c*deriv_p;
                the_f_m_pml[ind] += m_change;
                the_f[ind] += m_change + p_change;
              }
          } else {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
                const double em = the_f_m_pml[ind];
                const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                const double deriv_p = f_p[ind+stride_p]-f_p[ind];
                const double m_change = decay_m[ind]*(c*m_deriv_m - C_m[ind]*em);
                const double p_change = inveps[ind]*c*deriv_p;
                the_f_m_pml[ind] += m_change;
                the_f[ind] += m_change + p_change;
              }
          }
        } else {
          if (s3==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                  const double em = the_f_m_pml[ind];
                  const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                  const double deriv_p = f_p[ind+stride_p]-f_p[ind];
                  const double m_change = decay_m[ind]*(c*m_deriv_m - C_m[ind]*em);
                  const double p_change = inveps[ind]*c*deriv_p;
                  the_f_m_pml[ind] += m_change;
                  the_f[ind] += m_change + p_change;
                }
          } else {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                  const double em = the_f_m_pml[ind];
                  const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                  const double deriv_p = f_p[ind+stride_p]-f_p[ind];
                  const double m_change = decay_m[ind]*(c*m_deriv_m - C_m[ind]*em);
                  const double p_change = inveps[ind]*c*deriv_p;
                  the_f_m_pml[ind] += m_change;
                  the_f[ind] += m_change + p_change;
                }
          }
        }
      }
    } else {
      if (have_p_pml) {
        if (n3==1) {
          if (s2==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
                const double ep = the_f_p_pml[ind];
                const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                const double deriv_p = f_p[ind+stride_p]-f_p[ind];
                const double m_change = inveps[ind]*c*m_deriv_m;
                const double p_change = decay_p[ind]*(c*deriv_p - C_p[ind]*ep);
                the_f_p_pml[ind] += p_change;
                the_f[ind] += m_change + p_change;
              }
          } else {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
                const double ep = the_f_p_pml[ind];
                const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                const double deriv_p = f_p[ind+stride_p]-f_p[ind];
                const double m_change = inveps[ind]*c*m_deriv_m;
                const double p_change = decay_p[ind]*(c*deriv_p - C_p[ind]*ep);
                the_f_p_pml[ind] += p_change;
                the_f[ind] += m_change + p_change;
              }
          }
        } else {
          if (s3==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                  const double ep = the_f_p_pml[ind];
                  const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                  const double deriv_p = f_p[ind+stride_p]-f_p[ind];
                  const double m_change = inveps[ind]*c*m_deriv_m;
                  const double p_change = decay_p[ind]*(c*deriv_p - C_p[ind]*ep);
                  the_f_p_pml[ind] += p_change;
                  the_f[ind] += m_change + p_change;
                }
          } else {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                  const double ep = the_f_p_pml[ind];
                  const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                  const double deriv_p = f_p[ind+stride_p]-f_p[ind];
                  const double m_change = inveps[ind]*c*m_deriv_m;
                  const double p_change = decay_p[ind]*(c*deriv_p - C_p[ind]*ep);
                  the_f_p_pml[ind] += p_change;
                  the_f[ind] += m_change + p_change;
                }
          }
        }
      } else {
        if (n3==1) {
          if (s2==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
                const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                const double deriv_p = f_p[ind+stride_p]-f_p[ind];
                const double m_change = inveps[ind]*c*m_deriv_m;
                const double p_change = inveps[ind]*c*deriv_p;
                the_f[ind] += inveps[ind]*c*(m_deriv_m + deriv_p);
              }
          } else {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
                const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                const double deriv_p = f_p[ind+stride_p]-f_p[ind];
                const double m_change = inveps[ind]*c*m_deriv_m;
                const double p_change = inveps[ind]*c*deriv_p;
                the_f[ind] += inveps[ind]*c*(m_deriv_m + deriv_p);
              }
          }
        } else {
          if (s3==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                  const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                  const double deriv_p = f_p[ind+stride_p]-f_p[ind];
                  const double m_change = inveps[ind]*c*m_deriv_m;
                  const double p_change = inveps[ind]*c*deriv_p;
                  the_f[ind] += inveps[ind]*c*(m_deriv_m + deriv_p);
                }
          } else {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                  const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                  const double deriv_p = f_p[ind+stride_p]-f_p[ind];
                  const double m_change = inveps[ind]*c*m_deriv_m;
                  const double p_change = inveps[ind]*c*deriv_p;
                  the_f[ind] += inveps[ind]*c*(m_deriv_m + deriv_p);
                }
          }
        }
      }
    }
  } else {
    if (have_m_pml) {
      if (n3==1) {
        if (s2==1) {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
              const double em = the_f[ind];
              const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
              const double m_change = decay_m[ind]*(c*m_deriv_m - C_m[ind]*em);
              the_f[ind] += m_change;
            }
        } else {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
              const double em = the_f[ind];
              const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
              const double m_change = decay_m[ind]*(c*m_deriv_m - C_m[ind]*em);
              the_f[ind] += m_change;
            }
        }
      } else {
        if (s3==1) {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0; i2<n2; i2++)
              for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                const double em = the_f[ind];
                const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                const double m_change = decay_m[ind]*(c*m_deriv_m - C_m[ind]*em);
                the_f[ind] += m_change;
              }
        } else {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0; i2<n2; i2++)
              for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                const double em = the_f[ind];
                const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                const double m_change = decay_m[ind]*(c*m_deriv_m - C_m[ind]*em);
                the_f[ind] += m_change;
              }
        }
      }
    } else {
      if (n3==1) {
        if (s2==1) {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
              const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
              const double m_change = inveps[ind]*c*m_deriv_m;
              the_f[ind] += inveps[ind]*c*m_deriv_m;
            }
        } else {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
              const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
              const double m_change = inveps[ind]*c*m_deriv_m;
              the_f[ind] += inveps[ind]*c*m_deriv_m;
            }
        }
      } else {
        if (s3==1) {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0; i2<n2; i2++)
              for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                const double m_change = inveps[ind]*c*m_deriv_m;
                the_f[ind] += inveps[ind]*c*m_deriv_m;
              }
        } else {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0; i2<n2; i2++)
              for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                const double m_deriv_m = f_m[ind]-f_m[ind+stride_m];
                const double m_change = inveps[ind]*c*m_deriv_m;
                the_f[ind] += inveps[ind]*c*m_deriv_m;
              }
        }
      }
    }
  }
} else {
  if (have_p_pml) {
    if (n3==1) {
      if (s2==1) {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
            const double ep = the_f[ind];
            const double deriv_p = f_p[ind+stride_p]-f_p[ind];
            const double p_change = decay_p[ind]*(c*deriv_p - C_p[ind]*ep);
            the_f[ind] += p_change;
          }
      } else {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
            const double ep = the_f[ind];
            const double deriv_p = f_p[ind+stride_p]-f_p[ind];
            const double p_change = decay_p[ind]*(c*deriv_p - C_p[ind]*ep);
            the_f[ind] += p_change;
          }
      }
    } else {
      if (s3==1) {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0; i2<n2; i2++)
            for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
              const double ep = the_f[ind];
              const double deriv_p = f_p[ind+stride_p]-f_p[ind];
              const double p_change = decay_p[ind]*(c*deriv_p - C_p[ind]*ep);
              the_f[ind] += p_change;
            }
      } else {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0; i2<n2; i2++)
            for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
              const double ep = the_f[ind];
              const double deriv_p = f_p[ind+stride_p]-f_p[ind];
              const double p_change = decay_p[ind]*(c*deriv_p - C_p[ind]*ep);
              the_f[ind] += p_change;
            }
      }
    }
  } else {
    if (n3==1) {
      if (s2==1) {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
            const double deriv_p = f_p[ind+stride_p]-f_p[ind];
            const double p_change = inveps[ind]*c*deriv_p;
            the_f[ind] += inveps[ind]*c*deriv_p;
          }
      } else {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
            const double deriv_p = f_p[ind+stride_p]-f_p[ind];
            const double p_change = inveps[ind]*c*deriv_p;
            the_f[ind] += inveps[ind]*c*deriv_p;
          }
      }
    } else {
      if (s3==1) {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0; i2<n2; i2++)
            for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
              const double deriv_p = f_p[ind+stride_p]-f_p[ind];
              const double p_change = inveps[ind]*c*deriv_p;
              the_f[ind] += inveps[ind]*c*deriv_p;
            }
      } else {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0; i2<n2; i2++)
            for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
              const double deriv_p = f_p[ind+stride_p]-f_p[ind];
              const double p_change = inveps[ind]*c*deriv_p;
              the_f[ind] += inveps[ind]*c*deriv_p;
            }
      }
    }
  }
}
