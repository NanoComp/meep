if (have_m) {
  if (have_p) {
    if (have_m_pml) {
      if (have_p_pml) {
        if (n3==1) {
          if (s2==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
                const double hm = the_f_m_pml[ind];
                const double hp = the_f_p_pml[ind];
                const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
                the_f_p_pml[ind] += decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
                the_f_m_pml[ind] += decay_m[ind]*(c*deriv_m - C_m[ind]*hm);
                the_f[ind] += (decay_m[ind]*(c*deriv_m - C_m[ind]*hm)) +
                   decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
              }
          } else { // not s2==1
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
                const double hm = the_f_m_pml[ind];
                const double hp = the_f_p_pml[ind];
                const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
                the_f_p_pml[ind] += decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
                the_f_m_pml[ind] += decay_m[ind]*(c*deriv_m - C_m[ind]*hm);
                the_f[ind] += (decay_m[ind]*(c*deriv_m - C_m[ind]*hm)) +
                   decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
              }
          }
        } else { // not n3==1
          if (s3==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                  const double hm = the_f_m_pml[ind];
                  const double hp = the_f_p_pml[ind];
                  const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                  const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
                  the_f_p_pml[ind] += decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
                  the_f_m_pml[ind] += decay_m[ind]*(c*deriv_m - C_m[ind]*hm);
                  the_f[ind] += (decay_m[ind]*(c*deriv_m - C_m[ind]*hm)) +
                     decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
                }
          } else { // not s3==1
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                  const double hm = the_f_m_pml[ind];
                  const double hp = the_f_p_pml[ind];
                  const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                  const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
                  the_f_p_pml[ind] += decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
                  the_f_m_pml[ind] += decay_m[ind]*(c*deriv_m - C_m[ind]*hm);
                  the_f[ind] += (decay_m[ind]*(c*deriv_m - C_m[ind]*hm)) +
                     decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
                }
          }
        }
      } else { // not have_p_pml
        if (n3==1) {
          if (s2==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
                const double hm = the_f_m_pml[ind];
                const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
                the_f_m_pml[ind] += decay_m[ind]*(c*deriv_m - C_m[ind]*hm);
                the_f[ind] += (decay_m[ind]*(c*deriv_m - C_m[ind]*hm)) +
                   c*m_deriv_p;
              }
          } else { // not s2==1
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
                const double hm = the_f_m_pml[ind];
                const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
                the_f_m_pml[ind] += decay_m[ind]*(c*deriv_m - C_m[ind]*hm);
                the_f[ind] += (decay_m[ind]*(c*deriv_m - C_m[ind]*hm)) +
                   c*m_deriv_p;
              }
          }
        } else { // not n3==1
          if (s3==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                  const double hm = the_f_m_pml[ind];
                  const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                  const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
                  the_f_m_pml[ind] += decay_m[ind]*(c*deriv_m - C_m[ind]*hm);
                  the_f[ind] += (decay_m[ind]*(c*deriv_m - C_m[ind]*hm)) +
                     c*m_deriv_p;
                }
          } else { // not s3==1
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                  const double hm = the_f_m_pml[ind];
                  const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                  const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
                  the_f_m_pml[ind] += decay_m[ind]*(c*deriv_m - C_m[ind]*hm);
                  the_f[ind] += (decay_m[ind]*(c*deriv_m - C_m[ind]*hm)) +
                     c*m_deriv_p;
                }
          }
        }
      }
    } else { // not have_m_pml
      if (have_p_pml) {
        if (n3==1) {
          if (s2==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
                const double hp = the_f_p_pml[ind];
                const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
                the_f_p_pml[ind] += decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
                the_f[ind] += c*deriv_m + decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
              }
          } else { // not s2==1
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
                const double hp = the_f_p_pml[ind];
                const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
                the_f_p_pml[ind] += decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
                the_f[ind] += c*deriv_m + decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
              }
          }
        } else { // not n3==1
          if (s3==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                  const double hp = the_f_p_pml[ind];
                  const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                  const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
                  the_f_p_pml[ind] += decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
                  the_f[ind] += c*deriv_m + decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
                }
          } else { // not s3==1
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                  const double hp = the_f_p_pml[ind];
                  const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                  const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
                  the_f_p_pml[ind] += decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
                  the_f[ind] += c*deriv_m + decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
                }
          }
        }
      } else { // not have_p_pml
        if (n3==1) {
          if (s2==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
                const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
                the_f[ind] += c*(deriv_m + m_deriv_p);
              }
          } else { // not s2==1
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
                const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
                the_f[ind] += c*(deriv_m + m_deriv_p);
              }
          }
        } else { // not n3==1
          if (s3==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                  const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                  const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
                  the_f[ind] += c*(deriv_m + m_deriv_p);
                }
          } else { // not s3==1
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                  const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                  const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
                  the_f[ind] += c*(deriv_m + m_deriv_p);
                }
          }
        }
      }
    }
  } else { // not have_p
    if (have_m_pml) {
      if (n3==1) {
        if (s2==1) {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
              const double hm = the_f[ind];
              const double deriv_m = f_m[ind]-f_m[ind-stride_m];
              the_f[ind] += decay_m[ind]*(c*deriv_m - C_m[ind]*hm);
            }
        } else { // not s2==1
          for (int i1=0; i1<n1; i1++)
            for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
              const double hm = the_f[ind];
              const double deriv_m = f_m[ind]-f_m[ind-stride_m];
              the_f[ind] += decay_m[ind]*(c*deriv_m - C_m[ind]*hm);
            }
        }
      } else { // not n3==1
        if (s3==1) {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0; i2<n2; i2++)
              for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                const double hm = the_f[ind];
                const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                the_f[ind] += decay_m[ind]*(c*deriv_m - C_m[ind]*hm);
              }
        } else { // not s3==1
          for (int i1=0; i1<n1; i1++)
            for (int i2=0; i2<n2; i2++)
              for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                const double hm = the_f[ind];
                const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                the_f[ind] += decay_m[ind]*(c*deriv_m - C_m[ind]*hm);
              }
        }
      }
    } else { // not have_m_pml
      if (n3==1) {
        if (s2==1) {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
              const double deriv_m = f_m[ind]-f_m[ind-stride_m];
              the_f[ind] += c*deriv_m;
            }
        } else { // not s2==1
          for (int i1=0; i1<n1; i1++)
            for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
              const double deriv_m = f_m[ind]-f_m[ind-stride_m];
              the_f[ind] += c*deriv_m;
            }
        }
      } else { // not n3==1
        if (s3==1) {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0; i2<n2; i2++)
              for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                the_f[ind] += c*deriv_m;
              }
        } else { // not s3==1
          for (int i1=0; i1<n1; i1++)
            for (int i2=0; i2<n2; i2++)
              for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                const double deriv_m = f_m[ind]-f_m[ind-stride_m];
                the_f[ind] += c*deriv_m;
              }
        }
      }
    }
  }
} else { // not have_m
  if (have_p_pml) {
    if (n3==1) {
      if (s2==1) {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
            const double hp = the_f[ind];
            const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
            the_f[ind] += decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
          }
      } else { // not s2==1
        for (int i1=0; i1<n1; i1++)
          for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
            const double hp = the_f[ind];
            const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
            the_f[ind] += decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
          }
      }
    } else { // not n3==1
      if (s3==1) {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0; i2<n2; i2++)
            for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
              const double hp = the_f[ind];
              const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
              the_f[ind] += decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
            }
      } else { // not s3==1
        for (int i1=0; i1<n1; i1++)
          for (int i2=0; i2<n2; i2++)
            for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
              const double hp = the_f[ind];
              const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
              the_f[ind] += decay_p[ind]*(c*m_deriv_p - C_p[ind]*hp);
            }
      }
    }
  } else { // not have_p_pml
    if (n3==1) {
      if (s2==1) {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
            const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
            the_f[ind] += c*m_deriv_p;
          }
      } else { // not s2==1
        for (int i1=0; i1<n1; i1++)
          for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
            const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
            the_f[ind] += c*m_deriv_p;
          }
      }
    } else { // not n3==1
      if (s3==1) {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0; i2<n2; i2++)
            for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
              const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
              the_f[ind] += c*m_deriv_p;
            }
      } else { // not s3==1
        for (int i1=0; i1<n1; i1++)
          for (int i2=0; i2<n2; i2++)
            for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
              const double m_deriv_p = f_p[ind-stride_p]-f_p[ind];
              the_f[ind] += c*m_deriv_p;
            }
      }
    }
  }
}
