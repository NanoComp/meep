if (have_m) {
  if (have_p) {
    if (have_m_pml) {
      if (have_p_pml) {
        if (n3==1) {
          if (s2==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
                const double em = the_f[ind] - the_f_pml[ind];
                const double ep = the_f_pml[ind];
                the_f_pml[ind] += decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind])))
                   - C_p[ind]*ep);
                the_f[ind] += (decay_m[ind]*((c*((f_m[ind]-f_m[ind+stride_m])))
                   - C_m[ind]*em)) + decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind])))
                      - C_p[ind]*ep);
              }
          } else {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
                const double em = the_f[ind] - the_f_pml[ind];
                const double ep = the_f_pml[ind];
                the_f_pml[ind] += decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind])))
                   - C_p[ind]*ep);
                the_f[ind] += (decay_m[ind]*((c*((f_m[ind]-f_m[ind+stride_m])))
                   - C_m[ind]*em)) + decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind])))
                      - C_p[ind]*ep);
              }
          }
        } else {
          if (s3==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                  const double em = the_f[ind] - the_f_pml[ind];
                  const double ep = the_f_pml[ind];
                  the_f_pml[ind] += decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind])))
                     - C_p[ind]*ep);
                  the_f[ind] += (decay_m[ind]*((c*((f_m[ind]-f_m[ind+stride_m])))
                     - C_m[ind]*em)) + decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind])))
                        - C_p[ind]*ep);
                }
          } else {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                  const double em = the_f[ind] - the_f_pml[ind];
                  const double ep = the_f_pml[ind];
                  the_f_pml[ind] += decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind])))
                     - C_p[ind]*ep);
                  the_f[ind] += (decay_m[ind]*((c*((f_m[ind]-f_m[ind+stride_m])))
                     - C_m[ind]*em)) + decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind])))
                        - C_p[ind]*ep);
                }
          }
        }
      } else {
        if (n3==1) {
          if (s2==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
                const double em = the_f_pml[ind];
                the_f_pml[ind] += decay_m[ind]*((c*((f_m[ind]-f_m[ind+stride_m])))
                   - C_m[ind]*em);
                the_f[ind] += (decay_m[ind]*((c*((f_m[ind]-f_m[ind+stride_m])))
                   - C_m[ind]*em)) + inveps[ind]*(c*((f_p[ind+stride_p]-f_p[ind])));
              }
          } else {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
                const double em = the_f_pml[ind];
                the_f_pml[ind] += decay_m[ind]*((c*((f_m[ind]-f_m[ind+stride_m])))
                   - C_m[ind]*em);
                the_f[ind] += (decay_m[ind]*((c*((f_m[ind]-f_m[ind+stride_m])))
                   - C_m[ind]*em)) + inveps[ind]*(c*((f_p[ind+stride_p]-f_p[ind])));
              }
          }
        } else {
          if (s3==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                  const double em = the_f_pml[ind];
                  the_f_pml[ind] += decay_m[ind]*((c*((f_m[ind]-f_m[ind+stride_m])))
                     - C_m[ind]*em);
                  the_f[ind] += (decay_m[ind]*((c*((f_m[ind]-f_m[ind+stride_m])))
                     - C_m[ind]*em)) + inveps[ind]*(c*((f_p[ind+stride_p]-f_p[ind])));
                }
          } else {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                  const double em = the_f_pml[ind];
                  the_f_pml[ind] += decay_m[ind]*((c*((f_m[ind]-f_m[ind+stride_m])))
                     - C_m[ind]*em);
                  the_f[ind] += (decay_m[ind]*((c*((f_m[ind]-f_m[ind+stride_m])))
                     - C_m[ind]*em)) + inveps[ind]*(c*((f_p[ind+stride_p]-f_p[ind])));
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
                const double ep = the_f_pml[ind];
                the_f_pml[ind] += decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind])))
                   - C_p[ind]*ep);
                the_f[ind] += (inveps[ind]*(c*((f_m[ind]-f_m[ind+stride_m]))))
                   + decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind]))) - C_p[ind]*ep);
              }
          } else {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
                const double ep = the_f_pml[ind];
                the_f_pml[ind] += decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind])))
                   - C_p[ind]*ep);
                the_f[ind] += (inveps[ind]*(c*((f_m[ind]-f_m[ind+stride_m]))))
                   + decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind]))) - C_p[ind]*ep);
              }
          }
        } else {
          if (s3==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                  const double ep = the_f_pml[ind];
                  the_f_pml[ind] += decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind])))
                     - C_p[ind]*ep);
                  the_f[ind] += (inveps[ind]*(c*((f_m[ind]-f_m[ind+stride_m]))))
                     + decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind]))) - C_p[ind]*ep);
                }
          } else {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                  const double ep = the_f_pml[ind];
                  the_f_pml[ind] += decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind])))
                     - C_p[ind]*ep);
                  the_f[ind] += (inveps[ind]*(c*((f_m[ind]-f_m[ind+stride_m]))))
                     + decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind]))) - C_p[ind]*ep);
                }
          }
        }
      } else {
        if (n3==1) {
          if (s2==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
                the_f[ind] += c*inveps[ind]*(((f_m[ind]-f_m[ind+stride_m]))
                   + (f_p[ind+stride_p]-f_p[ind]));
              }
          } else {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
                the_f[ind] += c*inveps[ind]*(((f_m[ind]-f_m[ind+stride_m]))
                   + (f_p[ind+stride_p]-f_p[ind]));
              }
          }
        } else {
          if (s3==1) {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                  the_f[ind] += c*inveps[ind]*(((f_m[ind]-f_m[ind+stride_m]))
                     + (f_p[ind+stride_p]-f_p[ind]));
                }
          } else {
            for (int i1=0; i1<n1; i1++)
              for (int i2=0; i2<n2; i2++)
                for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                  the_f[ind] += c*inveps[ind]*(((f_m[ind]-f_m[ind+stride_m]))
                     + (f_p[ind+stride_p]-f_p[ind]));
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
              the_f[ind] += decay_m[ind]*((c*((f_m[ind]-f_m[ind+stride_m])))
                 - C_m[ind]*em);
            }
        } else {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
              const double em = the_f[ind];
              the_f[ind] += decay_m[ind]*((c*((f_m[ind]-f_m[ind+stride_m])))
                 - C_m[ind]*em);
            }
        }
      } else {
        if (s3==1) {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0; i2<n2; i2++)
              for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                const double em = the_f[ind];
                the_f[ind] += decay_m[ind]*((c*((f_m[ind]-f_m[ind+stride_m])))
                   - C_m[ind]*em);
              }
        } else {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0; i2<n2; i2++)
              for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                const double em = the_f[ind];
                the_f[ind] += decay_m[ind]*((c*((f_m[ind]-f_m[ind+stride_m])))
                   - C_m[ind]*em);
              }
        }
      }
    } else {
      if (n3==1) {
        if (s2==1) {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
              the_f[ind] += c*inveps[ind]*((f_m[ind]-f_m[ind+stride_m]));
            }
        } else {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
              the_f[ind] += c*inveps[ind]*((f_m[ind]-f_m[ind+stride_m]));
            }
        }
      } else {
        if (s3==1) {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0; i2<n2; i2++)
              for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
                the_f[ind] += c*inveps[ind]*((f_m[ind]-f_m[ind+stride_m]));
              }
        } else {
          for (int i1=0; i1<n1; i1++)
            for (int i2=0; i2<n2; i2++)
              for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
                the_f[ind] += c*inveps[ind]*((f_m[ind]-f_m[ind+stride_m]));
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
            the_f[ind] += decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind])))
               - C_p[ind]*ep);
          }
      } else {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
            const double ep = the_f[ind];
            the_f[ind] += decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind])))
               - C_p[ind]*ep);
          }
      }
    } else {
      if (s3==1) {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0; i2<n2; i2++)
            for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
              const double ep = the_f[ind];
              the_f[ind] += decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind])))
                 - C_p[ind]*ep);
            }
      } else {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0; i2<n2; i2++)
            for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
              const double ep = the_f[ind];
              the_f[ind] += decay_p[ind]*((c*((f_p[ind+stride_p]-f_p[ind])))
                 - C_p[ind]*ep);
            }
      }
    }
  } else {
    if (n3==1) {
      if (s2==1) {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=1) {
            the_f[ind] += c*inveps[ind]*((f_p[ind+stride_p]-f_p[ind]));
          }
      } else {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+=s2) {
            the_f[ind] += c*inveps[ind]*((f_p[ind+stride_p]-f_p[ind]));
          }
      }
    } else {
      if (s3==1) {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0; i2<n2; i2++)
            for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=1) {
              the_f[ind] += c*inveps[ind]*((f_p[ind+stride_p]-f_p[ind]));
            }
      } else {
        for (int i1=0; i1<n1; i1++)
          for (int i2=0; i2<n2; i2++)
            for (int i3=0, ind=i1*s1+i2*s2; i3<n3; i3++, ind+=s3) {
              the_f[ind] += c*inveps[ind]*((f_p[ind+stride_p]-f_p[ind]));
            }
      }
    }
  }
}
