import StepGen
import Monad ( liftM )

main = putStr $ gencode $ job

job = docode [
    loop_fields $ finish_polarizations,
    swap_polarizations
  ]

{- Half-step polarization energy.

The energy change associated with a polarization field is equal to dP*E.
This means E must be known at H time.  To acheive this we do the update in
two steps, once with E before it is updated, and once after (with the same
dP, of course).

-}

swap_polarizations = docode [
    doline "// The polarizations got switched...",
    doexp "polarization *temp = olpol",
    doexp "olpol = pol",
    doexp "pol = temp"
  ]

finish_polarizations =
    whether_or_not "is_real" $ loop_new_and_old_polarizations $
    docode [prepare_pb_vars,
            whether_or_not "fac" $ whether_or_not "fac > 0" $ loop_points $
            docode [loop_complex half_step_polarization_energy,
                    update_polarization_saturation,
                    loop_complex step_polarization_itself
                   ]
           ]

half_step_polarization_energy =
    doexp $ "np->energy[ec][i] += 0.5*(np->P[ec]["<<cmp<<"][i] - "<<
              "op->P[ec]["<<cmp<<"][i])*f[ec]["<<cmp<<"][i]"

prepare_pb_vars =
    docode [doexp "const double fac = np->saturation_factor",
            doexp "const double g = op->pb->gamma;",
            doexp "const double om = op->pb->omeganot;",
            doexp "const double funinv = 1.0/(1+0.5*g);"]

update_polarization_saturation = if_ "fac" $
    ifelse_ "fac > 0" (doexp "np->s[ec][i] = max(-np->energy[ec][i]*fac, 0.0)")
                      (doexp "np->s[ec][i] = np->energy[ec][i]*fac")

step_polarization_itself =
    doexp $ "op->P[ec]["<<cmp<<"][i] = funinv*((2-om*om)*np->P[ec]["<<cmp<<"][i] + "<<
            "(0.5*g-1)*op->P[ec]["<<cmp<<"][i] + np->s[ec][i]*f[ec]["<<cmp<<"][i])"
{- Stuff below is more sort of general-use functions -}

get_cmp_part num = ("cmp")|?| ("real("<<num<<")") |:| ("imag("<<num<<")")

loop_polarizations job =
    if_ "pol" $ doblock "for (polarization *p = pol; p; p = p->next)" job
loop_new_and_old_polarizations job = if_ "pol" $ doblock
    "for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next)" job
loop_fields = doblock "FOR_E_AND_D(ec,dc) if (f[ec][0])"
cmp = ("cmp")|?|"1"|:|"0"
loop_complex job =
    ifelse_ "is_real" realjob (for_true_false "cmp" $ docode [job])
        where realjob = declare "cmp" False $ docode [job]
loop_points = doblock "for (int i=0;i<ntot;i++)"
