import StepGen
import Monad ( liftM )

main = putStr $ gencode $ job

job = with_d_minus_p $ whether_or_not "e_sources" $
      loop_fields $ docode [
        prepare_polarizations,
        loop_complex $ docode [
          calc_d_minus_p,
          update_e_from_d_minus_p
        ]
      ]

d_minus_p = ("have_nonzero_polarization")|?|"d_minus_p[i]"|:|"f[dc][cmp][i]"

{- Here is where we compute the polarization -}

with_d_minus_p job =
    with_or_withot_polarization $ ifelse_ "have_nonzero_polarization" (
      docode [doexp "double *d_minus_p = new double[ntot]",
              job,
              doexp "delete[] d_minus_p"]
    ) (
      job
    )

with_or_withot_polarization job =
    ifelse_ "pol" p_job $ ifelse_ "e_sources" p_job $ nop_job
    where p_job = declare "have_nonzero_polarization" True job
          nop_job = declare "have_nonzero_polarization" False job

calc_d_minus_p = if_ "have_nonzero_polarization" $
    docode [
      -- First we look at the contribution from polarization fields.
      loop_points $ docode [
        doexp $ "d_minus_p[i]" |=| "f[dc][cmp][i]",
        loop_polarizations $ doexp $ d_minus_p |-=| "p->P[ec][cmp][i]"],
      -- The following code calculates the polarization from sources.
      loop_sources "e_sources" "s" $
        doexp $ "d_minus_p[s->i]" |-=| get_cmp_part "s->get_dipole_now()*s->A[ec]"
    ]

{- Half-step polarization energy.

The energy change associated with a polarization field is equal to dP*E.
This means E must be known at H time.  To acheive this we do the update in
two steps, once with E before it is updated, and once after (with the same
dP, of course).

-}

prepare_polarizations =
    loop_points $ loop_new_and_old_polarizations $
    docode [prepare_polarization_energy,
            loop_complex half_step_polarization_energy]

prepare_polarization_energy = doexp "np->energy[ec][i] = op->energy[ec][i]"

half_step_polarization_energy =
    doexp "np->energy[ec][i] += 0.5*(np->P[ec][cmp][i] - op->P[ec][cmp][i]) * f[ec][cmp][i]"

{- Here is where we compute E from D - P -}

update_e_from_d_minus_p =
    whether_or_not "p_here" $
      docode [
        --doexp "const double *the_inveps = ma->inveps[ec][component_direction(ec)]",
        --loop_points $ doexp $ "f[ec][cmp][i]" |=| "the_inveps[i]" |*| d_minus_p
        loop_points $ doexp $ "f[ec][cmp][i]" |=|
               "ma->inveps[ec][component_direction(ec)][i]" |*| d_minus_p
      ]

{- Stuff below is more sort of general-use functions -}

get_cmp_part num = ("cmp==0")|?| ("real("<<num<<")") |:| ("imag("<<num<<")")

loop_polarizations job =
    if_ "pol" $ doblock "for (polarization *p = pol; p; p = p->next)" job
loop_new_and_old_polarizations job = if_ "pol" $ doblock
    "for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next)" job
loop_sources start svar job =
    if_ start $ doblock ("for (src *"++svar++" = "++start++"; "++
                               svar++"; "++svar++" = "++svar++"->next)") job
loop_fields = doblock "FOR_E_AND_D(ec,dc) if (f[ec][0])"
loop_complex = doblock "DOCMP"
--loop_complex job =
--    ifelse_ "is_real" realjob (for_true_false "cmp==0" $ docode [initcmp, job])
--        where realjob = declare "cmp==0" True $
--                        docode [initcmp, job]
--              initcmp = doexp $ "const int cmp" |=| ("cmp==0")|?|"0"|:|"1";
loop_points = doblock "for (int i=0;i<ntot;i++)"
