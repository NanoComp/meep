import StepGen
import Complex
import YeeLattice

main = putStr $ gencode $ job

job = consider_electric_polarizations $ docode [
    loop_electric_fields $ finish_polarizations,
    swap_polarizations
  ]

swap_polarizations = docode [
    doline "// The polarizations got switched...",
    doexp "polarization *temp = olpol",
    doexp "olpol = pol",
    doexp "pol = temp"
  ]

finish_polarizations =
    whether_or_not "is_real" $ loop_new_and_old_polarizations $
    docode [doexp "const double om_psqr = fabs(op->pb->saturated_sigma)",
            doexp "const double oo_ep_om_psqr = 0.5/(om_psqr)",
            doexp "const double oo_e_sat = 1.0/op->pb->energy_saturation",
            doexp "const double om_sqr = op->pb->omeganot*op->pb->omeganot",
            doexp "const double g = op->pb->gamma",
            doexp "const double funinv = 1.0/(1+0.5*g)",
            ifelse_ "np->saturation_factor"
            (docode [loop_points $ loop_complex half_step_polarization_energy,
                     loop_inner $ step_saturable_polarization])
            (loop_points $
             docode [loop_complex half_step_polarization_energy,
                     loop_complex step_polarization_itself])
           ]

{- Half-step polarization energy.

The energy change associated with a polarization field is equal to dP*E.
This means E must be known at H time.  To acheive this we do the update in
two steps, once with E before it is updated, and once after (with the same
dP, of course).

-}

half_step_polarization_energy =
    doexp $ "np->energy[ec][i] += (0.5)*(np->P[ec]["<<cmp<<"][i] - "<<
              "op->P[ec]["<<cmp<<"][i])*f[ec]["<<cmp<<"][i]"

step_saturable_polarization =
  docode
  [doexp $ "const double energy_here" |=|
               sum_over_components (\c i-> "np->energy["<<c<<"]["<<i<<"]"),
   doexp $ "const double P_sqr" |=|
               sum_over_components (\c i-> sqr_complex $
               ("np->P["<<c<<"]["<<cmp<<"]["<<i<<"]") |+| ("op->P["<<c<<"]["<<cmp<<"]["<<i<<"]")),
   doexp $ "const double Pdot_sqr" |=|
               sum_over_components (\c i-> sqr_complex $
               ("np->P["<<c<<"]["<<cmp<<"]["<<i<<"]") |-| ("op->P["<<c<<"]["<<cmp<<"]["<<i<<"]")),
   doexp $ "const double energy_P_here" |=| ("om_sqr" |*| "P_sqr" |+|
                                             "(1.0/(Courant*Courant))" |*| "Pdot_sqr")|*|"oo_ep_om_psqr",
   doexp $ "np->s[ec][i] = -(energy_here - energy_P_here)*(om_psqr*oo_e_sat)",
   loop_complex $ step_polarization_itself
  ]

step_polarization_itself =
    doexp $ "op->P[ec]["<<cmp<<"][i] = funinv*((2-om_sqr)*np->P[ec]["<<cmp<<"][i] + "<<
            "(0.5*g-1)*op->P[ec]["<<cmp<<"][i] + np->s[ec][i]*f[ec]["<<cmp<<"][i])"

{- The following code is for a stochastically bouncing polarization. -}
{-
stochastically_step_polarization =
    doblock "" [doexp $ "const double velA" |=| new_p |-| old_p,
                doexp $ "const double impulse = -velA*g + 0.0*gaussian()",
                doexp $ "const double velB" |=| ("impulse" |+| new_p |-| old_p)
                                                |-| "om_sqr"|*| new_p,
                doexp $ "const double p_A" |=| "0.5"|*|(old_p |+| new_p),
                doexp $ "const double p_B" |=| new_p |+| "0.5"|*|"velB",
                doexp $ "const double old_energy = velA*velA + (om_sqr)*p_A*p_A",
                doexp $ "const double new_energy = velB*velB + (om_sqr)*p_B*p_B",
                doexp $ "const double thermal_energy = 0;//np->temperature[ec][i]*exponential()",
                ifelse_ "new_energy < old_energy + thermal_energy"
                (doexp $ old_p |=| "2-om_sqr"|*| new_p |+| "impulse" |+|
                                   ("np->s[ec][i]*f[ec]["<<cmp<<"][i]") |-| old_p)
                (doexp $ old_p |=| "2-om_sqr"|*| new_p |+|
                                   ("np->s[ec][i]*f[ec]["<<cmp<<"][i]") |-| old_p)
               ]

guessed_p = new_p |+| "impulse"
new_p = "np->P[ec]["<<cmp<<"][i]"
old_p = "op->P[ec]["<<cmp<<"][i]"
-}

{- Stuff below is more sort of general-use functions -}

loop_polarizations job =
    if_ "pol" $ doblock "for (polarization *p = pol; p; p = p->next)" job
loop_new_and_old_polarizations job = if_ "pol" $ doblock
    "for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next)" job