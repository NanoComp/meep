import StepGen
import Complex
import YeeLattice

main = putStr $ gencode $ job

job = with_d_minus_p $ whether_or_not "e_sources" $ consider_electric_polarizations $
      docode [
        loop_electric_fields $ docode [
          prepare_polarizations,
          for_complex $ calc_d_minus_p
        ],
	calc_d_minus_p_sources,
        loop_electric_fields $ regardless_of_inveps $ for_complex $ update_e_from_d_minus_p
      ]

d_minus_p c i = ("have_nonzero_polarization")|?|("d_minus_p["<<c<<"]["<<cmp<<"]["<<i<<"]")
              |:|("f[(component)("<<c<<"+10)]["<<cmp<<"]["<<i<<"]")

{- Here is where we compute the polarization -}

with_d_minus_p job =
    with_or_withot_polarization $ ifelse_ "have_nonzero_polarization" (
      docode [doexp "double *d_minus_p[5][2]",
              for_complex $ doblock "FOR_ELECTRIC_COMPONENTS(ec)" $
                ifelse_ "f[ec]"
                (doexp $ "d_minus_p[ec]["<<cmp<<"] = new double[v.ntot()]")
                (doexp $ "d_minus_p[ec]["<<cmp<<"] = 0"),
              job,
              for_complex $ doblock "FOR_ELECTRIC_COMPONENTS(ec)" $
                doexp $ "delete[] d_minus_p[ec]["<<cmp<<"]"]
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
        doexp $ d_minus_p "ec" "i" |=| "f[dc]["<<cmp<<"][i]",
        loop_polarizations $ doexp $ d_minus_p "ec" "i" |-=| "p->P[ec]["<<cmp<<"][i]"]
    ]

calc_d_minus_p_sources = if_ "have_nonzero_polarization" $
    -- The following code calculates the polarization from sources.
    loop_sources "e_sources" "spt" $
      doblock "if (f[spt->c][0])" $ 
	docode [
	  doexp "const complex<double> A = spt->dipole()",
          for_complex $
	    doexp $ d_minus_p "spt->c" "spt->i" |-=| get_cmp_part "A"
	]

{- Half-step polarization energy.

The energy change associated with a polarization field is equal to dP*E.
This means E must be known at H time.  To acheive this we do the update in
two steps, once with E before it is updated, and once after (with the same
dP, of course).

-}

prepare_polarizations =
    whether_or_not "is_real" $ loop_points $ loop_new_and_old_polarizations $
    docode [doexp "np->energy[ec][i] = op->energy[ec][i]",
            loop_complex half_step_polarization_energy]

half_step_polarization_energy =
    doexp $ "np->energy[ec][i] += (0.5/(4*pi))*(np->P[ec]["<<cmp<<"][i] - op->P[ec]["<<cmp<<"][i])"
              <<" * f[ec]["<<cmp<<"][i]"

{- Here is where we compute E from D - P -}

update_e_from_d_minus_p =
    whether_or_not "p_here" $ loop_inner $
        doexp $ ("f[ec]["<<cmp<<"][i]") |=|
               sum_over_components_with_prefactor (\c -> inveps c "i") (\c i-> d_minus_p c i)

regardless_of_inveps job =
    declare "s->inveps[ec][d_ec]" True $
    using_symmetry consider_inv
    where consider_inv "1D" = job
          consider_inv "2DTM" = job
          consider_inv "2DTE" = whether_or_not "s->inveps[ec][d_1]" job
          consider_inv _ = whether_or_not "s->inveps[ec][d_1]" $
                           whether_or_not "s->inveps[ec][d_2]" job

inveps :: String -> String -> Expression
inveps "ec" i = ("s->inveps[ec][d_ec]")
              |?| ("s->inveps[ec][d_ec]["++i++"]") |:| "0"
inveps "ec_1" i = ("s->inveps[ec][d_1]")
                |?| ("s->inveps[ec][d_1]["++i++"]") |:| "0"
inveps "ec_2" i = ("s->inveps[ec][d_2]")
                |?| ("s->inveps[ec][d_2]["++i++"]") |:| "0"
inveps c i = error $ "inveps can't use "++c

{- Stuff below is more sort of general-use functions -}

loop_polarizations job =
    if_ "pol" $ doblock "for (polarization *p = pol; p; p = p->next)" job
loop_new_and_old_polarizations job = if_ "pol" $ doblock
    "for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next)" job
loop_sources start svar job =
    if_ start $ doblock ("for (src_pt *"++svar++" = "++start++"; "++
                               svar++"; "++svar++" = "++svar++"->next)") job