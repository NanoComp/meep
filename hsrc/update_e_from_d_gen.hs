import StepGen
import Complex
import YeeLattice
import System ( getArgs )

main = do args <- getArgs
          putStr $ gencode $ job args

job ["prepare"]
    = if_ "have_d_minus_p" $ whether_or_not "pol" $
      consider_electric_polarizations $
      doblock "FOR_E_AND_D(ec,dc) if (f[ec][0])" $
      docode
      [
       prepare_polarizations,
       for_complex $ calc_d_minus_p
      ]

job ["sources"]
    = if_ "have_d_minus_p" $ if_ "e_sources" $
      consider_electric_polarizations $ calc_d_minus_p_sources

job ["update"] = whether_or_not "have_d_minus_p" $ consider_electric_polarizations $
                 loop_electric_fields $ regardless_of_inveps $ update_e_from_d_minus_p

job _ = error "Must provide one argument:  prepare, sources or update_e"

d_minus_p c i = ("have_d_minus_p")|?|("d_minus_p["<<c<<"]["<<cmp<<"]["<<i<<"]")
              |:|("f[(component)("<<c<<"+10)]["<<cmp<<"]["<<i<<"]")

{- Here is where we compute the polarization -}

calc_d_minus_p = if_ "have_d_minus_p" $
    docode [
      -- First we look at the contribution from polarization fields.
      loop_points $ docode [
        doexp $ d_minus_p "ec" "i" |=| "f[dc]["<<cmp<<"][i]",
        loop_polarizations $ doexp $ d_minus_p "ec" "i" |-=| "p->P[ec]["<<cmp<<"][i]"]
    ]

calc_d_minus_p_sources = if_ "have_d_minus_p" $
    -- The following code calculates the polarization from sources.
    loop_sources "e_sources" "sv" $
      doblock "if (f[sv->c][0])" $ 
        doblock "for (int j=0; j<sv->npts; j++)" $ 
	  docode
          [
	   doexp "const complex<double> A = sv->dipole(j)",
           for_complex $
	     doexp $ d_minus_p "sv->c" "sv->index[j]" |-=| get_cmp_part "A"
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
    ifelse_ ("s->kerr[ec]") update_e_from_d_minus_p_nonlinear update_e_from_d_minus_p_linear

update_e_from_d_minus_p_linear =
    whether_or_not "p_here" $ for_complex $ loop_inner $
     doexp $ ("f[ec]["<<cmp<<"][i]") |=|
         sum_over_components_with_prefactor (\c -> inveps c "i") (\c i-> d_minus_p c i)

update_e_from_d_minus_p_nonlinear =
    whether_or_not "p_here" $ whether_or_not "is_real" $ loop_inner $ 
    docode
    [
     doexp $ "const double dsqr_here" |=|
         sum_over_components (\c i-> sqr_complex (d_minus_p c i)),
     doexp $ "const double frac_change" |=| "calc_nonlinear_inveps(dsqr_here, s->inveps[ec][d_ec][i], s->kerr[ec][i])",
     loop_complex $ doexp $ ("f[ec]["<<cmp<<"][i]") |=| "frac_change" |*|
         sum_over_components_with_prefactor (\c -> inveps c "i") (\c i-> d_minus_p c i)
    ]

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
    if_ start $ doblock ("for (src_vol *"++svar++" = "++start++"; "++
                               svar++"; "++svar++" = "++svar++"->next)") job
