import StepGen
import Monad ( liftM )

main = putStr $ gencode $ job

job = docode [
    loop_fields $ finish_polarizations,
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
    docode [doexp "const double fac = np->saturation_factor",
            doexp "const double g = op->pb->gamma",
            doexp "const double om = op->pb->omeganot",
            doexp "const double invomsqr = 1.0/(om*om)",
            doexp "const double funinv = 1.0/(1+0.5*g)",
            ifelse_ "fac"
            (whether_or_not "fac > 0" $
             docode [loop_points $ loop_complex half_step_polarization_energy,
                     loop_inner $ step_saturable_polarization])
            (loop_points $
             docode [loop_complex half_step_polarization_energy,
                     --loop_complex stochastically_step_polarization
                     loop_complex step_polarization_itself])
           ]

{- Half-step polarization energy.

The energy change associated with a polarization field is equal to dP*E.
This means E must be known at H time.  To acheive this we do the update in
two steps, once with E before it is updated, and once after (with the same
dP, of course).

-}

half_step_polarization_energy =
    doexp $ "np->energy[ec][i] += 0.5*(np->P[ec]["<<cmp<<"][i] - "<<
              "op->P[ec]["<<cmp<<"][i])*f[ec]["<<cmp<<"][i]"

step_saturable_polarization = if_ "fac" $
  docode
  [ifelse_ "fac > 0" (doexp $ "np->s[ec]["<<i<<"] = max(-np->energy[ec]["<<i<<"]*fac, 0.0)")
                     (doexp $ "np->s[ec]["<<i<<"] = np->energy[ec]["<<i<<"]*fac"),
   loop_complex $ doexp $
   "op->P[ec]["<<cmp<<"]["<<i<<"] = funinv*((2-om*om)*np->P[ec]["<<cmp<<"]["<<i<<"] + "<<
   "(0.5*g-1)*op->P[ec]["<<cmp<<"]["<<i<<"] + np->s[ec]["<<i<<"]*f[ec]["<<cmp<<"]["<<i<<"])"
  ]

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

{- The following bit loops over "owned" points -}

stride d = ("stride_any_direction["++d++"]") |?|
           (("stride_any_direction["++d++"]==1")|?|"1"|:|("stride_any_direction["++d++"]"))
           |:| "0"
ind d = ("num_any_direction["++d++"]==1")|?|"0"|:|"i"++d
num d = ("num_any_direction["++d++"]==1")|?|"1"|:|("num_any_direction["++d++"]")

i = ind "X" |*| stride "X" |+| ind "Y" |*| stride "Y" |+|
    ind "Z" |*| stride "Z" |+| ind "R" |*| stride "R"

loop_inner job =
    consider_directions ["Z","R","X","Y"] $
    lad ["X","Y","R","Z"] job
    where lad [] job = job
          lad (d:ds) job = loop_direction d $ lad ds job

consider_directions [] x = x
consider_directions (d:ds) x =
    ifelse_ ("stride_any_direction["++d++"]")
                (whether_or_not ("num_any_direction["++d++"]==1") $
                 ifelse_ ("stride_any_direction["++d++"]==1")
                 (havent_got_stride_one ds $ special_case d $ the_rest)
                 (special_case d $ the_rest))
                (declare ("stride_any_direction["++d++"]==1") False $
                 declare ("num_any_direction["++d++"]==1") False the_rest)
    where the_rest = consider_directions ds x
          special_case "R" x = declare "stride_any_direction[X]" False $
                               declare "stride_any_direction[Y]" False x
          special_case "Y" x = declare "stride_any_direction[X]" True x
          special_case _ x = x

havent_got_stride_one [] x = x
havent_got_stride_one (d:ds) x =
    declare ("stride_any_direction["++d++"]==1") False $ havent_got_stride_one ds x

loop_direction d job =
    ifelse_ ("stride_any_direction["++d++"]")
    (ifelse_ ("num_any_direction["++d++"]==1")
     job
     (doblock ("for (int i"<<d<<"=0; "<<
               "i"<<d<<"<"<<num d<<"; "<<
               (("i"<<d) |+=| stride d)<<")") job))
    job
