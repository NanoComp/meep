import StepGen
import Monad ( liftM )

main = putStr $ gencode $ job

job = whether_or_not "pol" $ loop_fields $
      docode [
        find_inveps,
        loop_complex $ loop_points $
        docode [
          find_polarization,
          doexp $ "f[ec][cmp][i]" |=|
                      "the_inveps[i]" |*| ("f[dc][cmp][i]" |-| the_polarization)
        ]
      ]

the_polarization = ("pol")|?|"the_polarization"|:|"0"

find_inveps =
    doexp "const double *the_inveps = ma->inveps[ec][component_direction(ec)]"

find_polarization =
    if_ "pol" $
    docode [doexp "double the_polarization = 0.0",
            doexp "FOR_POLARIZATIONS(pol, p) the_polarization += p->P[ec][cmp][i]"]

loop_fields job = docode [doline "FOR_E_AND_D(ec,dc) if (f[ec][0]) {",
                          indent job,
                          doline "}"]

loop_complex job = docode [doline "DOCMP {",
                           indent job,
                           doline "}"]

loop_points job = docode [doline "for (int i=0;i<ntot;i++) {",
                          indent job,
                          doline "}"]

indent = liftM (map ("  "++))
