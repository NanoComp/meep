\begin{code}
module YeeLattice ( loop_points, loop_inner,
                    consider_electric_polarizations,
                    loop_fields, ec,
                    loop_other_electric_components, sum_other_components, other_c, other_d,
                  ) where

import List ( (\\) )

import StepGen
\end{code}

\begin{code}
loop_points :: Code -> Code
loop_points = doblock "for (int i=0;i<ntot;i++)"

consider_electric_polarizations x =
    figure_out_fields $ cep ["Er", "Ey", "Ex", "Ep", "Ez"] x
    where cep [] x = x
          cep (c:cs) x = ifelse_ ("f["++c++"][0]") (cep cs x)
                         (declare ("other="++c) False $ declare ("field="++c) False $ cep cs x)

figure_out_fields x = ifelse_ "f[Er][0]" (am_cylindrical x) nc
    where nc = not_cylindrical $ ifelse_ "f[Ey][0]" hey ney
          hey = declare "stride_any_direction[X]==1" False $
                is_translatable ["X","Y"] $
                declare "f[Ex][0]" True $
                ifelse_ "f[Ez][0]" (ifelse_ "stride_any_direction[Z]" am3d am2d)
                                   am2dte
          ney = ifelse_ "f[Ex][0]" am1d am2dtm
          am1d = not_translatable ["X","Y"] $ is_translatable ["Z"] $
                 declare "f[Ey][0]" False $ declare "f[Ez][0]" False $
                 declare "f[Ex][0]" True $
                 primary_direction "Z" $ comment "Am in 1D" x
          am2d = primary_direction "Y" $ comment "Am in 2D" x
          am2dtm = primary_direction "Y" $ declare "f[Ez][0]" True $
                   not_translatable ["Z"] $ comment "Am in 2D TM" x
          am2dte = primary_direction "Y" $ declare "f[Ez][0]" False $
                   comment "Am in 2D TE" x
          am3d = primary_direction "Z" $ comment "Am in 3D" x

not_cylindrical x = declare "f[Er][0]" False $ declare "f[Ep][0]" False $
                    not_translatable ["R"] $ x
am_cylindrical x = not_translatable ["X","Y"] $ is_translatable ["R","Z"] $
    declare "f[Er][0]" True $ declare "f[Ep][0]" True $ declare "f[Ez][0]" True $
    declare "f[Ex][0]" False $ declare "f[Ey][0]" False $
    primary_direction "Z" $ comment "Am in cylindrical coordinates." x

not_translatable [] x = x
not_translatable (d:ds) x =
    declare ("stride_any_direction["++d++"]") False $ not_translatable ds x
is_translatable [] x = x
is_translatable (d:ds) x =
    declare ("stride_any_direction["++d++"]") True $ is_translatable ds x

primary_direction d x =
    declare ("stride_any_direction["++d++"]==1") True $
    not_primary_direction (["X","Y","Z","R"]\\[d]) x
not_primary_direction [] x = x
not_primary_direction (d:ds) x =
    declare ("stride_any_direction["++d++"]==1") False $ not_primary_direction ds x

--loop_fields x = doblock "FOR_ELECTRIC_COMPONENTS(ec) if (f[ec][0])" x
--ec = "ec"
loop_fields x =
    for_any_one_of ["field=Ex","field=Ey","field=Ez","field=Ep","field=Er"] x
ec = "field=Ex" |?| "Ex" |:| "field=Ey" |?| "Ey" |:| "field=Ez" |?| "Ez" |:|
     "field=Er" |?| "Er" |:| "field=Ep" |?| "Ep" |:| "aaack bug!"

sum_other_components e =
    sum_for_any_one_of ["other=Ex","other=Ey","other=Ez","other=Ep","other=Er"] e
loop_other_electric_components x =
    for_any_one_of ["other=Ex","other=Ey","other=Ez","other=Ep","other=Er"] x
other_c = "other=Ex" |?| "Ex" |:| "other=Ey" |?| "Ey" |:| "other=Ez" |?| "Ez" |:|
          "other=Er" |?| "Er" |:| "Ep"
other_d = "other=Ex" |?| "X" |:| "other=Ey" |?| "Y" |:| "other=Ez" |?| "Z" |:|
          "other=Er" |?| "R" |:| "P"

{- The following bit loops over "owned" points -}

stride d = ("stride_any_direction["++d++"]") |?|
           (("stride_any_direction["++d++"]==1")|?|"1"|:|("stride_any_direction["++d++"]"))
           |:| "0"
ind d = ("num_any_direction["++d++"]==1")|?|"0"|:|"i"++d
num d = ("num_any_direction["++d++"]==1")|?|"1"|:|("num_any_direction["++d++"]")

define_i =
    doexp $ "const int i" |=| "yee_idx"
                          |+| ind "X" |*| stride "X" |+| ind "Y" |*| stride "Y"
                          |+| ind "Z" |*| stride "Z" |+| ind "R" |*| stride "R"

loop_inner job =
    consider_directions ["Z","R","X","Y"] $
    lad ["X","Y","R","Z"] job
    where lad [] job = docode [define_i, job]
          lad (d:ds) job = loop_direction d $ lad ds job

consider_directions [] x = x
consider_directions (d:ds) x =
    ifelse_ ("stride_any_direction["++d++"]")
                (whether_or_not ("num_any_direction["++d++"]==1") $
                 ifelse_ ("stride_any_direction["++d++"]==1")
                 (not_primary_direction ds $ special_case d $ the_rest)
                 (special_case d $ the_rest))
                (declare ("stride_any_direction["++d++"]==1") False $
                 declare ("num_any_direction["++d++"]==1") False the_rest)
    where the_rest = consider_directions ds x
          special_case "R" x = not_translatable ["X","Y"] x
          special_case "Y" x = is_translatable ["X"] x
          special_case _ x = x

loop_direction d job =
    ifelse_ ("stride_any_direction["++d++"]")
    (ifelse_ ("num_any_direction["++d++"]==1")
     job
     (doblock ("for (int i"<<d<<"=0; "<<
               "i"<<d<<"<"<<num d<<"; "<<
               (("i"<<d) |+=| stride d)<<")") job))
    job
\end{code}
