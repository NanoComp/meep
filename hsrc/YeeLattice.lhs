\begin{code}
module YeeLattice ( loop_points, loop_inner,
                    consider_electric_polarizations, consider_all_directions,
                    loop_electric_fields, sum_over_components,
                    sum_over_components_with_prefactor,
                    using_symmetry,
                    sum_D_function_on_eps_point,
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
          cep (c:cs) x = whether_or_not ("f["++c++"][0]") $ cep cs x

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
                 declare "num_any_direction[Z]==1" False $
                 primary_direction "Z" $
                 declare "1D" True $ comment "Am in 1D" $ x
          am2d = primary_direction "Y" $
                 declare "2D" True $ comment "Am in 2D" x
          am2dtm = primary_direction "Y" $ declare "f[Ez][0]" True $
                   not_translatable ["Z"] $
                   declare "2DTM" True $ comment "Am in 2D TM" x
          am2dte = primary_direction "Y" $ declare "f[Ez][0]" False $
                   not_translatable ["Z"] $
                   declare "2DTE" True $ comment "Am in 2D TE" x
          am3d = primary_direction "Z" $
                 declare "3D" True $ comment "Am in 3D" x

not_cylindrical x = declare "f[Er][0]" False $ declare "f[Ep][0]" False $
                    not_translatable ["R"] $ x
am_cylindrical x = not_translatable ["X","Y"] $ is_translatable ["R","Z"] $
    declare "f[Er][0]" True $ declare "f[Ep][0]" True $ declare "f[Ez][0]" True $
    declare "f[Ex][0]" False $ declare "f[Ey][0]" False $
    primary_direction "Z" $
    declare "CYLINDRICAL" True $ comment "Am in cylindrical coordinates." x

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

loop_electric_fields x =
    doblock "FOR_E_AND_D(ec,dc) if (f[ec][0])" $
    docode [doexp $ "const int yee_idx = v.yee_index(ec)",
            doexp $ "const int d_ec = component_direction(ec)",
            doexp $ "const int s_ec = stride_any_direction[d_ec]",
            using_symmetry $ define_others,
            x]
    where dc n "CYLINDRICAL" = "(direction)(((d_ec-2)+"++show n++")%3+2)"
          dc n "2DTE" = "(direction)((d_ec+"++show n++")%2)"
          dc n _ = "(direction)((d_ec+"++show n++")%3)"
          define_others "1D" = doexp ""
          define_others "2DTE" = define_other 1
          define_others "2DTM" = doexp ""
          define_others _ = docode [define_other 1, define_other 2]
          define_other n =
              docode [doexp $ "const direction d_"<<show n<<
                              " = " << with_symmetry (dc n),
                      doexp $ "const component ec_"<<show n<<
                              " = direction_component(ec,d_" << show n<<")",
                      doexp $ "const component dc_"<<show n<<
                              " = direction_component(dc,d_" << show n<<")",
                      doexp $ "const int s_"++show n++
                              " = stride_any_direction[d_"++show n++"]"]

sum_over_components :: (String -> String -> Expression) -> Expression
sum_over_components = sum_over_gen_component "ec"

sum_over_gen_component :: String -> (String -> String -> Expression) -> Expression
sum_over_gen_component c e = with_symmetry sumit
    where sumit "1D" = e c "i"
          sumit "2DTM" = e c "i"
          sumit "2DTE" = e c "i" |+|
              "0.25"|*|(e c_1 "i" |+|
                        e c_1 "i-s_1" |+|
                        e c_1 "i+s_ec" |+|
                        e c_1 "i-s_1+s_ec")
          sumit _ =  e c "i" |+|
              "0.25"|*|(e c_1 "i" |+|
                        e c_1 "i-s_1" |+|
                        e c_1 "i+s_ec" |+|
                        e c_1 "i-s_1+s_ec") |+|
              "0.25"|*|(e c_2 "i" |+|
                        e c_2 "i-s_2" |+|
                        e c_2 "i+s_ec" |+|
                        e c_2 "i-s_2+s_ec")
          c_1 = c ++ "_1"
          c_2 = c ++ "_2"

sum_over_components_with_prefactor :: (String -> Expression)
                                   -> (String -> String -> Expression) -> Expression
sum_over_components_with_prefactor pre e = with_symmetry sumit
    where sumit "1D" = pre "ec" |*| e "ec" "i"
          sumit "2DTM" = pre "ec" |*| e "ec" "i"
          sumit "2DTE" = pre "ec" |*| e "ec" "i" |+|
              "0.25"|*| pre "ec_1" |*| (e "ec_1" "i" |+|
                                        e "ec_1" "i-s_1" |+|
                                        e "ec_1" "i+s_ec" |+|
                                        e "ec_1" "i-s_1+s_ec")
          sumit _ =  pre "ec" |*| e "ec" "i" |+|
              "0.25"|*| pre "ec_1" |*| (e "ec_1" "i" |+|
                                        e "ec_1" "i-s_1" |+|
                                        e "ec_1" "i+s_ec" |+|
                                        e "ec_1" "i-s_1+s_ec") |+|
              "0.25"|*| pre "ec_2" |*| (e "ec_2" "i" |+|
                                        e "ec_2" "i-s_2" |+|
                                        e "ec_2" "i+s_ec" |+|
                                        e "ec_2" "i-s_2+s_ec")

{- The following bit loops over "owned" points -}

declare_strides :: Code -> Code
declare_strides job = using_symmetry ds
    where ds "1D" = docode [doexp $ s "Z", job]
          ds "3D" = docode $ map (doexp . s) ["X","Y","Z"] ++ [job]
          ds "CYLINDRICAL" = docode $ map (doexp . s) ["R","Z"] ++ [job]
          ds _ = docode $ map (doexp . s) ["X","Y"] ++ [job]
          s d = ("const int s" ++ d) |=| actual_stride d

actual_stride :: String -> Expression
actual_stride d = ("stride_any_direction["++d++"]") |?|
                  (("stride_any_direction["++d++"]==1")|?|"1"
                   |:|("stride_any_direction["++d++"]"))
                      |:| "0"
stride d = ("stride_any_direction["++d++"]") |?|
           (("stride_any_direction["++d++"]==1")|?|"1"
            |:|("s" ++ d)) |:| "0"

ind d = casedef ["num_any_direction["++d++"]==1"] (\_ -> "0") $ "i"++d
num d = casedef ["num_any_direction["++d++"]==1"] (\_ -> "1") $ "num_any_direction["++d++"]"

define_i =
    doexp $ "const int i" |=| "yee_idx"
                          |+| ind "X" |*| stride "X" |+| ind "Y" |*| stride "Y"
                          |+| ind "Z" |*| stride "Z" |+| ind "R" |*| stride "R"

loop_inner :: Code -> Code
loop_inner job = consider_all_directions $ declare_strides $ lad ["X","Y","R","Z"]
    where lad [] = docode [define_i, job]
          lad (d:ds) = loop_direction d $ lad ds

consider_all_directions job =
    look_for_short $ consider_directions ["Z","R","X","Y"] job
    where look_for_short job =
              ifelse_ "num_any_direction[Z]==1"
              (declare "num_any_direction[X]==1" False $
               declare "num_any_direction[Y]==1" False $
               declare "num_any_direction[R]==1" False $ job) $
              ifelse_ "num_any_direction[X]==1"
              (declare "num_any_direction[Y]==1" False $
               declare "num_any_direction[R]==1" False $ job) $
              ifelse_ "num_any_direction[Y]==1"
              (declare "num_any_direction[R]==1" False $ job) $
              whether_or_not "num_any_direction[R]==1" job

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
    (ifdefelse ("num_any_direction["++d++"]==1")
     job
     (doblock ("for (int i"<<d<<"=0; "<<
               "i"<<d<<"<"<<num d<<"; "<<
               ("i"<<d<<"++")<<")") job))
    job

sum_D_function_on_eps_point :: (String -> String -> Expression) -> Expression
sum_D_function_on_eps_point e = with_symmetry sumit
    where sumit "1D" = "0.5" |*| two_D "Dx" "Z"
          sumit "2DTE" = "0.5" |*| (two_D "Dy" "X" |+| two_D "Dx" "Y")
          sumit "2DTM" = "0.25" |*| (four_D "Dz" "X" "Y")
          sumit "CYLINDRICAL" = "0.5" |*| (two_D "Dz" "R" |+| two_D "Dr" "Z") |+| 
                                "0.25" |*| (four_D "Dp" "R" "Z")
          sumit "2D" = "0.5" |*| (two_D "Dy" "X" |+| two_D "Dx" "Y")
                     |+| "0.25" |*| four_D "Dz" "X" "Y"
          sumit "3D" = "0.25" |*| (four_D "Dz" "X" "Y" |+|
                                   four_D "Dx" "Y" "Z" |+|
                                   four_D "Dy" "X" "Z")
          two_D f d = e f "i" |+| e f ("i+s" ++ d)
          four_D f d1 d2 = e f "i" |+| e f ("i+s" ++ d1)
                       |+| e f ("i+s" ++ d2) |+| e f ("i+s"++d1++"+s"++d2) 
\end{code}

\begin{code}
with_symmetry e = casedef ["1D","2D","2DTE","2DTM","3D","CYLINDRICAL"] e $
                  error "aaack"
using_symmetry x = docode [ifdef "CYLINDRICAL" $ x "CYLINDRICAL",
                           ifdef "1D" $ x "1D",
                           ifdef "2D" $ x "2D",
                           ifdef "2DTE" $ x "2DTE",
                           ifdef "2DTM" $ x "2DTM",
                           ifdef "3D" $ x "3D"]
\end{code}
