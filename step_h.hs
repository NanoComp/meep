import StepGen
import Monad ( liftM )

main = putStr $ gencode $ job

job = in_any_case $ loop_over $
      docode [store_h_minus,
              store_h_plus,
              update_f_pml,
              update_f
             ]

store_h_minus = if_ "have_m_pml" (
                  ifelse_ "have_p" (
                    ifelse_ "have_p_pml" (
                      doexp "const double hm = the_f[ind] - the_f_pml[ind]"
                    ) (
                      doexp "const double hm = the_f_pml[ind]"
                    )
                  ) (
                    doexp "const double hm = the_f[ind]"
                  )
                )
store_h_plus = if_ "have_p_pml" (
                 ifelse_ "have_m" (
                   doexp "const double hp = the_f_pml[ind]"
                 ) (
                   doexp "const double hp = the_f[ind]"
                 )
               )
update_f_pml :: Code
update_f_pml = ifelse_ "have_p_pml" (
                 if_ "have_m" (
                   doexp ("the_f_pml[ind]" |+=| p_update)
                 )
               ) (
                 if_ "have_m_pml" (
                   if_ "have_p" (
                     doexp ("the_f_pml[ind]" |+=| m_update)
                   )
                 )
               )

update_f :: Code
update_f = doexp $ "the_f[ind]" |+=| m_p_update

m_update = decay_m |*| ("c" |*| deriv_m |-| sig_m |*| "hm")
p_update = decay_p |*| ("c" |*| m_deriv_p |-| sig_p |*| "hp")
m_p_update = 
    ("have_p_pml")
        |?| (m_update |+| p_update)
        |:| ("have_m_pml")
            |?| (m_update |+| p_update)
            |:| ("c" |*| (deriv_m |+| m_deriv_p))

decay_p = ("have_p_pml") |?| "decay_p[ind]" |:| "1"
decay_m = ("have_m_pml") |?| "decay_m[ind]" |:| "1"
sig_m = ("have_m_pml") |?| "C_m[ind]" |:| "0"
sig_p = ("have_p_pml") |?| "C_p[ind]" |:| "0"
deriv_m = ("have_m") |?| "f_m[ind]-f_m[ind-stride_m]" |:| "0"
m_deriv_p = ("have_p") |?| "f_p[ind-stride_p]-f_p[ind]" |:| "0"

if_have_p_else x y = ifelse_ "have_p" x
                     (declare "have_p_pml" False $ declare "have_m" True y)
if_have_m_else x y = ifelse_ "have_m" x
                     (declare "have_m_pml" False $ declare "have_p" True y)
any_have_p x = if_have_p_else x x
any_have_m x = if_have_m_else x x

in_any_case x = any_have_m $ any_have_p $
                whether_or_not "have_m_pml" $
                whether_or_not "have_p_pml" x

loop_over inside =
    whether_or_not "n3==1" $
    whether_or_not "s2==1" $
    whether_or_not "s3==1" $
    ifelse_ "n3==1" (
      docode [doline "for (int i1=0; i1<n1; i1++)",
              doline $ "  for (int i2=0, ind=i1*s1; i2<n2; i2++,ind+="<<s2<<") {",
              liftM (map ("    "++)) inside,
              doline "  }"]
    ) (
      docode [doline "for (int i1=0; i1<n1; i1++)",
              doline "  for (int i2=0; i2<n2; i2++)",
              doline $ "    for (int i3=0, ind=i1*s1+i2*s2; i3<n3;"<<
                               " i3++, ind+="<<s3<<") {",
              liftM (map ("      "++)) inside,
              doline "    }"]
    )

s2 = ("s2==1")|?|"1"|:|"s2"
s3 = ("s3==1")|?|"1"|:|"s3"
