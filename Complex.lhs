\begin{code}
module Complex ( get_cmp_part, cmp, loop_complex, sqr_complex ) where

import StepGen
\end{code}

\begin{code}
get_cmp_part num = ("cmp")|?| ("imag("<<num<<")") |:| ("real("<<num<<")")

cmp = ("is_real")|?|"0"|:|(("cmp")|?|"1"|:|"0")
loop_complex job =
    ifelse_ "is_real" realjob (for_true_false "cmp" $ docode [job])
        where realjob = declare "cmp" False $ docode [job]

sqr_complex :: Expression -> Expression
sqr_complex e = ("is_real") |?| (e |*| e)
                            |:| (sum_true_false "cmp" $ e |*| e)
\end{code}
