\begin{code}
module Complex ( get_cmp_part, cmp, loop_complex ) where

import StepGen
\end{code}

\begin{code}
get_cmp_part num = ("cmp")|?| ("real("<<num<<")") |:| ("imag("<<num<<")")

cmp = ("cmp")|?|"1"|:|"0"
loop_complex job =
    ifelse_ "is_real" realjob (for_true_false "cmp" $ docode [job])
        where realjob = declare "cmp" False $ docode [job]
\end{code}
