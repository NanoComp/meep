import System

main = do
  args <- getArgs
  if length args /= 1
     then exitWith $ ExitFailure 1
     else return ()
  cs <- readFile $ head args
  putStr $ unlines $ ignore $ lines cs

preproc :: [String] -> [String]
preproc ("\\end{verbatim}":ss) = ignore ss
preproc ("\\end{comment}":ss) = ignore ss
preproc (s:ss) = s : preproc ss
preproc [] = []

ignore :: [String] -> [String]
ignore ("\\begin{verbatim}":ss) = preproc ss
ignore ("\\begin{comment}":ss) = preproc ss
ignore (_:ss) = ignore ss
ignore [] = []
