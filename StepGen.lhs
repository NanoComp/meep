\begin{code}
module StepGen ( Code, Expression, gencode,
                 doexp, docode, doline,
                 if_, ifelse_, (|?|), (|:|),
                 whether_or_not, declare, for_loop,
                 (|+|), (|-|), (|*|), (|+=|), (<<),
               ) where

import Data.FiniteMap
import Monad ( liftM )
\end{code}

\begin{code}
newtype CC a = CC (FiniteMap String Bool -> (FiniteMap String Bool, a))

instance Monad CC where
    (CC a) >>= b = CC f where f st = let (st', x) = a st
                                         CC bb = b x
                                         in bb st'
    return a = CC $ \st -> (st, a)
type Expression = CC String
type Code = CC [String]

gencode :: Code -> String
gencode (CC code) = unlines $ snd $ code emptyFM
\end{code}

\begin{code}
for_loop :: CODE a => String -> String -> a -> Code
for_loop i num inside =
    ifelse_ (num++"==1") (
      docode [doexp $ "const int "<<i<<" = 0",
              docode inside]
    ) (
      docode [doline $ "for (int "<<i<<"=0; "<<i<<"<"<<num<<"; "<<i<<"++) {",
              liftM (map ("  "++)) (docode inside),
              doline "}"]
    )

if_ :: CODE a => String -> a -> Code
if_ b thendo = ifelseCC b (docode thendo) (return []) $ \thethen _ ->
               if is_empty thethen then return []
               else docode [doline $ "if (" ++ b ++ ") {",
                            return $ map ("  "++) thethen,
                            doline "}"
                           ]

ifelse_ :: (CODE a, CODE b) => String -> a -> b -> Code
ifelse_ b thendo elsedo =
    ifelseCC b (docode thendo) (docode elsedo) $ \thethen theelse ->
    case (is_empty thethen, is_empty theelse) of
    (True,True) -> return []
    (True,_) -> return theelse
    (_,True) -> return thethen
    _ -> docode [doline $ "if (" ++ b ++ ") {",
                 return $ map ("  "++) thethen,
                 doline "} else {",
                 return $ map ("  "++) theelse,
                 doline "}"
                ]
is_empty [] = True
is_empty _ = False
whether_or_not :: CODE a => String -> a -> Code
whether_or_not s x = ifelse_ s x x
(|?|) :: (EXPRESSION a, EXPRESSION b) => String -> a -> b -> Expression
infixl 8 |?|
(|?|) b x y = ifelseCC b (expression x) (expression y) $ \xx yy ->
              return $ "("++b++") ? "++p xx++" : "++p yy
(|:|) :: a -> a
infixr 0 |:|
(|:|) a = a
declare :: CODE a => String -> Bool -> a -> Code
declare s b x = do check <- istrueCC s
                   if check == Just (not b)
                      then return []
                      else do setCC s b
                              docode x
\end{code}

\begin{code}
(<<) :: (EXPRESSION a, EXPRESSION b) => a -> b -> Expression
infixr 8 <<
(<<) x y = do xe <- expression x
              ye <- expression y
              return $ xe++ye

(|-|) :: (EXPRESSION a, EXPRESSION b) => a -> b -> Expression
infixr 5 |-|
(|-|) x y = do xe <- expression x
               ye <- expression y
               if xe == "0"
                  then if ye == "0" then return "0"
                       else return $ " - "++ye
                  else if ye == "0" then return xe
                       else return $ padd xe++" - "++padd ye
(|+|) :: (EXPRESSION a, EXPRESSION b) => a -> b -> Expression
infixr 5 |+|
(|+|) x y = do xe <- expression x
               ye <- expression y
               if xe == "0" then return ye
                  else if ye == "0" then return xe
                       else return $ padd xe++" + "++p ye
(|+=|) :: (EXPRESSION a, EXPRESSION b) => a -> b -> Expression
infixr 4 |+=|
x |+=| y = do xe <- expression x
              ye <- expression y
              if ye == "0" then return ""
                 else return $ xe++" += "++ye
(|*|) :: (EXPRESSION a, EXPRESSION b) => a -> b -> Expression
infixr 6 |*|
(|*|) x y = do xe <- expression x
               ye <- expression y
               case (xe,ye) of
                 ("0",_) -> return "0"
                 (_,"0") -> return "0"
                 ("1",_) -> return ye
                 (_,"1") -> return xe
                 _ -> return $ padd xe++" * "++padd ye
p s | '?' `elem` s = "("++s++")"
p s = s
padd s | '+' `elem` s = "("++s++")"
padd s | '-' `elem` s = "("++s++")"
padd s | '?' `elem` s = "("++s++")"
padd s = s
\end{code}

\begin{code}
class EXPRESSION s where
    expression :: s -> Expression
instance EXPRESSION Expression where
    expression x = x
instance EXPRESSION String where
    expression s = return s

doline :: EXPRESSION a => a -> Code
doline e = do o <- expression e
              return [o]
doexp :: EXPRESSION a => a -> Code
doexp e = do o <- expression e
             if o == "" then return []
                        else return $ linebreak $ o++";"
linebreak s | length (dropWhile (==' ') s) < 70 = [s]
linebreak s = let spaces = takeWhile (==' ') s
                  body = dropWhile (==' ') s in -- (
    case break (==' ') $ drop 55 body of
    (_,"") -> [s]
    ("",_) -> [s]
    (first,second) -> (spaces++take 55 body++first) :
                      linebreak (spaces++"  "++second)
               
doexps :: EXPRESSION a => [a] -> Code
doexps es = do os <- sequence $ map expression es
               return $ map (++";") os
\end{code}

\begin{code}

class CODE s where
    docode :: s -> Code
instance CODE Code where
    docode x = x
instance CODE a => CODE [a] where
    docode co = do c <- sequence $ map docode co
                   return $ concat c
\end{code}

\begin{code}
istrueCC :: String -> CC (Maybe Bool)
istrueCC s = CC $ \st -> (st, lookupFM st s)

setCC :: String -> Bool -> CC ()
setCC s b = CC $ \st -> (addToFM (delFromFM st s) s b, ())

ifelseorCC :: String -> CC a -> CC a -> CC a -> CC a
ifelseorCC s thendo elsedo eitherdo = do istrue <- istrueCC s
                                         case istrue of
                                           Just True -> thendo
                                           Just False -> elsedo
                                           Nothing -> eitherdo

ifelseCC :: Eq a => String -> CC a -> CC a -> (a -> a -> CC a) -> CC a
ifelseCC s thendo elsedo doif =
    do istrue <- istrueCC s
       case istrue of
         Just True -> thendo
         Just False -> elsedo
         Nothing -> do tt <- tryWith s True thendo
                       ee <- tryWith s False elsedo
                       if tt == ee then return tt
                                   else doif tt ee

tryWith :: String -> Bool -> CC a -> CC a
tryWith s b (CC x) = CC $ \st -> (st, snd $ x (addToFM (delFromFM st s) s b))
\end{code}
