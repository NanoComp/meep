data Bools = Have { have_m :: Maybe Bool, have_p :: Maybe Bool,
                    have_m_pml :: Maybe Bool, have_p_pml :: Maybe Bool }

main = putStr $ unlines $ loop `given` all_unknown

loop :: B2 [String]
loop = in_any_case
            (loop_over
             (store_h_minus +||+ store_h_plus +||+ update_f_pml +||+ update_f))
store_h_minus = if_have_m_pml
                (if_have_p
                 (if_have_p_pml
                  ["const double hm = the_f[ind] - the_f_pml[ind];"]
                  `else_do`
                  ["const double hm = the_f_pml[ind];"]
                 )
                 `else_do`
                 ["const double hm = the_f[ind];"])
                `else_do` nothing
store_h_plus = if_have_p_pml
               (if_have_m
                ["const double hp = the_f_pml[ind];"]
                `else_do`
                ["const double hp = the_f[ind];"])
               `else_do` nothing
update_f_pml = if_have_any_pml
                 (if_have_p_pml
                  (if_have_m
                   ["the_f_pml[ind]" |+=| p_update +|+ ";"]
                   `else_do` nothing)
                  (if_have_p
                   ["the_f_pml[ind]" |+=| m_update +|+ ";"]
                   `else_do` nothing))
               `else_do` nothing
update_f :: [B2 String]
update_f = ["the_f[ind]" |+=| m_p_update +|+";"]

m_update :: B2 String
p_update :: B2 String
m_update = paren(decay_m |*| paren("c" |*| deriv_m |-| sig_m |*| "hm"))
p_update = paren(decay_p |*| paren("c" |*| m_deriv_p |-| sig_p |*| "hp"))
m_p_update :: B2 String
m_p_update = (m_update |+| p_update) `or_if_no_pml`
             ("c" |*| paren(deriv_m |+| m_deriv_p))

in_any_case job = either_way if_have_m $
                  either_way if_have_p $
                  either_way if_have_p_pml $
                  either_way if_have_m_pml job

decay_p = "decay_p[ind]" `have_p_pml_or` "1"
decay_m = "decay_m[ind]" `have_m_pml_or` "1"
sig_m = "C_m[ind]" `have_m_pml_or` "0"
sig_p = "C_p[ind]" `have_p_pml_or` "0"
deriv_m = "(f_m[ind]-f_m[ind-stride_m])" `have_m_or` "0"
m_deriv_p = "(f_p[ind-stride_p]-f_p[ind])" `have_p_or` "0"

if_have_any_pml :: (Fancy [String] a, Fancy [String] b) =>
                   a -> b -> B2 [String]
if_have_any_pml y n = if_have_m_pml y
                      `else_do` if_have_p_pml y
                                `else_do` n

nothing :: B2 [String]
nothing = B2 (\_ -> [])

(|=|) :: (Fancy String a, Fancy String b) => a -> b -> B2 String
infixr 5 |=|
(|=|) x y = B2 f where f b = case y `given` b of
                             "0" -> ""
                             ystr -> case x `given` b of
                                     xstr -> xstr++" = "++ystr

(|+=|) :: (Fancy String a, Fancy String b) => a -> b -> B2 String
infixr 5 |+=|
(|+=|) x y = B2 f where f b = case y `given` b of
                              "0" -> ""
                              ystr -> case x `given` b of
                                      xstr -> xstr++" += "++ystr

(|*|) :: (Fancy String a, Fancy String b) => a -> b -> B2 String
infixr 8 |*|
(|*|) x y = B2 f where f b = case x `given` b of
                             "0" -> "0"
                             "1" -> y `given` b
                             xstr -> case y `given` b of
                                     "0" -> "0"
                                     "1" -> xstr
                                     ystr -> xstr++"*"++ystr

(|-|) :: (Fancy String a, Fancy String b) => a -> b -> B2 String
infixr 7 |-|
(|-|) x y = B2 f where f b = case x `given` b of
                             "0" -> y `given` b
                             xstr -> case y `given` b of
                                     "0" -> xstr
                                     ystr -> xstr++" - "++ystr

(|+|) :: (Fancy String a, Fancy String b) => a -> b -> B2 String
infixr 7 |+|
(|+|) x y = B2 f where f b = case x `given` b of
                             "0" -> y `given` b
                             xstr -> case y `given` b of
                                     "0" -> xstr
                                     ystr -> xstr++" + "++ystr

paren :: Fancy String f => f -> B2 String
paren x = B2 f where f b = if has_space (x `given` b)
                           then "("++ x `given` b ++")"
                           else x `given` b
has_space [] = False
has_space (' ':_) = True
has_space (_:cs) = has_space cs

(+|+) :: (Fancy String f, Fancy String g) => f -> g -> B2 String
infixr 3 +|+
(+|+) x y = B2 f where f b = (x `given` b) ++ (y `given` b)

(+||+) :: (Fancy [String] f, Fancy [String] g) => f -> g -> B2 [String]
infixr 3 +||+
(+||+) x y = B2 f where f b = (x `given` b) ++ (y `given` b)

all_unknown = Have Nothing Nothing Nothing Nothing

either_way ifthen job = ifthen job job
or_if_no_pml :: (Fancy String f, Fancy String g) => f -> g -> B2 String
or_if_no_pml x y = x `have_m_pml_or` (x `have_p_pml_or` y)
have_m_or = concise_if_else "have_m" have_m (\v b -> b { have_m = Just v })
have_p_or = concise_if_else "have_p" have_p (\v b -> b { have_p = Just v })
if_have_m = gen_if_else "have_m" have_m (\v b -> b { have_m = Just v })
if_have_p = gen_if_else "have_p" have_p (\v b -> b { have_p = Just v })
have_m_pml_or = concise_if_else "have_m_pml"
                (\b -> case have_m b of
                       Just False -> Just False
                       _ -> have_m_pml b)
                (\v b -> if v then b { have_m_pml = Just True,
                                       have_m = Just True }
                              else b { have_m_pml = Just False })
have_p_pml_or = concise_if_else "have_p_pml"
                (\b -> case have_p b of
                       Just False -> Just False
                       _ -> have_p_pml b)
                (\v b -> if v then b { have_p_pml = Just True,
                                       have_p = Just True }
                              else b { have_p_pml = Just False })
if_have_m_pml = gen_if_else "have_m_pml"
                (\b -> case have_m b of
                       Just False -> Just False
                       _ -> have_m_pml b)
                (\v b -> if v then b { have_m_pml = Just True,
                                       have_m = Just True }
                              else b { have_m_pml = Just False })
if_have_p_pml = gen_if_else "have_p_pml"
                (\b -> case have_p b of
                       Just False -> Just False
                       _ -> have_p_pml b)
                (\v b -> if v then b { have_p_pml = Just True,
                                       have_p = Just True }
                              else b { have_p_pml = Just False })

else_do a b = a b
infixr 1 `else_do`

gen_if_else :: (Fancy [String] a, Fancy [String] b) =>
               String -> (Bools -> Maybe Bool) -> (Bool -> Bools -> Bools)
            -> a -> b -> B2 [String]
gen_if_else name get put thendo elsedo = B2 gie where
  gie b | get b == Just True = thendo `given` b
        | get b == Just False = elsedo `given` b
        | otherwise =
          case thendo `given` (put True b) of
          [";"] -> case elsedo `given` (put False b) of
                   [";"] -> [] :: [String]
                   [] -> [] :: [String]
                   ed -> ["if (!"++name++") {"] ++
                           map ("  "++) ed ++
                         ["}"]
          [] -> case elsedo `given` (put False b) of
                [";"] -> [] :: [String]
                [] -> [] :: [String]
                ed -> ["if (!"++name++") {"] ++
                        map ("  "++) ed ++
                      ["}"]
          td -> case elsedo `given` (put False b) of
                [";"] -> ["if ("++name++") {"] ++
                           map ("  "++) td ++
                         ["}"]
                [] -> ["if ("++name++") {"] ++
                        map ("  "++) td ++
                      ["}"]
                ed | ed /= td -> ["if ("++name++") {"] ++
                                   map ("  "++) td ++
                                 ["} else {"] ++
                                   map ("  "++) ed ++
                                 ["}"]
                   | otherwise -> ed

concise_if_else :: (Fancy String a, Fancy String b) =>
                   String -> (Bools -> Maybe Bool) -> (Bool -> Bools -> Bools)
                -> a -> b -> B2 String
concise_if_else name get put thendo elsedo = B2 $ cie where
  cie b
      | get b == Just True = thendo `given` b
      | get b == Just False = elsedo `given` b
      | otherwise = "("++name++")?" ++ thendo `given` (put True b) ++
                                ":" ++ elsedo `given` (put False b)

cif bool inside = ["if ("++bool++") {"] ++ map ("  "++) inside ++ ["}"]

ifelse bool thendo elsedo =
    ["if ("++bool++") {"] ++ map ("  "++) thendo ++
    ["} else {"] ++ map ("  "++) elsedo ++ ["}"]

loop_over inside = B2 $ lo where
  lo b =
    case inside `given` b of
    [] -> []
    [";"] -> []
    ins -> ["for (int i1 = 0; i1 < n1; i1++)",
            " for (int i2 = 0; i2 < n2; i2++)",
            "  for (int i3 = 0, ind = i1*s1+i2*s2; i3 < n3; i3++, ind += s3) {"]
           ++ map ("   "++) ins ++ [" }"]

class (Fancy s) a where
    given :: a -> Bools -> s

instance (Fancy a) a where
    a `given` _ = a
instance (Fancy [a]) [B2 a] where
    xs `given` b = map (`given` b) xs
newtype (B2 a) = B2 (Bools -> a)
instance (Fancy a) (B2 a) where
    (B2 f) `given` b = f b
