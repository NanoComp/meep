module Main (main) where

import IO
import System
import FiniteMap
import RegexString
import List
import Monad

main = do
  argv <- getArgs;
  if length argv /= 2
   then do
       progname <- getProgName
       putStr $ "Usage: "++progname++" <file> <name>\n"
       exitWith $ ExitFailure 2
   else return () -- null operation
  (fname:bname:rest) <- getArgs;
  bls <- liftM (parse_band_file bname) (readFile fname)
  writeFile (fname++".xmgr") $ write_bands bls
  system $ "gracebat -hdevice EPS -printfile "++fname++".eps -hardcopy "++fname++".xmgr"
  exitWith ExitSuccess

newtype WaveVector = K Double
    deriving (Eq)
instance Show WaveVector where
    show (K k) = show k
data BandLine = BandLine WaveVector Int Int Double Double
                deriving (Eq,Show)

parse_band_file :: String -> String -> [BandLine]
get_freq :: [BandLine] -> Int -> Int -> WaveVector -> Double
get_ms :: [BandLine] -> [Int]
get_ns :: [BandLine] -> [Int]
get_ks :: Int -> Int -> [BandLine] -> [WaveVector]
plot_band :: Int -> Int -> [BandLine] -> String

write_bands :: [BandLine] -> String
write_bands bls = foldl (++) header [plot_band m n bls | m <- get_ms bls, n <- get_ns bls]

plot_band m n bls = case get_ks m n bls of
  [] -> ""
  ks ->
    color (1+m) setstart++
    "@target G0.S"++show setstart++"\n@type xy\n"++
    (foldl (++) "" $ map (\ k-> show k++" "++show (get_freq bls m n k)++"\n") ks)
    where setstart = m*1000+n
color c set = "@    s"++show set++" line color "++show c++"\n"++
              "@    s"++show set++" symbol color "++show c++"\n"

get_m (BandLine _ m _ _ _) = m
get_ms bls = nub $ map get_m bls

get_n (BandLine _ _ n _ _) = n
get_ns bls = nub $ map get_n bls

get_ks m n (BandLine k m' n' _ _ : bls) | n == n' && m == m' = nub $ k : get_ks m n bls
get_ks m n (_:bls) = get_ks m n bls
get_ks _ _ [] = []

get_freq (BandLine k' m' n' f _:bls) m n k
    | k == k' && m == m' && n == n' = f
get_freq (_:bls) m n k = get_freq bls m n k

parse_band_file name cs = pbf name $ lines cs
pbf _ [] = []
pbf name (l:ls) =
  case parse_band_line name l of
  Just bl -> bl : pbf name ls
  Nothing -> pbf name ls

parse_freq_line name l =
  case matchRegex 
           (mkRegex (name++nonspaces_items 5)) l of
  Just [sk,sm,sn,sf,sd] ->
      case (read sk, read sm, read sn, read sf, read sd) of
      (k,m,n,f,d) -> Just $ BandLine (K k) m n f d
  _ -> Nothing
parse_band_line :: String -> String -> Maybe BandLine
parse_band_line name l = parse_freq_line name l

nonspaces_items 0 = ""
nonspaces_items n = " ([^ ]+)" ++ nonspaces_items (n-1)

header = "# Grace project file
     #
     @page size 792, 612
     @page scroll 5%
     @page inout 5%
     @default symbol size 0.330000
"
