module Main (main) where

import IO
import System
import FiniteMap
import RegexString
import Monad
import List

modename fname kname mname = fname++"-"++kname++"-"++mname

main = do
  argv <- getArgs;
  if length argv /= 4
   then do
       progname <- getProgName
       putStr $ "Usage: "++progname++" <file> <name> <k> <m>\n"
       exitWith $ ExitFailure 2
   else return () -- null operation
  [fname,bname,kname,mname] <- getArgs;
  bls <- liftM (parse_band_file bname) $ readFile fname
  case (K $ read kname, read mname :: Int, modename fname kname mname) of
       (k,m,out) -> do
       writeFile (out++".xmgr") $
                 (foldl (++) header [plot_mode k m n bls | n <- get_ns k m bls])
       system $ "gracebat -printfile "++out++".ps -hardcopy "++out++".xmgr"
       exitWith ExitSuccess

newtype Radius = R Double
    deriving (Eq,Show,Ord)
newtype WaveVector = K Double
    deriving (Eq,Show)
data Field = Field Double Double Double Double Double Double
             deriving (Eq,Show) 
data BandLine = BandLine WaveVector Int Int Double Double
              | FieldLine WaveVector Int Int Radius Field
                deriving (Eq,Show)

parse_band_file :: String -> String -> [BandLine]
get_freq :: WaveVector -> Int -> Int -> [BandLine] -> Double
get_mode :: WaveVector -> Int -> Int -> [BandLine] -> [(Radius, Field)]
get_ms :: WaveVector -> [BandLine] -> [Int]
get_ns :: WaveVector -> Int -> [BandLine] -> [Int]
plot_mode :: WaveVector -> Int -> Int -> [BandLine] -> String

get_dr :: [(Radius, Field)] -> Double
get_dr mo = case (fst $mo!!1) of R r -> r
get_rmax :: [(Radius, Field)] -> Double
get_rmax mo = case head $ reverse $ sort (map fst mo) of R r -> r

er (R r, Field er ep ez hr hp hz) = (r, er)
ep (R r, Field er ep ez hr hp hz) = (r, ep)
ez (R r, Field er ep ez hr hp hz) = (r, ez)
hr (R r, Field er ep ez hr hp hz) = (r, hr)
hp (R r, Field er ep ez hr hp hz) = (r, hp)
hz (R r, Field er ep ez hr hp hz) = (r, hz)
shiftr dr (r,f) = (r+dr*0.5,f)
power (Field er ep ez hr hp hz) = er*er+ep*ep+ez*ez+hr*hr+hp*hp+hz*hz

plot_mode k m n bls =
    graphheader n++
    red    (0+setstart)++
    star   (0+setstart)++
    red    (1+setstart)++
    circle (1+setstart)++
    blue   (2+setstart)++
    square (2+setstart)++
    blue   (3+setstart)++
    star   (3+setstart)++
    dashed (3+setstart)++
    blue   (4+setstart)++
    circle (4+setstart)++
    dashed (4+setstart)++
    red    (5+setstart)++
    square (5+setstart)++
    dashed (5+setstart)++
    "@target G"++show n++".S"++show (0+setstart)++"\n@type xy\n"++
               plot_field (map (shiftr dr.er) fs)++
    "@target G"++show n++".S"++show (1+setstart)++"\n@type xy\n"++
               plot_field (map ep fs)++
    "@target G"++show n++".S"++show (2+setstart)++"\n@type xy\n"++
               plot_field (map ez fs)++
    "@target G"++show n++".S"++show (3+setstart)++"\n@type xy\n"++
               plot_field (map hr fs)++
    "@target G"++show n++".S"++show (4+setstart)++"\n@type xy\n"++
               plot_field (map (shiftr dr.hp) fs)++
    "@target G"++show n++".S"++show (5+setstart)++"\n@type xy\n"++
               plot_field (map (shiftr dr.hz) fs)
    where setstart = m*1000
          fs = norm_mode $ get_mode k m n bls
          rmax = get_rmax fs
          dr = get_dr fs
plot_field ((r,f):ps) = show r++" "++show f++"\n"++plot_field ps
plot_field [] = "\n"
circle set = "@    s"++show set++" symbol 1\n"
star   set = "@    s"++show set++" symbol 9\n"
square set = "@    s"++show set++" symbol 2\n"
red  set = "@    s"++show set++" line color 2
@    s"++show set++" symbol color 2\n"
blue set = "@    s"++show set++" line color 4
@    s"++show set++" symbol color 4\n"
dashed set = "@    s"++show set++" line linestyle 3\n"
len [] = 0.0
len (p:ps) = 1.0 + len ps
norm_mode ps = 
    case sqrt $ (sum $ map (power . snd) ps) / len ps of
    n -> map (\ (r,Field er ep ez hr hp hz) ->
              (r,Field (er/n) (ep/n) (ez/n) (hr/n) (hp/n) (hz/n))) ps

get_ms k (BandLine k' m _ _ _ : bls)
    | k == k' = m : get_ms k bls
get_ms k (_:bls) = get_ms k bls
get_ms _ [] = []

get_ns k m (BandLine k' m' n _ _ : bls)
    | k == k' && m == m' = n : get_ns k m bls
get_ns k m (_:bls) = get_ns k m bls
get_ns _ _ [] = []

get_mode k m n (FieldLine k' m' n' r f:bls)
    | k == k' && m == m' && n == n' = (r, f) : get_mode k m n bls
get_mode k m n (_:bls) = get_mode k m n bls
get_mode _ _ _ [] = []

get_freq k m n (BandLine k' m' n' f _:bls)
    | k == k' && m == m' && n == n' = f
get_freq k m n (_:bls) = get_freq k m n bls

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
parse_field_line name l =
  case matchRegex
           (mkRegex (name++"-fields"++nonspaces_items 16)) l of
  Just [sk,sm,sn,sr,ser,ier,sep,iep,sez,iez,shr,ihr,shp,ihp,shz,ihz] ->
      case (read sk,read sm,read sn,read sr,
            read ser,read sep,read sez,
            read shr,read shp,read shz) of
      (k,m,n,r,er,ep,ez,hr,hp,hz) ->
          Just $ FieldLine (K k) m n (R r) (Field er ep ez hr hp hz)
  _ -> Nothing
parse_band_line :: String -> String -> Maybe BandLine
parse_band_line name l =
    case parse_freq_line name l of
    Just x -> Just x
    Nothing -> parse_field_line name l

nonspaces_items 0 = ""
nonspaces_items n = " ([^ ]+)" ++ nonspaces_items (n-1)

header = "# Grace project file
     #
     @page size 792, 612
     @page scroll 5%
     @page inout 5%
     @default symbol size 0.330000
"

graphheader intn = 
    let n = show (intn :: Int)
        numdown = (fromInteger $ floor (fromIntegral intn/2.0)) :: Double
        xmin = show $ 0.15+(fromIntegral intn/2.0 - numdown)
        xmax = show $ 0.65+(fromIntegral intn/2.0 - numdown)
        ymin = show $ 0.75-numdown * 0.2
        ymax = show $ 0.95-numdown * 0.2
    in "\n@g"++n++" on
     @g"++n++" hidden false
     @g"++n++" type XY
     @g"++n++" stacked false
     @g"++n++" bar hgap 0.000000
     @g"++n++" fixedpoint off
     @g"++n++" fixedpoint type 0
     @g"++n++" fixedpoint xy 0.000000, 0.000000
     @g"++n++" fixedpoint format general general
     @g"++n++" fixedpoint prec 6, 6
     @with g"++n++"
     @    znorm 1
     @    world xmin 0.0
     @    world xmax 1.0
     @    world ymin -2.0
     @    world ymax 2.0
     @    view xmin "++xmin++"
     @    view xmax "++xmax++"
     @    view ymin "++ymin++"
     @    view ymax "++ymax++"
     @    autoscale onread none
          "
