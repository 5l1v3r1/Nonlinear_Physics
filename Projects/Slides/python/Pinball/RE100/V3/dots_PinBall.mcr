#!MC 1200

# VARIABLES



$!VarSet |alfa1| = 0
$!VarSet |alfa2| = 0
$!VarSet |alfa3| = 0
$!VarSet  |T| = 0



$!LOOP 100 
$!VarSet |MFBD| = './' 
$!VarSet |FileBase|='|MFBD|TEC_U_|LOOP%03d|'
$!ReadDataSet '"|FileBase|.plt"'
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN2D
  VARNAMELIST = '"X" "Y"  "U" "V" "P" "VORT"'


$!LOOP |NUMZONES|


$!VarSet  |T| = (|T| + 0.1)



$!if |T| < 50.0

$!VarSet |C1| = 0
$!VarSet |C2| = 0
$!VarSet |C3| = 0

$!elseif  |T| < 100

$!VarSet |C1| = -2
$!VarSet |C2| =  2
$!VarSet |C3| =  0

#C1 = -2.0 ; C2 = 2.0 ; C3 = 0.0 ;

$!elseif |T| < 150

$!VarSet |C1| =  2
$!VarSet |C2| = -2
$!VarSet |C3| =  0

#C1 = 2.0 ; C2 = -2.0 ; C3 = 0.0 ;

$!elseif  |T| < 200

$!VarSet |C1| =  2
$!VarSet |C2| = -2
$!VarSet |C3| =  1

#C1 = 2.0 ; C2 = -2.0 ; C3 = 1.0 ;

$!elseif   |T| < 250 

$!VarSet |C1| =  2
$!VarSet |C2| = -2
$!VarSet |C3| = -1

#C1 = 2.0 ; C2 = -2.0 ; C3 = -1.0 ;

$!else

$!VarSet |C1| =  0
$!VarSet |C2| =  0
$!VarSet |C3| =  0

#C1 = 0.0 ; C2 = 0.0 ; C3 = 0.0 ;

$!endif



$!VarSet |alfa1| = (|alfa1| + |C1| * 0.2  )
$!VarSet |alfa2| = (|alfa2| + |C2| * 0.2  ) 
$!VarSet |alfa3| = (|alfa3| + |C3| * 0.2  )  



$!VarSet |X1| = (0.4*cos(|alfa1|) +  0.00)
$!VarSet |Y1| = (0.4*sin(|alfa1|) +  0.75)
$!ATTACHGEOM 
  GEOMTYPE = CIRCLE
  ANCHORPOS
    {
    X = |X1|
    Y = |Y1|
    }
  RAWDATA
0.05

$!VarSet |X1| = (0.4*cos(|alfa2|) +  0.00)
$!VarSet |Y1| = (0.4*sin(|alfa2|) -  0.75)
$!ATTACHGEOM 
  GEOMTYPE = CIRCLE
  ANCHORPOS
    {
    X = |X1|
    Y = |Y1|
    }
  RAWDATA
0.05

$!VarSet |X1| = (0.4*cos(|alfa3|)    -1.5 * sqrt(3./4.) )
$!VarSet |Y1| = (0.4*sin(|alfa3|) +                0.0 )
$!ATTACHGEOM 
  GEOMTYPE = CIRCLE
  ANCHORPOS
    {
    X = |X1|
    Y = |Y1|
    }
  RAWDATA
0.05


$!ATTACHTEXT
XYPOS
{
X = 73.0
Y = 5.5
}
TEXT = "TIME = |ZONENAME[|loop|]|"
TEXTSHAPE
{
HEIGHT = 12
}



$!ACTIVEFIELDMAPS = [|LOOP|]
#####$!ACTIVEFIELDMAPS += [|NUMZONES|]


$!TWODAXIS PRESERVEAXISSCALE = YES
$!TWODAXIS YDETAIL{RANGEMAX = 6}
$!TWODAXIS YDETAIL{RANGEMIN = -6}
$!TWODAXIS XDETAIL{RANGEMAX = 20}
$!TWODAXIS XDETAIL{RANGEMIN = -6}
$!REDRAWALL
$!GLOBALCONTOUR 1  VAR = 5
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 1
  APPROXNUMVALUES = 15
$!GLOBALCONTOUR 1  VAR = 6
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 1
  APPROXNUMVALUES = 15
$!GLOBALCONTOUR 1  COLORMAPFILTER{COLORMAPDISTRIBUTION = CONTINUOUS}
$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN = -10}}
$!GLOBALCONTOUR 1  COLORMAPFILTER{CONTINUOUSCOLOR{CMAX = 10}}
$!FIELDLAYERS SHOWCONTOUR = YES
$!REDRAWALL

$!VIEW PUSH
$!REDRAWALL 

$!EXPORTSETUP EXPORTFORMAT = PNG
$!EXPORTSETUP IMAGEWIDTH = 1400
$!EXPORTSETUP EXPORTFNAME = '|FileBase|_|LOOP%04d|.png'
$!EXPORT
  EXPORTREGION = CURRENTFRAME

$!ENDLOOP


$!EndLoop





