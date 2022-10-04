TITL ''@ water.res
CELL  0.71073 16.193 16.193 11.2421 90 90 120
ZERR 6  0.00150  0.00150  0.00110  0.00000  0.00000  0.00000
LATT -1
SYMM  -y, +x-y, +z
SYMM  +y, +x, 1/2-z
SYMM  -x+y, -x, +z
SYMM  -x, -x+y, 1/2-z
SYMM  +x-y, -y, 1/2-z
SYMM  -x, -y, -z
SYMM  +y, -x+y, -z
SYMM  -y, -x, 1/2+z
SYMM  +x-y, +x, -z
SYMM  +x, +x-y, 1/2+z
SYMM  -x+y, +y, 1/2+z
SYMM  2/3+x, 1/3+y, 1/3+z
SYMM  1/3+x, 2/3+y, 2/3+z
SYMM  2/3-y, 1/3+x-y, 1/3+z
SYMM  1/3-y, 2/3+x-y, 2/3+z
SYMM  2/3+y, 1/3+x, 5/6-z
SYMM  1/3+y, 2/3+x, 1/6-z
SYMM  2/3-x+y, 1/3-x, 1/3+z
SYMM  1/3-x+y, 2/3-x, 2/3+z
SYMM  2/3-x, 1/3-x+y, 5/6-z
SYMM  1/3-x, 2/3-x+y, 1/6-z
SYMM  2/3+x-y, 1/3-y, 5/6-z
SYMM  1/3+x-y, 2/3-y, 1/6-z
SYMM  2/3-x, 1/3-y, 1/3-z
SYMM  1/3-x, 2/3-y, 2/3-z
SYMM  2/3+y, 1/3-x+y, 1/3-z
SYMM  1/3+y, 2/3-x+y, 2/3-z
SYMM  2/3-y, 1/3-x, 5/6+z
SYMM  1/3-y, 2/3-x, 1/6+z
SYMM  2/3+x-y, 1/3+x, 1/3-z
SYMM  1/3+x-y, 2/3+x, 2/3-z
SYMM  2/3+x, 1/3+x-y, 5/6+z
SYMM  1/3+x, 2/3+x-y, 1/6+z
SYMM  2/3-x+y, 1/3+y, 5/6+z
SYMM  1/3-x+y, 2/3+y, 1/6+z
SFAC Fe Cl O H 
UNIT 6.00012 0 36 0

REM ######################################################
REM This file exported by ShelXle is for information or 
REM visualiztion purposes only. You may run into trouble 
REM if you try to refine it against data.
REM DISABLE_REFINE
REM ######################################################

FVAR 1.00 0.31437 0.77327 

FE1   1   0.000000   0.000000   0.500000  10.16667   0.01569   0.01569 =
      0.02514   0.00000   0.00000   0.00785
O1    3   0.074199   0.116656   0.399075  11.00000   0.01652   0.01952 =
      0.03410   0.00449  -0.00042   0.00501
O11    3  -0.116656  -0.042457   0.399075  11.00000  0.01952   0.02602 = 
 0.03410  -0.00491  -0.00449   0.01451
rem  0.02153 0.04405 =
rem 0.0341 0.00449 -0.00042 0.02654

O12    3   0.042457  -0.074199   0.399075  11.00000   0.02602   0.01652 =
      0.03410   0.00042   0.00491   0.01151
O13    3  -0.074199  -0.116656   0.600925  11.00000   0.01652   0.01952 =
      0.03410   0.00449  -0.00042   0.00501
O14    3   0.116656   0.042457   0.600925  11.00000   0.01952   0.02602 =
      0.03410  -0.00491  -0.00449   0.01451
O15    3  -0.042457   0.074199   0.600925  11.00000   0.02602   0.01652 =
      0.03410   0.00042   0.00491   0.01151


EQIV $1  -y, +x-y, +z
EQIV $2  -x+y, -x, +z
EQIV $3  -x, -y, 1-z
EQIV $4  +y, -x+y, 1-z
EQIV $5  +x-y, +x, 1-z


BIND O1 FE1
BIND_$1 O1 FE1
BIND_$2 O1 FE1
BIND_$3 O1 FE1
BIND_$4 O1 FE1
BIND_$5 O1 FE1

HKLF 0 !don't refine this
END
