TITL post_end in P2(1)/c
REM Synthetic fixture for plan item B2a / D-10.
REM A bare SADI makes SHELXL emit the derived restraints *after* END,
REM using ACTA TABS part suffixes (^a = PART 1, ^b = PART 2) which
REM SHELXL itself cannot read back in. Those lines must be treated as
REM inert output: never atom-linked, never validated, never rewritten.
CELL 0.71073 10.5086 20.9035 20.5072 90 94.13 90
ZERR 4 0.0003 0.0005 0.0005 0 0.001 0
LATT 1
SYMM -X, 0.5+Y, 0.5-Z
SFAC C H O
UNIT 8 8 4
ACTA
SADI
SAME C1 > C4
WGHT 0.1
FVAR 1.0
PART 1
C1    1     0.10000  0.20000  0.30000  21.00000  0.05
C2    1     0.15000  0.25000  0.35000  21.00000  0.05
PART 2
C3    1     0.20000  0.30000  0.40000 -21.00000  0.05
C4    1     0.25000  0.35000  0.45000 -21.00000  0.05
PART 0
O1    3     0.30000  0.40000  0.50000  11.00000  0.05
HKLF 4
END

WGHT      0.0181     32.5755

REM SADI restraints, including those derived from SAME
REM This list may refer to PART numbers by a,b,c etc. but
REM this has not yet been implemented for input into SHELXL
SADI 0.0200  C1^a C2^a  C3^b C4^b
SADI 0.0400  C1^a C4^b  C2^a C3^b

REM Highest difference peak  0.356,  deepest hole -0.311,  1-sigma level  0.079
Q1    1   0.65660  0.40640  0.66150  11.00000  0.05    0.36
