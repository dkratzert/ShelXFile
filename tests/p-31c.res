REM Solution 1  R1  0.081,  Alpha = 0.0146  in P31c
REM Flack x = -0.072 ( 0.041 ) from Parsons' quotients
REM C19.667 N2.667 P4
TITL p-31c-neu in P-31c
CELL  0.71073  12.5067  12.5067  24.5615   90.000   90.000  120.000
ZERR   2   0.0043   0.0043   0.0085   0.000   0.000   0.000
LATT -1
SYMM -Y, X-Y, Z
SYMM -X+Y, -X, Z
SYMM Y, X, 1/2+Z
SYMM X-Y, -Y, 1/2+Z
SYMM -X, -X+Y, 1/2+Z
SFAC C H N P Cl
UNIT 120 186 14 12 12
TEMP -173.000
OMIT 0   0   2
L.S. 10
BOND $H
ACTA
CONF
EQIV $1 -y+1, x-y, z
HTAB N1 Cl1_$1
EQIV $2 -x+y, -x+1, z
HTAB C3 Cl1_$2
HTAB N1' Cl1_$1
HTAB N2 Cl2
EQIV $3 -x+y+1, -x+1, z
HTAB C14 Cl2_$3
HTAB N2' Cl2
HTAB C14' Cl2_$3
LIST 4
FMAP 2
PLAN 40
WGHT    0.034600    0.643600
FVAR       0.22604   0.76052   0.85152
CL1   5    0.629304    0.639920    0.624939    11.00000    0.01547    0.01791 =
         0.02428   -0.00077   -0.00018    0.00586
CL2   5    0.304944    0.322118    0.301971    11.00000    0.01480    0.01856 =
         0.02708    0.00328    0.00049    0.00870
N3    3    0.333333    0.666667    0.451910    10.33333    0.03861    0.03861 =
         0.08508    0.00000    0.00000    0.01931
C23   1    0.333333    0.666667    0.557385    10.33333    0.03402    0.03402 =
         0.08915    0.00000    0.00000    0.01701
AFIX 137
H23A  2    0.402538    0.658074    0.570685    10.33333   -1.50000
H23B  2    0.255536    0.597462    0.570685    10.33333   -1.50000
H23C  2    0.341926    0.744464    0.570685    10.33333   -1.50000
AFIX   0
C24   1    0.333333    0.666667    0.498258    10.33333    0.02354    0.02354 =
         0.08455    0.00000    0.00000    0.01177
REM  ##########
REM  Molekuel 1
REM  ##########
 
P1    4    0.342244    0.295314    0.631261    11.00000    0.01034    0.01241 =
         0.01613   -0.00002    0.00018    0.00416
SADI N1 P1 N1' P1
SIMU P1 > C3'
RIGU P1 > C3'
DELU P1 > C3'
EADP C2 C2'
EADP N1 N1'
EADP C3 C3'
DFIX 0.91 N1 H1 N1' H1'
FLAT 0.1 P1 N1 C3 H1
FLAT 0.1 P1 N1' C3' H1'
SADI H1 P1 H1' P1
SADI H1 N1 H1' N1'
PART 1
N1    3    0.226172    0.154594    0.640737    31.00000    0.01038    0.01176 =
         0.02529    0.00061   -0.00164    0.00435
H1    2    0.245554    0.097197    0.643752    31.00000   -1.30000
C1    1    0.000000    0.000000    0.700589    30.33333    0.01374    0.01374 =
         0.01788    0.00000    0.00000    0.00687
AFIX 137
H1A   2    0.035397   -0.049518    0.713889    30.33333   -1.50000
H1B   2   -0.084915   -0.035397    0.713889    30.33333   -1.50000
H1C   2    0.049518    0.084915    0.713889    30.33333   -1.50000
AFIX   0
C2    1    0.000000    0.000000    0.638319    30.33333    0.00953    0.00953 =
         0.01652    0.00000    0.00000    0.00476
C3    1    0.106350    0.120870    0.615617    31.00000    0.00828    0.01237 =
         0.02284    0.00192   -0.00063    0.00365
AFIX  23
H3A   2    0.112067    0.112325    0.575810    31.00000   -1.20000
H3B   2    0.088009    0.188258    0.621903    31.00000   -1.20000
AFIX   0
PART 2
SAME N1 > C3
N1'   3    0.227651    0.152870    0.641741   -31.00000    0.01038    0.01176 =
         0.02529    0.00061   -0.00164    0.00435
H1'   2    0.244709    0.096442    0.628554   -31.00000   -1.30000
C1'   1    0.000000    0.000000    0.583017   -30.33333    0.03377    0.03377 =
         0.01385    0.00000    0.00000    0.01689
AFIX 137
H1D   2    0.011195   -0.067640    0.569717   -30.33333   -1.50000
H1E   2    0.067640    0.078835    0.569717   -30.33333   -1.50000
H1F   2   -0.078835   -0.011195    0.569717   -30.33333   -1.50000
AFIX   0
C2'   1    0.000000    0.000000    0.645158   -30.33333    0.00953    0.00953 =
         0.01652    0.00000    0.00000    0.00476
C3'   1    0.107648    0.120397    0.666075   -31.00000    0.00828    0.01237 =
         0.02284    0.00192   -0.00063    0.00365
AFIX  23
H3C   2    0.090221    0.188045    0.658750   -31.00000   -1.20000
H3D   2    0.113405    0.114220    0.706037   -31.00000   -1.20000
AFIX   0
PART 0
 
C4    1    0.306898    0.405021    0.659149    11.00000    0.01911    0.01442 =
         0.02384   -0.00287    0.00182    0.00824
AFIX 137
H4A   2    0.233993    0.397954    0.640859    11.00000   -1.50000
H4B   2    0.377153    0.488398    0.653805    11.00000   -1.50000
H4C   2    0.290042    0.389252    0.698171    11.00000   -1.50000
AFIX   0
C5    1    0.470340    0.305628    0.667286    11.00000    0.01645    0.02213 =
         0.02460   -0.00207   -0.00391    0.00888
AFIX 137
H5A   2    0.453211    0.298639    0.706457    11.00000   -1.50000
H5B   2    0.544083    0.385295    0.659564    11.00000   -1.50000
H5C   2    0.484349    0.238639    0.655701    11.00000   -1.50000
AFIX   0
C6    1    0.377381    0.323976    0.560310    11.00000    0.01621    0.01254 =
         0.01794   -0.00086    0.00057    0.00394
C7    1    0.460211    0.294079    0.536081    11.00000    0.01889    0.02429 =
         0.02160    0.00001    0.00120    0.01102
AFIX  43
H7    2    0.506217    0.268862    0.558026    11.00000   -1.20000
AFIX   0
C8    1    0.475279    0.301270    0.479751    11.00000    0.02231    0.02599 =
         0.02566   -0.00029    0.00823    0.01145
AFIX  43
H8    2    0.532817    0.282415    0.463497    11.00000   -1.20000
AFIX   0
C9    1    0.407502    0.335454    0.447486    11.00000    0.02848    0.02243 =
         0.01793    0.00117    0.00514    0.00634
AFIX  43
H9    2    0.417393    0.338856    0.409062    11.00000   -1.20000
AFIX   0
C10   1    0.324163    0.365164    0.471163    11.00000    0.02491    0.02499 =
         0.02371    0.00157   -0.00258    0.01147
AFIX  43
H10   2    0.277551    0.389036    0.448883    11.00000   -1.20000
AFIX   0
C11   1    0.309483    0.359754    0.527451    11.00000    0.02372    0.02352 =
         0.02150   -0.00014    0.00054    0.01313
AFIX  43
H11   2    0.253148    0.380441    0.543576    11.00000   -1.20000
AFIX   0
REM  ##########
REM  Molekuel 2
REM  ##########
P2    4    0.637854    0.639076    0.312737    11.00000    0.01285    0.01096 =
         0.01797    0.00072   -0.00012    0.00658
SADI N2 P2 N2' P2
SIMU P2 > C14'
RIGU P2 > C14'
DELU P2 > C14'
EADP C13 C13'
EADP N2 N2'
DFIX 0.91 N2 H2 N2' H2'
FLAT 0.1 P2 N2 C14 H2
FLAT 0.1 P2 N2' C14' H2'
SADI H2 P2 H2' P2
SADI H2 N2 H2' N2'
PART 1
N2    3    0.600731    0.494473    0.306659    21.00000    0.01169    0.01129 =
         0.02580   -0.00020   -0.00113    0.00644
H2    2    0.523194    0.442256    0.300634    21.00000   -1.30000
C12   1    0.666667    0.333333    0.242844    20.33333    0.01892    0.01892 =
         0.01486    0.00000    0.00000    0.00946
AFIX 137
H12A  2    0.601767    0.348829    0.229544    20.33333   -1.50000
H12B  2    0.651171    0.252938    0.229544    20.33333   -1.50000
H12C  2    0.747062    0.398233    0.229544    20.33333   -1.50000
AFIX   0
C13   1    0.666667    0.333333    0.305063    20.33333    0.00903    0.00903 =
         0.01484    0.00000    0.00000    0.00451
C14   1    0.687009    0.456818    0.328297    21.00000    0.01317    0.01186 =
         0.02049   -0.00179   -0.00224    0.00769
AFIX  23
H14A  2    0.678268    0.449665    0.368394    21.00000   -1.20000
H14B  2    0.772526    0.522340    0.320039    21.00000   -1.20000
AFIX   0
PART 2
SAME N2 > C14
N2'   3    0.600564    0.495293    0.299217   -21.00000    0.01169    0.01129 =
         0.02580   -0.00020   -0.00113    0.00644
H2'   2    0.525333    0.440958    0.311051   -21.00000   -1.30000
C12'  1    0.666667    0.333333    0.361432   -20.33333    0.01781    0.01781 =
         0.01786    0.00000    0.00000    0.00891
AFIX 137
H12D  2    0.590714    0.328990    0.374732   -20.33333   -1.50000
H12E  2    0.738276    0.409286    0.374732   -20.33333   -1.50000
H12F  2    0.671009    0.261724    0.374732   -20.33333   -1.50000
AFIX   0
C13'  1    0.666667    0.333333    0.299495   -20.33333    0.00903    0.00903 =
         0.01484    0.00000    0.00000    0.00451
C14'  1    0.685815    0.456923    0.276937   -21.00000    0.01771    0.01365 =
         0.01618    0.00061    0.00330    0.00999
AFIX  23
H14C  2    0.771588    0.522362    0.284820   -21.00000   -1.20000
H14D  2    0.676051    0.450042    0.236874   -21.00000   -1.20000
AFIX   0
PART 0
 
C15   1    0.502784    0.648423    0.299794    11.00000    0.01803    0.01661 =
         0.03458    0.00038   -0.00408    0.01187
AFIX 137
H15A  2    0.472592    0.617989    0.263001    11.00000   -1.50000
H15B  2    0.522471    0.734514    0.302834    11.00000   -1.50000
H15C  2    0.438894    0.597914    0.326419    11.00000   -1.50000
AFIX   0
C16   1    0.755526    0.735949    0.266553    11.00000    0.02210    0.01480 =
         0.02134    0.00372    0.00485    0.00764
AFIX 137
H16A  2    0.831256    0.734138    0.274835    11.00000   -1.50000
H16B  2    0.771254    0.820682    0.269993    11.00000   -1.50000
H16C  2    0.729172    0.706624    0.229279    11.00000   -1.50000
AFIX   0
C17   1    0.688940    0.691498    0.381021    11.00000    0.02023    0.01524 =
         0.01557   -0.00067   -0.00127    0.00951
C18   1    0.811740    0.775620    0.393433    11.00000    0.01702    0.01825 =
         0.02742   -0.00005   -0.00088    0.00680
AFIX  43
H18   2    0.871214    0.811832    0.365196    11.00000   -1.20000
AFIX   0
C19   1    0.846488    0.806141    0.447733    11.00000    0.02643    0.01958 =
         0.03276   -0.00581   -0.00928    0.01050
AFIX  43
H19   2    0.929913    0.863860    0.456461    11.00000   -1.20000
AFIX   0
C20   1    0.760225    0.752877    0.488907    11.00000    0.03959    0.02466 =
         0.02435   -0.00554   -0.00733    0.02133
AFIX  43
H20   2    0.784691    0.774317    0.525753    11.00000   -1.20000
AFIX   0
C21   1    0.638202    0.668342    0.476723    11.00000    0.03528    0.03074 =
         0.02123    0.00234    0.00582    0.01795
AFIX  43
H21   2    0.579362    0.631004    0.505122    11.00000   -1.20000
AFIX   0
C22   1    0.602761    0.638784    0.423091    11.00000    0.02067    0.02078 =
         0.02371    0.00036    0.00306    0.00716
AFIX  43
H22   2    0.518887    0.581997    0.414676    11.00000   -1.20000
AFIX   0
HKLF 4
 
REM  p-31c-neu in P-31c
REM R1 =  0.0308 for    4999 Fo > 4sig(Fo)  and  0.0343 for all    5352 data
REM    287 parameters refined using    365 restraints
 
END  
     
WGHT      0.0348      0.6278 

REM Highest difference peak  0.224,  deepest hole -0.252,  1-sigma level  0.053
Q1    1   0.4074  0.2995  0.6529  11.00000  0.05    0.22
Q2    1   0.6988  0.6930  0.2844  11.00000  0.05    0.22
Q3    1   0.6041  0.5504  0.3254  11.00000  0.05    0.22
Q4    1   0.6769  0.6691  0.3478  11.00000  0.05    0.22
Q5    1   0.3605  0.3175  0.5940  11.00000  0.05    0.20
Q6    1   0.8270  0.7612  0.4237  11.00000  0.05    0.20
Q7    1   0.7509  0.7471  0.3914  11.00000  0.05    0.20
Q8    1   0.3555  0.3621  0.5501  11.00000  0.05    0.19
Q9    1   0.4416  0.3349  0.5451  11.00000  0.05    0.18
Q10   1   0.8179  0.8151  0.4169  11.00000  0.05    0.18
Q11   1   0.6601  0.6835  0.6536  11.00000  0.05    0.18
Q12   1   0.6445  0.6478  0.3989  11.00000  0.05    0.17
Q13   1   0.4557  0.2818  0.5083  11.00000  0.05    0.17
Q14   1   0.6550  0.6816  0.5952  11.00000  0.05    0.17
Q15   1   0.5339  0.6537  0.2284  11.00000  0.05    0.16
Q16   1   0.4164  0.3019  0.5464  11.00000  0.05    0.16
Q17   1   0.5877  0.4931  0.2723  11.00000  0.05    0.16
Q18   1   0.3259  0.3292  0.2650  11.00000  0.05    0.15
Q19   1   0.8494  0.8615  0.3688  11.00000  0.05    0.15
Q20   1   0.3351  0.3711  0.3327  11.00000  0.05    0.15
Q21   1   0.2410  0.1323  0.7123  11.00000  0.05    0.15
Q22   1   0.4796  0.4099  0.7164  11.00000  0.05    0.15
Q23   1   0.2910  0.1848  0.6176  11.00000  0.05    0.15
Q24   1   0.6196  0.4591  0.3825  11.00000  0.05    0.15
Q25   1   0.3189  0.4839  0.3973  11.00000  0.05    0.15
Q26   1   0.8096  0.6920  0.3679  11.00000  0.05    0.15
Q27   1   0.5753  0.6616  0.3039  11.00000  0.05    0.14
Q28   1   0.9225  0.8080  0.2662  11.00000  0.05    0.14
Q29   1   0.4167  0.2634  0.6893  11.00000  0.05    0.14
Q30   1   0.4788  0.2174  0.5458  11.00000  0.05    0.14
Q31   1   0.9991  0.8420  0.3674  11.00000  0.05    0.14
Q32   1   0.8971  0.8148  0.4136  11.00000  0.05    0.14
Q33   1   0.5976  0.4382  0.2218  11.00000  0.05    0.14
Q34   1   0.3133  0.3590  0.6334  11.00000  0.05    0.14
Q35   1   0.2423  0.7411  0.5639  11.00000  0.05    0.14
Q36   1   0.5629  0.6170  0.6068  11.00000  0.05    0.13
Q37   1   0.4419  0.4570  0.7147  11.00000  0.05    0.13
Q38   1   0.5666  0.4955  0.2257  11.00000  0.05    0.13
Q39   1   0.9806  0.7850  0.3765  11.00000  0.05    0.13
Q40   1   0.8809  0.9280  0.3680  11.00000  0.05    0.13
