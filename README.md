# ShelXFile

This is a full implementation of the SHELXL<sup>[[1](https://github.com/dkratzert/ShelXFile/blob/master/README.md#references)]</sup> file syntax. Additionally it is able to edit SHELX properties using Python.
The implementation is Python3-only and supports SHELXL after 2017 (You should not use old versions anyway).
ShelXFile may eventually become a new heart of DSR<sup>[[2](https://github.com/dkratzert/ShelXFile/blob/master/README.md#references)]</sup> and is already used as file parser in StructureFinder<sup>[[3](https://github.com/dkratzert/ShelXFile/blob/master/README.md#references)]</sup>.

ShelXFile always keeps the file order intact. Every SHELX instruction like DFIX or an atom is stored as an class object in the list ShelXlFile.\_reslist. When writing the ShelXlFile content to disk, it wites the \_reslist content to disk.

ShelXFile tries to detect all possible syntax errors that SHELXL would not like either. If ShelXFile.DEBUG is True, more output about syntax and other errors are printed out. Otherwise, the parser is quiet except for really severe errors like a missing unit cell.

Not every part of ShelXFile is complete, for example it will not recognize if you add restraints with atom names that are not in the SHELX file. Please help me improving it!

Examples:
```python

from shelxfile.shelx import ShelXlFile
shx = ShelXlFile('./tests/p21c.res')

shx.cell
[10.5086, 20.9035, 20.5072, 90.0, 94.13, 90.0]

shx.atoms
O1    3    0.074835    0.238436    0.402457   -31.00000    0.01579    0.03095    0.01852   -0.00468   -0.00210    0.01153
C1    1    0.028576    0.234542    0.337234   -31.00000    0.02311    0.03617    0.01096   -0.01000    0.00201    0.00356
C2    1    0.121540    0.194460    0.298291   -31.00000    0.02960    0.04586    0.01555   -0.00485   -0.00023    0.01102
...

a = shx.atoms.get_atom_by_name('F1_2')  # Atom F1 in residue 2

a
Atom ID: 273

str(a)
'F1    4    0.245205    0.192674    0.649231   -21.00000    0.05143    0.03826    0.03193   -0.00579   -0.01865   -0.00485'

a.to_isotropic()
str(a)
'F1    4    0.245205    0.192674   0.649231  -21.00000    0.04000'

a.position
273

str(shx._reslist[273])  # In regular code, do not access shx._reslist directly!
'F1    4    0.245205    0.192674   0.649231  -21.00000    0.04000'

a.name
'F1'

a.resinum
2

a.part
2

a.find_atoms_around(dist=2.0, only_part=1)
[Atom ID: 254, Atom ID: 256, Atom ID: 260]

[str(x) for x in a.find_atoms_around(dist=2.2, only_part=2)]
['C2   1    0.192984    0.140449    0.621265   -21.00000    0.04315    0.02747    0.02385    0.00686   -0.00757    0.00126', 
'F2    4    0.264027    0.090306    0.642441   -21.00000    0.06073    0.04450    0.03972    0.01630   -0.01260    0.01460', 
'F3    4    0.078582    0.131920    0.643529   -21.00000    0.05691    0.04955    0.03374    0.01040    0.01881    0.00375']

a.cart_coords
[1.617897551082389, 4.027560959000001, 13.279336538026433]

a.frac_coords
[0.245205, 0.192674, 0.649231]

a.occupancy
1.0

a.sfac_num
4

[x for x in a.find_atoms_around(dist=2.5, only_part=2)]
[Atom ID: 269, Atom ID: 271, Atom ID: 275, Atom ID: 277]

for x in a.find_atoms_around(dist=2.5, only_part=2):
    x.delete()
    
[x for x in a.find_atoms_around(dist=2.5, only_part=2)]
[]  # Atoms are now deleted from shx._reslist.

shx.restraints
SADI_CCF3 0.02 C1 C2 C1 C3 C1 C4
SADI_CCF3 0.02 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4
SADI_CCF3 0.04 C2 C3 C3 C4 C2 C4
SADI_CCF3 0.04 O1 C2 O1 C3 O1 C4
...

shx.restraints[1]
SADI_CCF3 0.02 C1 C2 C1 C3 C1 C4

shx.restraints[1].textline
'SADI_CCF3 0.02 C1 C2 C1 C3 C1 C4'

shx.restraints[1].residue_class
'CCF3'

shx.restraints[1].s
0.02
```
Writes current shx object to test.ins
All lines in ShelXlFile._reslist get wrapped after 79 characters with " =\n " as
specified by SHELXL during the file writing.
```python
shx.write_shelx_file('test.ins')
```
No matter if you loaded a .res or .ins file with refine(), SHELXL refines the structure of the ShelXlFile() object. 
The default for refine() are zero least squares cycles:
```python
shx.refine(2)
```
```
-------------------------------------------------------------------------------
 Running SHELXL with "/usr/local/bin/shelxl -b3000 /Users/daniel/GitHub/ShelXFile/tests/p21c" and "L.S. 2"
 wR2 =  0.1005 before cycle   1 for   10786 data and    945 /    945 parameters
 wR2 =  0.1005 before cycle   2 for   10786 data and    945 /    945 parameters
 wR2 =  0.1005 before cycle   3 for   10786 data and      0 /    945 parameters
 SHELXL Version 2018/3
```
## References
[1] http://shelx.uni-goettingen.de/, G. M. Sheldrick, Acta Cryst. (2015). C71, 3-8.
https://doi.org/10.1107/S2053229614024218

[2] https://github.com/dkratzert/DSR

[3] https://github.com/dkratzert/StructureFinder

