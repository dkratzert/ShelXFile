# ShelXFile

This is a full implementation of the SHELXL[1] file syntax. Additionally it is able to edit SHELX properties using Python.
The implementation is Python3-only and supports SHELXL after 2017 (You should not use old versions anyway).
ShelXFile may eventually become a new heart of DSR[2] and is already used as file parser in StructureFinder[3].

ShelXFile always keeps the file order intact. Every SHELX instruction like DFIX or an atom is stored as an class object in the list ShelXlFile.\_reslist. When writing the ShelXlFile content to disk, it wites the \_reslist content to disk.

Examples:
```python

from shelxfile.shelx import ShelXlFile
shx = ShelXlFile('./p21c.res')

shx.cell
[10.5086, 20.9035, 20.5072, 90.0, 94.13, 90.0]

shx.atoms
O1    3    0.074835    0.238436    0.402457   -31.00000    0.01579    0.03095      0.01852   -0.00468   -0.00210    0.01153
C1    1    0.028576    0.234542    0.337234   -31.00000    0.02311    0.03617      0.01096   -0.01000    0.00201    0.00356
C2    1    0.121540    0.194460    0.298291   -31.00000    0.02960    0.04586      0.01555   -0.00485   -0.00023    0.01102
...

a = shx.atoms.get_atom_by_name('F1_2')

a
ID: 255

str(a)
'F1    4    0.245205    0.192674    0.649231   -21.00000    0.05143    0.03826      0.03193   -0.00579   -0.01865   -0.00485'

a.to_isotropic()
str(a)
'F1    4    0.245205    0.192674   0.649231  -21.00000    0.04000'

a.name
'F1'

a.resinum
2

a.part
2

a.find_atoms_around(dist=2.0, only_part=1)
[ID: 236, ID: 238, ID: 242]

[str(x) for x in a.find_atoms_around(dist=2.2, only_part=2)]
['C2    1    0.192984    0.140449    0.621265   -21.00000    0.04315    0.02747      0.02385    0.00686   -0.00757    0.00126', 
'F2    4    0.264027    0.090306    0.642441   -21.00000    0.06073    0.04450      0.03972    0.01630   -0.01260    0.01460', 
'F3    4    0.078582    0.131920    0.643529   -21.00000    0.05691    0.04955      0.03374    0.01040    0.01881    0.00375']

a.cart_coords
[1.617897551082389, 4.027560959000001, 13.279336538026433]

a.frac_coords
[0.245205, 0.192674, 0.649231]

a.line_numbers
[255, 256]

a.occupancy
1.0

a.sfac_num
4

shx.write_shelx_file('test.ins')
(writes current shx object to test.ins)

shx.restraints
SADI_CCF3 0.02 C1 C2 C1 C3 C1 C4
SADI_CCF3 0.02 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4
SADI_CCF3 0.04 C2 C3 C3 C4 C2 C4
SADI_CCF3 0.04 O1 C2 O1 C3 O1 C4
...
```

[1] http://shelx.uni-goettingen.de/, G. M. Sheldrick, Acta Cryst. (2015). C71, 3-8. https://doi.org/10.1107/S2053229614024218

[2] https://github.com/dkratzert/DSR

[3] https://github.com/dkratzert/StructureFinder

