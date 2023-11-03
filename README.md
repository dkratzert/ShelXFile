# Shelxfile
<a href="https://repology.org/project/python:shelxfile/versions">
    <img src="https://repology.org/badge/vertical-allrepos/python:shelxfile.svg" alt="Packaging status" align="right">
</a>

[![Unit tests](https://github.com/dkratzert/ShelXFile/actions/workflows/unit-tests.yml/badge.svg)](https://github.com/dkratzert/ShelXFile/actions/workflows/unit-tests.yml)
![Contributions](https://img.shields.io/badge/contributions-welcome-blue)

This is a full implementation of the SHELXL[[1](https://github.com/dkratzert/Shelxfile/blob/master/README.md#references)] file syntax. Additionally it is able to edit SHELX properties using Python.
The implementation is Python3-only and supports SHELXL after 2017 (You should not use old versions anyway).
Shelxfile may eventually become a new heart of DSR[[2](https://github.com/dkratzert/Shelxfile/blob/master/README.md#references)] and is already used as file parser in StructureFinder[[3](https://github.com/dkratzert/Shelxfile/blob/master/README.md#references)].

Shelxfile always keeps the file order intact. Every SHELX instruction like DFIX or an atom is stored as an class object in the list Shelxfile.\_reslist. When writing the Shelxfile content to disk, it wites the \_reslist content to disk.

Shelxfile tries to detect all possible syntax errors that SHELXL would not like either. If Shelxfile.DEBUG is True, more output about syntax and other errors are printed out. Otherwise, the parser is quiet except for really severe errors like a missing unit cell.

Not every part of Shelxfile is complete, for example it will not recognize if you add restraints with atom names that are not in the SHELX file. Please help me improving it!

**Source Code**

[You can find the ShelXfile source code at GitHub](https://github.com/dkratzert/ShelXFile).

Examples:


```python
pip install shelxfile

>>> from shelxfile import Shelxfile
>>> shx = Shelxfile(verbose=True) # or debug=True, debug will halt on errors.
>>> shx.read_file('src/tests/resources/p21c.res')  # or .read_string() 
>>> shx.cell
CELL 0.71073 10.5086 20.9035 20.5072 90 94.13 90

>>> list(shx.cell)
[10.5086, 20.9035, 20.5072, 90.0, 94.13, 90.0]

>>> shx.cell.volume
4493.047384590458

>>> shx.cell.a
10.5086

>>> shx.to_cif('test.cif')  
# Writes a CIF file from the content of p21c.res

# You can overwrite any parameter in a shelx file:
>>> shx.plan
PLAN 20

>>> shx.plan.npeaks
20

>>> shx.plan.set('PLAN 30')
>>> shx.plan
PLAN 30

>>> shx.atoms
O1     3    0.074835    0.238436    0.402457   -31.00000    0.01579    0.03095      0.01852   -0.00468   -0.00210    0.01153
C1     1    0.028576    0.234542    0.337234   -31.00000    0.02311    0.03617      0.01096   -0.01000    0.00201    0.00356
C2     1    0.121540    0.194460    0.298291   -31.00000    0.02960    0.04586      0.01555   -0.00485   -0.00023    0.01102
F
...

>>> shx.atoms.hydrogen_atoms
[Atom ID: 81, Atom ID: 88, Atom ID: 95, ... ]

>>> shx.atoms.hydrogen_atoms[1].name
'H32'

>>> shx.atoms.n_hydrogen_atoms
24

# Atoms with a riding model e.g. hydrogen atom riding on a carbon atom:
>>> shx.atoms.riding_atoms
[Atom ID: 81, Atom ID: 88, Atom ID: 95, ... ]

>>> a = shx.atoms.get_atom_by_name('F1_2')  # Atom F1 in residue 2
>>> a
Atom ID: 258  # <- The Atom ID is the index number in the Shelxfile._reslist list

>>> str(a)
'F1    4    0.245205    0.192674    0.649231   -21.00000    0.05143    0.03826    0.03193   -0.00579   -0.01865   -0.00485'

>>> a.to_isotropic()
>>> str(a)
'F1    4    0.245205    0.192674   0.649231  -21.00000    0.04000'

>>> a.position  # position in the SHELX .res file
273

>>> str(shx._reslist[273])  # In regular code, do not access shx._reslist directly!
'F1    4    0.245205    0.192674   0.649231  -21.00000    0.04000'

>>> a.name
'F1'

>>> a.element
'F'

# Introduce a new element
>>> a.element = 'Na'
>>> shx.sfac_table
SFAC C  H  O  F  Al  Ga  Na

>>> a.resinum
2

>>> a.part
2

>>> shx.sfac2elem(4)
'F'

>>> shx.elem2sfac('F')
4

>>> a.find_atoms_around(dist=2.0, only_part=1)
[Atom ID: 254, Atom ID: 256, Atom ID: 260]  # Found some atoms 

>>> [str(x) for x in a.find_atoms_around(dist=2.2, only_part=2)]
['C2     1    0.192984    0.140449    0.621265   -21.00000    0.04315    0.02747      0.02385    0.00686   -0.00757    0.00126', 
'F2     4    0.264027    0.090306    0.642441   -21.00000    0.06073    0.04450      0.03972    0.01630   -0.01260    0.01460', 
'F3     4    0.078582    0.131920    0.643529   -21.00000    0.05691    0.04955      0.03374    0.01040    0.01881    0.00375']

>>> a.cart_coords
[1.617897551082389, 4.027560959000001, 13.279336538026433]

>>> a.frac_coords
[0.245205, 0.192674, 0.649231]

>>> a.occupancy
1.0

>>> a.sfac_num
4

>>>[x for x in a.find_atoms_around(dist=2.5, only_part=2)]
[Atom ID: 254, Atom ID: 256, Atom ID: 260, Atom ID: 262]

>>> for x in a.find_atoms_around(dist=2.5, only_part=2):
...    x.delete()

>>> [x for x in a.find_atoms_around(dist=2.5, only_part=2)]
[]  # Atoms are now deleted from shx._reslist.

>>> shx.restraints
SADI_CCF3 0.02 C1 C2 C1 C3 C1 C4
SADI_CCF3 0.02 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4
SADI_CCF3 0.04 C2 C3 C3 C4 C2 C4
SADI_CCF3 0.04 O1 C2 O1 C3 O1 C4
SADI_CCF3 0.04 F1 F2 F2 F3 F3 F1 F4 F5 F5 F6 F6 F4 F7 F8 F8 F9 F9 F7
...

>>> shx.restraints[1]
SADI_CCF3 0.02 C1 C2 C1 C3 C1 C4

str(shx.restraints[1])
'SADI_CCF3 0.02 C1 C2 C1 C3 C1 C4'

shx.restraints[1].residue_class
'CCF3'

# The residue class 'CCF3' has three residues with these numbers: 
>>> shx.restraints[1].residue_number
[4, 1, 2]

# The esd of the SADI restraint:
shx.restraints[1].s
0.02

```

Writes current shx object to test.ins
All lines in Shelxfile._reslist get wrapped after 79 characters with " =\n " as
specified by SHELXL during the file writing.

```python
>>> shx.write_shelx_file('test.ins')
```
No matter if you loaded a .res or .ins file with refine(), SHELXL refines the structure of the Shelxfile() object:

```python
>>> shx.insert_anis()
>>> shx.refine(2)

 Running SHELXL with "/usr/local/bin/shelxl -b3000 /Users/daniel/GitHub/Shelxfile/tests/p21c" and "L.S. 2"
 wR2 =  0.1143 before cycle   1 for   10786 data and    945 /    945 parameters
 wR2 =  0.1025 before cycle   2 for   10786 data and    945 /    945 parameters
 wR2 =  0.1006 before cycle   3 for   10786 data and      0 /    945 parameters
 SHELXL Version 2018/3
```
```python
# Symmcards that are implied by lattice symmetry are generated on-the-fly:
>>> shx.symmcards
| 1  0  0|   | 0.0|
| 0  1  0| + | 0.0|
| 0  0  1|   | 0.0|

|-1  0  0|   | 0.0|
| 0 -1  0| + | 0.0|
| 0  0 -1|   | 0.0|

|-1  0  0|   | 0.0|
| 0  1  0| + | 0.5|
| 0  0 -1|   | 0.5|

| 1  0  0|   | 0.0|
| 0 -1  0| + |-0.5|
| 0  0  1|   |-0.5|

# Complete or "grow" structures with higher symmetry: 
>>> shx = Shelxfile('./tests/p-31c.res')
len(shx.atoms)
88
p = shx.grow()
len(p)
208

# The (bond) angle between three atoms:
>>> at1 = shx.atoms.get_atom_by_name('O1_4')
>>> at2 = shx.atoms.get_atom_by_name('C1_4')
>>> at3 = shx.atoms.get_atom_by_name('C2_4')
>>> shx.atoms.angle(at1, at2, at3)
109.688123

# The torsion angle between four atoms:
>>> at1 = shx.atoms.get_atom_by_name('O1')
>>> at2 = shx.atoms.get_atom_by_name('C1')
>>> at3 = shx.atoms.get_atom_by_name('C2')
>>> at4 = shx.atoms.get_atom_by_name('F1')
>>> shx.atoms.torsion_angle(at1, at2, at3, at4)
74.095731

```

and many more...

## References
[1] http://shelx.uni-goettingen.de/, G. M. Sheldrick, Acta Cryst. (2015). C71, 3-8.
https://doi.org/10.1107/S2053229614024218

[2] https://github.com/dkratzert/DSR

[3] https://github.com/dkratzert/StructureFinder

