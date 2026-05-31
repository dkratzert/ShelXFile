# Shelxfile
<a href="https://repology.org/project/python:shelxfile/versions">
    <img src="https://repology.org/badge/vertical-allrepos/python:shelxfile.svg" alt="Packaging status" align="right">
</a>

[![Unit tests](https://github.com/dkratzert/ShelXFile/actions/workflows/unit-tests.yml/badge.svg)](https://github.com/dkratzert/ShelXFile/actions/workflows/unit-tests.yml)
![Contributions](https://img.shields.io/badge/contributions-welcome-blue)

This is a full implementation of the SHELXL[[1](#references)] file syntax. Additionally it is able to edit SHELX properties using Python.
The implementation is Python3-only and supports SHELXL after 2017 (you should not use old versions anyway).
Shelxfile is used as file parser in StructureFinder[[3](#references)].

Shelxfile always keeps the file order intact. Every SHELX instruction like DFIX or an atom is stored as a class object in the list `Shelxfile._reslist`. When writing the Shelxfile content to disk, it writes the `_reslist` content to disk.

Shelxfile tries to detect all possible syntax errors that SHELXL would not like either. Use `verbose=True` during initialization for more output about syntax and other errors. Use `debug=True` to halt on errors. Otherwise, the parser is quiet except for really severe errors like a missing unit cell.

Not every part of Shelxfile is complete, for example it will not recognize if you add restraints with atom names that are not in the SHELX file. Please help me improving it!

**Source Code**

[You can find the ShelXfile source code at GitHub](https://github.com/dkratzert/ShelXFile).


## Installation

```shell
pip install shelxfile
```


## Quick Start

```python
from shelxfile import Shelxfile

shx = Shelxfile(verbose=True)  # or debug=True, debug will halt on errors.
shx.read_file('tests/resources/p21c.res')  # or shx.read_string('...')
```


## Examples

### Unit Cell

```python
>>> shx.cell
CELL 0.71073 10.5086 20.9035 20.5072 90 94.13 90

>>> list(shx.cell)
[10.5086, 20.9035, 20.5072, 90.0, 94.13, 90.0]

>>> shx.cell.volume
4493.047384590458

>>> shx.cell.a
10.5086
```

### CIF Export

```python
>>> shx.to_cif('test.cif')
# Writes a CIF file from the content of p21c.res
# Optionally pass a custom template: shx.to_cif('test.cif', template='my_template.tmpl')
```

### Modifying SHELX Instructions

You can overwrite any parameter in a SHELX file:

```python
>>> shx.plan
PLAN 20

>>> shx.plan.npeaks
20

>>> shx.plan.set('PLAN 30')
>>> shx.plan
PLAN 30
```

### Atoms

```python
>>> shx.atoms
O1     3    0.074835    0.238436    0.402457   ...
C1     1    0.028576    0.234542    0.337234   ...
C2     1    0.121540    0.194460    0.298291   ...
...

>>> len(shx.atoms)
148

>>> shx.atoms.number
148

>>> shx.atoms.hydrogen_atoms
[Atom ID: 134, Atom ID: 141, Atom ID: 148, ...]

>>> shx.atoms.hydrogen_atoms[1].name
'H32'

>>> shx.atoms.n_hydrogen_atoms
24

# Atoms with a riding model (e.g. hydrogen atom riding on a carbon atom):
>>> shx.atoms.riding_atoms
[Atom ID: 134, Atom ID: 141, Atom ID: 148, ...]

# Q-peaks in the file:
>>> shx.atoms.q_peaks
[Atom ID: 328, Atom ID: 329, ...]

# Number of anisotropic/isotropic atoms:
>>> shx.atoms.n_anisotropic_atoms
124
>>> shx.atoms.n_isotropic_atoms
24

# Residue numbers present:
>>> shx.atoms.residues
[0, 1, 2, 3, 4]
```

### Working with Individual Atoms

```python
>>> a = shx.atoms.get_atom_by_name('F1_2')  # Atom F1 in residue 2
>>> a
Atom ID: 258  # The Atom ID is the index number in Shelxfile._reslist

>>> str(a)
'F1    4    0.245205    0.192674    0.649231   -21.00000    0.05143    0.03826      0.03193   -0.00579   -0.01865   -0.00485'

>>> a.name
'F1'

>>> a.element
'F'

>>> a.resinum
2

>>> a.part.n
2

>>> a.sfac_num
4

>>> a.occupancy
0.51904

>>> a.frac_coords
(0.245205, 0.192674, 0.649231)

>>> a.cart_coords
(1.617897551082389, 4.027560959000001, 13.279336538026431)

>>> a.is_hydrogen
False

>>> a.is_isotropic
False

>>> a.atomid  # position in the SHELX .res file (_reslist index)
258

>>> str(shx._reslist[a.atomid])  # In regular code, do not access shx._reslist directly!
'F1    4    0.245205    0.192674    0.649231   -21.00000    0.05143    ...'
```

### Displacement Parameters

```python
>>> c = shx.atoms.get_atom_by_name('C1')

# Equivalent isotropic U (trace of U_cart / 3, IUCr definition):
>>> c.ueq
0.027...

# For riding hydrogen atoms the Uiso is a multiple of the pivot atom's Ueq:
>>> h = shx.atoms.get_atom_by_name('H34')
>>> h.pivot.name       # the carbon H34 rides on
'C34'
>>> h.Uiso             # = 1.2 × C34.ueq  (encoded as -1.2 in the .res file)
0.02956...
>>> h.Uiso == h.pivot.Uiso * 1.2
True

# Full anisotropic U-value chain (numpy arrays):
>>> c.ucif              # 3×3 U(cif) matrix  [U11 U12 U13 / U12 U22 U23 / U13 U23 U33]
array(...)
>>> c.ustar             # U(star) = N @ U(cif) @ N.T,  N = diag(a*, b*, c*)
array(...)
>>> c.u_cart            # U(cart) = A @ U(star) @ A.T,  A = orthogonalisation matrix
array(...)
```

### Connectivity Table

```python
# Pairwise bond connectivity for all atoms in the asymmetric unit:
>>> conn = shx.atoms.conntable   # tuple of (i, j) index pairs
>>> conn[0]
(0, 1)
```

### Modifying Atoms

```python
# Make an atom isotropic:
>>> a.to_isotropic()
>>> str(a)
'F1    4    0.245205    0.192674   0.649231  -21.00000    0.04000'

# Introduce a new element (automatically updates the SFAC table):
>>> a.element = 'Na'
>>> shx.sfac_table
SFAC C  H  O  F  Al  Ga  Na
```

### Adding and Deleting Atoms

```python
# Add a new atom:
>>> shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3], element='C')

# Delete atoms around a given atom:
>>> for x in a.find_atoms_around(dist=2.5, only_part=2):
...     x.delete()
```

### Finding Nearby Atoms

```python
>>> a.find_atoms_around(dist=2.0, only_part=1)
[Atom ID: 239, Atom ID: 241, Atom ID: 245]

>>> [str(x) for x in a.find_atoms_around(dist=2.2, only_part=2)]
['C2     1    0.192984    0.140449    ...', 'F2     4    ...', 'F3     4    ...']
```

### SFAC Table Lookups

```python
>>> shx.sfac2elem(4)
'F'

>>> shx.elem2sfac('F')
4
```

### Restraints

```python
>>> shx.restraints[1]
SADI_CCF3 0.02 C1 C2 C1 C3 C1 C4

>>> str(shx.restraints[1])
'SADI_CCF3 0.02 C1 C2 C1 C3 C1 C4'

>>> shx.restraints[1].residue_class
'CCF3'

# The residue class 'CCF3' has three residues with these numbers:
>>> shx.restraints[1].residue_number
[4, 1, 2]

# The esd of the SADI restraint:
>>> shx.restraints[1].s
0.02
```

### Distances and Angles

```python
# Distance between two atoms (by name):
>>> shx.atoms.distance('O1', 'C1')
1.3505645511659556

# Bond angle between three atoms:
>>> at1 = shx.atoms.get_atom_by_name('O1_4')
>>> at2 = shx.atoms.get_atom_by_name('C1_4')
>>> at3 = shx.atoms.get_atom_by_name('C2_4')
>>> shx.atoms.angle(at1, at2, at3)
109.68812347

# Torsion angle between four atoms:
>>> at1 = shx.atoms.get_atom_by_name('O1')
>>> at2 = shx.atoms.get_atom_by_name('C1')
>>> at3 = shx.atoms.get_atom_by_name('C2')
>>> at4 = shx.atoms.get_atom_by_name('F1')
>>> shx.atoms.torsion_angle(at1, at2, at3, at4)
74.09573117980119
```

### Symmetry Cards

Symmetry cards that are implied by lattice symmetry are generated on-the-fly:

```python
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
```

### Growing Structures

Complete or "grow" structures with higher symmetry:

```python
>>> shx2 = Shelxfile()
>>> shx2.read_file('tests/resources/p-31c.res')
>>> len(shx2.atoms)
88
>>> p = shx2.grow()
>>> len(p)
208
```

### Packing the Unit Cell

`pack()` applies all symmetry operations to the asymmetric unit and folds every
position back into `[0, 1)` fractional coordinates, removing duplicates.
Unlike `grow()`, it does not stitch molecular fragments together — it simply
fills one unit cell.

```python
>>> shx2 = Shelxfile()
>>> shx2.read_file('tests/resources/p-31c.res')
>>> len(shx2.atoms)          # asymmetric unit
88
>>> packed = shx2.pack()
>>> len(packed)              # full unit cell (Z × asymm unit, minus special positions)
304
>>> all(0.0 <= a.x < 1.0 and 0.0 <= a.y < 1.0 and 0.0 <= a.z < 1.0 for a in packed)
True
```

The result is a plain list of `Atom` objects — the `Shelxfile` object itself is
not modified.  Q-peaks can be included with `shx.pack(with_qpeaks=True)`.

### Writing Files

Writes the current `shx` object to a SHELX file.
All lines in `Shelxfile._reslist` get wrapped after 79 characters with `" =\n "` as
specified by SHELXL during the file writing.

```python
>>> shx.write_shelx_file('test.ins')
```

### Optional C++ Acceleration

The SDM (Shortest Distance Matrix), which underlies `grow()` and `pack()`, ships
with an optional C++ extension (`sdm_cpp`) compiled with pybind11 and OpenMP.
When present it is used automatically and can give a **5–10× speedup** on large
structures; if it is absent the pure-Python fallback is used silently.

```bash
# Build the extension (requires a C++17 compiler):
pip install pybind11
pip install -e . --no-build-isolation

# macOS only — OpenMP support:
brew install libomp
```

```python
from shelxfile.shelx.sdm import HAS_CPP
print(HAS_CPP)   # True when the extension is installed
```

### Sum Formula

```python
# Sum formula based on UNIT instruction:
>>> shx.sum_formula
'C0.25 H0.5 O0.75 F1 AL1.25 GA1.5'

# Exact sum formula from all atom occupancies:
>>> shx.sum_formula_exact
'C34 H24 O4 F36 Al1 Ga1'
```

### Residuals from .res File

These values are parsed from `REM` lines written by SHELXL into `.res` files:

```python
>>> shx.R1
0.04

>>> shx.wr2
0.1005

>>> shx.goof
1.016

>>> shx.space_group
'P2(1)/c'

>>> shx.wavelength
0.71073
```

### Refinement (Requires SHELXL)

No matter if you loaded a `.res` or `.ins` file, `refine()` runs SHELXL on the `Shelxfile` object:

```python
>>> shx.insert_anis()
>>> shx.refine(2)

 Running SHELXL with "/usr/local/bin/shelxl -b3000 ..." and "L.S. 2"
 wR2 =  0.1143 before cycle   1 for   10786 data and    945 /    945 parameters
 wR2 =  0.1025 before cycle   2 for   10786 data and    945 /    945 parameters
 wR2 =  0.1006 before cycle   3 for   10786 data and      0 /    945 parameters
 SHELXL Version 2018/3
```


## References
[1] http://shelx.uni-goettingen.de/, G. M. Sheldrick, Acta Cryst. (2015). C71, 3-8.
https://doi.org/10.1107/S2053229614024218

[2] https://github.com/dkratzert/DSR

[3] https://github.com/dkratzert/StructureFinder

