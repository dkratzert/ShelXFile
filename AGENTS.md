# AGENTS.md — ShelXFile Codebase Guide

* Do not write code immediately. Plan your approach first.
* Focus on understanding the architecture and conventions before making changes.
* When in doubt, refer back to this guide or ask for clarification.
* When adding new features, ensure they fit within the existing architecture and follow established conventions.
* When writing tests, use the existing test files in `tests` as fixtures instead of creating new ones, unless absolutely necessary.

## Project Overview
ShelXFile is a Python library for parsing, editing, and writing SHELXL crystallographic structure files (`.res`/`.ins`). The public API is a single class: `from shelxfile import Shelxfile`.

## Architecture: The `_reslist`

**Everything** revolves around `Shelxfile._reslist` — an ordered list that mirrors the `.res` file line-by-line. During parsing, string lines are replaced **in-place** with card objects (e.g. `CELL`, `Atom`, `DFIX`). When writing, `str()` is called on each item and written back to disk. This preserves the original file order exactly.

- **Read**: `shx.read_file('file.res')` → populates `_reslist`, then `parse_cards()` replaces strings with objects.
- **Write**: `shx.write_shelx_file('out.ins')` → iterates `_reslist`, skips indices in `shx.delete_on_write`.
- **Never access `_reslist` directly in user/test code** — use the API (`shx.atoms`, `shx.restraints`, etc.).

## Module Map

| Path | Responsibility |
|---|---|
| `shelxfile/shelx/shelx.py` | `Shelxfile` class — parser, high-level API (`grow()`, `pack()`, …) |
| `shelxfile/shelx/cards.py` | All SHELX instruction classes (`CELL`, `DFIX`, `Restraint`, …) |
| `shelxfile/atoms/atom.py` | `Atom` class — fractional/Cartesian coords, occupancy, SFAC, U-value chain (`ucif`, `ustar`, `u_cart`, `ueq`, `Uiso`) |
| `shelxfile/atoms/atoms.py` | `Atoms` container — iteration, lookup, geometry methods, `conntable` property |
| `shelxfile/shelx/sdm.py` | SDM (Shortest Distance Matrix) — `calc_sdm()`, `packer()`, `pack_unit_cell()`; optional C++ fast path via `sdm_cpp` |
| `shelxfile/misc/misc.py` | Parse error classes, `wrap_line`, `build_conntable`, `frac_to_cart`, `cart_to_frac` |
| `shelxfile/misc/dsrmath.py` | `Array`, `OrthogonalMatrix`, crystallographic math; also re-exports `frac_to_cart` and `cart_to_frac` |
| `shelxfile/misc/elements.py` | Element data tables, `get_radius_from_element()` |
| `shelxfile/refine/refine.py` | Thin wrapper that calls the external `shelxl` binary |
| `shelxfile/cif/cif_write.py` | CIF export using a Jinja-style template |
| `shelxfile/version.py` | **Single source of version**: `VERSION = '23'` |

## Key Conventions

* When making major changes. Add them to the documentation in README.md.

### Card parsing pattern
- Singleton cards use `self._assign_card(CardClass(self, spline), line_num)` → stored on a named attribute (`shx.cell`, `shx.wght`).
- List cards use `self._append_card(self.restraints, CardClass(self, spline), line_num)` → appended to a list attribute.

### Atom identification
`Shelxfile.is_atom(line)` returns `True` when: first token not in `SHX_CARDS`, ≥5 tokens, field[1] has no `.` (SFAC integer), and coords are all `≤ 4.0`.

### Atom naming
`atom.fullname` = `"C1_0"` (name + `_` + residue number). Residue 0 is the default. Use `shx.atoms.get_atom_by_name('F1_2')` to look up atom F1 in residue 2.

### SFAC/occupancy encoding
SHELX encodes occupancy as `31.00000` (fvar=3, occ=1.0). Use `atom.occupancy` (resolved float) not `atom.sof` (raw SHELX value). `shx.sfac2elem(n)` / `shx.elem2sfac('C')` convert between element symbols and SFAC indices.

### Parser error modes
```python
Shelxfile()            # silent — best for production use
Shelxfile(verbose=True)  # prints warnings
Shelxfile(debug=True)    # halts on first error (for development)
```

### C++ SDM acceleration (optional)
`shelxfile/shelx/sdm.py` tries `import sdm_cpp` first (pybind11 extension in `src/sdm_cpp/`). Falls back to pure Python silently. To build the extension:
```bash
pip install pybind11
pip install -e . --no-build-isolation
# macOS OpenMP: brew install libomp (detected automatically by setup.py)
```
The module-level `HAS_CPP: bool` flag reflects whether the extension is loaded.  
Both `calc_sdm()` and the pure-Python path produce identical results; the C++ version uses OpenMP to parallelise the outer atom loop.

### SDM methods
- `calc_sdm()` — builds the shortest-distance matrix and assigns `molindex` to every atom (Union-Find).
- `packer(sdm, need_symm)` — grows molecules by applying symmetry operations collected by `calc_sdm()`.
- `pack_unit_cell(symmop_indices=None, cart_tolerance=0.2, with_qpeaks=False)` — applies all (or selected) symmetry operations, folds positions to `[0, 1)`, deduplicates, and returns `list[Atom]`. Does **not** require `calc_sdm()` to have been called first.

### Atom U-value chain
`atom.uvals` stores the raw SHELXL values `[U11, U22, U33, U23, U13, U12]` (or `[Uiso, 0, 0, 0, 0, 0]`).  
Derived properties (all computed on-demand using numpy):
- `atom.ucif` — symmetric 3×3 U(cif) matrix
- `atom.ustar` — U(star) = N @ U(cif) @ N.T, N = diag(a\*, b\*, c\*)
- `atom.u_cart` — U(cart) = A @ U(star) @ A.T (A = orthogonalisation matrix)
- `atom.ueq` — equivalent isotropic U = trace(U_cart) / 3 (IUCr definition)
- `atom.Uiso` — for riding hydrogens (`uvals[0] < 0`): `abs(uvals[0]) × pivot.ueq`; otherwise equals `ueq`

### Connectivity table
`shx.atoms.conntable` returns a tuple of `(i, j)` index pairs (into `all_atoms`) representing covalent bonds, computed by `build_conntable()` in `misc/misc.py`.  Disorder parts, negative-part/symmgen fragments, and H–H pairs are excluded.

### Coordinate conversion utilities
`frac_to_cart(frac, cell)` and `cart_to_frac(cart, cell)` live in `shelxfile/misc/misc.py` and are **also re-exported** from `shelxfile/misc/dsrmath.py`, so either import path works:
```python
from shelxfile.misc.misc import frac_to_cart, cart_to_frac
from shelxfile.misc.dsrmath import frac_to_cart, cart_to_frac  # same functions
```

## Developer Workflows

### Running tests
```bash
pytest tests/          # all tests
pytest tests/test_shelx.py  # single file
```
Tests run from the project root; resource files are referenced as `'tests/resources/p21c.res'` (relative paths). The primary fixture file is `tests/resources/p21c.res`.

### Linting / type checking
```bash
ruff check shelxfile/
ty check shelxfile/
```

### Dependency management
The project uses `uv` (`uv.lock` present). Dev deps (`pytest`, `ruff`, `ty`) live under `[dependency-groups] dev` in `pyproject.toml`. Metadata is in `pyproject.toml`; `setup.py` handles only the optional C++ extension.

## Integration Points
- **SHELXL binary** (`shx.refine()`): looks for `shelxl` or `xl` on `PATH`; writes a `.ins` file, runs the binary, then reloads the `.res` output.
- **CIF export** (`shx.to_cif()`): uses `shelxfile/cif/cif_template.tmpl`; a custom template path can be passed.
- **DSR integration**: `REM DSR PUT/REPLACE` lines are collected into `shx.dsrlines` / `shx.dsrline_nums` for use by the DSR fragment-fitting tool.
- **Include files**: `+filename` lines in `.res` files are inlined during `read_file()`; recursive inclusion is detected and raises `ValueError`.

