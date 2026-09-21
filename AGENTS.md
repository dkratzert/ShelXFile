# AGENTS.md — ShelXFile Codebase Guide

## General Principles
- **No Assumptions:** If information or code is missing, ask for it. Do not make assumptions.
- **Refuse over Guessing:** Prefer refusing over guessing. If you don't know how to complete the code, say you don't know.
- **Require Specifications:** If the user asks you to write code, ask for a detailed specification first. Do not write code until you have a detailed specification.
- **Quality Standards:** Follow the rules of 'Refactoring' and 'Clean Code' as described by Martin Fowler and Robert C. Martin, respectively.

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
| `src/shelxfile/shelx/shelx.py` | `Shelxfile` class — parser, high-level API (`grow()`, `pack()`, …) |
| `src/shelxfile/shelx/cards.py` | All SHELX instruction classes (`CELL`, `DFIX`, `Restraint`, …) |
| `src/shelxfile/atoms/atom.py` | `Atom` class — fractional/Cartesian coords, occupancy, SFAC, U-value chain (`ucif`, `ustar`, `u_cart`, `ueq`, `Uiso`) |
| `src/shelxfile/atoms/atoms.py` | `Atoms` container — iteration, lookup, geometry methods, `conntable` property |
| `src/shelxfile/shelx/sdm.py` | SDM (Shortest Distance Matrix) — `calc_sdm()`, `packer()`, `pack_unit_cell()`; optional C++ fast path via `sdm_cpp` |
| `src/shelxfile/edit/document.py` | `ShelxDocument` — the editing façade and **the only object that touches `_reslist`** |
| `src/shelxfile/edit/card_meta.py` | Card classification enums (`CardLifetime`, `AtomGrouping`, `AtomListSemantics`, `AfixDependency`) and the `AtomReferencingCard` mixin |
| `src/shelxfile/edit/cascade.py` | `CascadeEngine` — works out the full consequences of a deletion, then applies them in one go |
| `src/shelxfile/edit/graph.py` | `AtomRestraintGraph` — atom ↔ card links, rebuilt on demand |
| `src/shelxfile/edit/token_resolver.py` | Atom-token grammar: `>`, `<`, `LAST`, `$element`, `_+`/`_-`, `_$n`, residue scoping |
| `src/shelxfile/edit/eqiv_factory.py` | `EqivFactory` — reuse-or-mint `EQIV` symmetry operations |
| `src/shelxfile/edit/line_map.py` | `render()` — text plus the `_reslist` index each line came from |
| `src/shelxfile/edit/reports.py` | `DeletionReport`, `EditReport`, `RenameReport` and their reasons |
| `src/shelxfile/gui/editor_widget.py` | Optional Qt widget. A **pure view** over `ShelxDocument`; holds no SHELXL knowledge |
| `src/shelxfile/misc/misc.py` | Parse error classes, `wrap_line`, `multiline_test`, `build_conntable`, `frac_to_cart`, `cart_to_frac` |
| `src/shelxfile/misc/dsrmath.py` | `Array`, `OrthogonalMatrix`, crystallographic math; also re-exports `frac_to_cart` and `cart_to_frac` |
| `src/shelxfile/misc/elements.py` | Element data tables, `get_radius_from_element()` |
| `src/shelxfile/refine/refine.py` | Thin wrapper that calls the external `shelxl` binary |
| `src/shelxfile/cif/cif_write.py` | CIF export using a Jinja-style template |
| `src/shelxfile/version.py` | **Single source of version**: `VERSION = '23'` |

## Key Conventions

* When making major changes. Add them to the documentation in README.md.

### Card parsing pattern
- Singleton cards use `self._assign_card(CardClass(self, spline), line_num)` → stored on a named attribute (`shx.cell`, `shx.wght`).
- List cards use `self._append_card(self.restraints, CardClass(self, spline), line_num)` → appended to a list attribute.

### Atom identification
`Shelxfile.is_atom(line)` returns `True` when: first token not in `SHX_CARDS`, ≥5 tokens, field[1] has no `.` (SFAC integer), and the coordinates are plausible. A coordinate is plausible when it is a plain fractional value (≤ 4.0) **or** decodes from SHELXL's `10*m + p` form — *"to fix any atom parameter, add 10"*, so `10.666600` is 0.6666 held fixed and `21.000000` is `1.0 × fv2`. The one refused code is `m = 1` with `abs(p) ≥ 1`, because `11.00000` is overwhelmingly a fixed occupancy that slid into a coordinate column. `Atom.coordinates_as_written` re-encodes on write, so a constraint is never silently released.

### Editing (the `edit` layer)
`ShelxDocument` owns every mutation. `Atom.delete()` deliberately does **not** cascade; use `doc.delete_atoms()` for cleanup. Cascades are one-directional (**D-9**): deleting an atom may remove cards, but removing a card never deletes atoms — `delete_restraint_with_atoms()` must be asked for by name.

Two rules are easy to get wrong and are enforced by tests:
* **An empty atom list is not a dead card.** A bare `ISOR`/`SIMU`/`ANIS` means *"all non-hydrogen atoms"* (`GLOBAL_WHEN_EMPTY`); a bare `SADI`/`CONN`/`HTAB` is a different instruction (`DIRECTIVE_WHEN_EMPTY`); `BUMP`/`DEFS` never take atoms at all (`NEVER_NAMES_ATOMS`). So a card is *removed*, never emptied.
* **`atom_semantics` consults the atoms as parsed**, not the current list, so a card edited down to nothing is still recognised as one that named atoms.

Card serialization uses a dirty flag: an untouched card echoes its source line verbatim, and only a card whose atoms changed is regenerated. This applies to both `Restraint` and `Command` — a card that fails to regenerate is left naming a deleted atom.

### Layering (D-8)
`gui → edit → shelx/atoms`, never the reverse. Nothing under `edit/` may import Qt; nothing under `gui/` may touch `_reslist` or make SHELXL parsing decisions. `tests/test_layering.py` enforces both.

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
`src/shelxfile/shelx/sdm.py` tries `from shelxfile import sdm_cpp` first (pybind11 extension built from `src/sdm_cpp/`, installed as `shelxfile/sdm_cpp.*.pyd|.so`). Falls back to pure Python silently. To build the extension:
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
- `atom.ueq` — equivalent isotropic U. Anisotropic atoms: trace(U_cart) / 3 (IUCr definition); isotropic atoms and Q-peaks: `uvals[0]`; riding atoms with a negative `uvals[0]` (SHELXL `-factor`): `abs(uvals[0]) × pivot.ueq`
- `atom.Uiso` — alias of `atom.ueq`

### Connectivity table
`shx.atoms.conntable` returns a tuple of `(i, j)` index pairs (into `all_atoms`) representing covalent bonds, computed by `build_conntable()` in `misc/misc.py`.  Disorder parts, negative-part/symmgen fragments, and H–H pairs are excluded.

### Coordinate conversion utilities
`frac_to_cart(frac, cell)` and `cart_to_frac(cart, cell)` live in `src/shelxfile/misc/misc.py` and are **also re-exported** from `src/shelxfile/misc/dsrmath.py`, so either import path works:
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

Card behaviour is covered by a **registry-driven gate**: `tests/card_catalog.py` holds one hand-written sample per atom-referencing card, each justified by a quoted sentence from the SHELXL manual, and `tests/test_card_coverage.py` fails until a newly added card class has one. Corpus frequency measures prevalence, not importance — a rare card gets the same scrutiny as a common one.

### Corpus tests (opt-in)
Marked `corpus` and skipped unless given data:
```bash
pytest -m corpus --corpus /path/to/structures     # or $SHELXFILE_CORPUS
pytest -m corpus --corpus PATH --corpus-sample 500  # quick subset
```
`tests/test_corpus.py` covers I1 (parse stability) and I2 (`dumps()` idempotence) as **ratchets** — the budgets in `tests/resources/corpus_expectations.json` may only ever be lowered. `tests/test_corpus_edits.py` deletes an atom from every structure and checks I3 (edit locality), I4 (report completeness) and I5 (no semantic escalation). A full sweep takes about six minutes and is held to `--corpus-time-budget`. No corpus content is ever committed, only aggregate counts.

### Linting / type checking
```bash
ruff check src/shelxfile/
ty check src/shelxfile/
```

### Dependency management
The project uses `uv` (`uv.lock` present). Dev deps (`pytest`, `ruff`, `ty`) live under `[dependency-groups] dev` in `pyproject.toml`. Metadata is in `pyproject.toml`; `setup.py` handles only the optional C++ extension.

## Integration Points
- **SHELXL binary** (`shx.refine()`): looks for `shelxl` or `xl` on `PATH`; writes a `.ins` file, runs the binary, then reloads the `.res` output.
- **CIF export** (`shx.to_cif()`): uses `src/shelxfile/cif/cif_template.tmpl`; a custom template path can be passed.
- **DSR integration**: `REM DSR PUT/REPLACE` lines are collected into `shx.dsrlines` / `shx.dsrline_nums` for use by the DSR fragment-fitting tool.
- **Include files**: `+filename` lines in `.res` files are inlined during `read_file()`; recursive inclusion is detected and raises `ValueError`.

