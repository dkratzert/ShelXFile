# AGENTS.md — ShelXFile Codebase Guide

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
| `shelxfile/shelx/shelx.py` | `Shelxfile` class — parser, high-level API |
| `shelxfile/shelx/cards.py` | All SHELX instruction classes (`CELL`, `DFIX`, `Restraint`, …) |
| `shelxfile/atoms/atom.py` | `Atom` class — fractional/Cartesian coords, occupancy, SFAC |
| `shelxfile/atoms/atoms.py` | `Atoms` container — iteration, lookup, geometry methods |
| `shelxfile/shelx/sdm.py` | SDM (Shortest Distance Matrix) — structure growing, bonding |
| `shelxfile/misc/misc.py` | Parse error classes, `wrap_line`, `build_conntable` |
| `shelxfile/misc/dsrmath.py` | `Array`, `OrthogonalMatrix`, crystallographic math |
| `shelxfile/refine/refine.py` | Thin wrapper that calls the external `shelxl` binary |
| `shelxfile/cif/cif_write.py` | CIF export using a Jinja-style template |
| `shelxfile/version.py` | **Single source of version**: `VERSION = '23'` |

## Key Conventions

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

