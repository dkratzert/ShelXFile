# Plan: Linking Atoms and Restraints (auto-cleanup on atom deletion)

Status: **Draft / not yet implemented** — written 2026-09-21.

## Problem statement

Currently, deleting an `Atom` (`Atom.delete()`) only removes it from
`Atoms.all_atoms` / `Shelxfile._reslist`. `Restraint` subclasses
(`DFIX`, `SADI`, `SAME`, `FLAT`, `ISOR`, `SIMU`, `DELU`, `RIGU`, `CHIV`,
`NCSY`, `EADP`, `EXYZ`, ...) store referenced atoms only as **plain
strings** in `self.atoms: list[str]`. There is no link back from an
`Atom` to the `Restraint`s that mention it, so deleting an atom leaves
dangling/invalid atom names in restraints. The only existing
consistency check, `Shelxfile._assign_atoms_to_restraints()`, only
*warns* about bad atoms — it does not clean anything up.

Additionally, `RESI`, `AFIX`, and `PART` cards act as **brackets**
around a run of atoms in `_reslist`. If an atom deletion empties such
a bracket, the bracket itself becomes meaningless and should be
removed too (except for the reset markers `AFIX 0`, `PART 0`,
`RESI 0`, which must never be removed).

## Chosen architecture: Option B — Central `AtomRestraintGraph`

Rejected alternative (Option A): storing a `restraints: set[Restraint]`
attribute directly on `Atom` and mutating restraints directly from
`Atom.delete()`. Rejected because it couples `Atom` and `Restraint`
tightly and doesn't cleanly support a GUI layer (Fastmolwidget)
observing/reporting changes.

Chosen: a **central, decoupled registry/service**,
`AtomRestraintGraph`, owned by `Shelxfile`, that:

- Maintains a derived (cache-like) mapping between atom fullnames and
  the restraints that reference them.
- Is rebuilt (`rebuild()`) after parsing and after any bulk
  atom/restraint mutation (`add_atom`, restraint parsing, etc.).
- Is invoked by `Atom.delete()` to clean up affected restraints and
  empty brackets, returning a summary of what was removed so a GUI can
  react (e.g. via a callback/Qt signal) without ShelXFile itself
  depending on Qt.

## Restraint atom-grouping semantics

Restraints differ in how their flat `atoms: list[str]` should be
interpreted when an atom is removed. Verified against the official
SHELXL manual (https://shelx.uni-goettingen.de/shelxl_html.php):

| Restraint | Grouping mode | Notes / min atoms |
|---|---|---|
| `DFIX`, `DANG`, `SADI` | `pairs` | Atoms may be given as literal pairs **or** as a range `A > B` which SHELXL expands to sequential pairs. Drop the whole resolved pair containing the deleted atom. Min 1 pair (2 atoms). |
| `SAME` | `same_fragments` | Two residue-correlated fragment lists (template + target), matched positionally. Removing an atom must also remove its positional counterpart in the other fragment. |
| `ISOR`, `SIMU`, `DELU`, `RIGU`, `FLAT`, `CHIV`, `EADP`, `EXYZ`, `NCSY` | `flat` | Simple list; drop just the deleted atom. Restraint auto-removed if it falls below a per-class `MIN_ATOMS` (e.g. `FLAT`/`CHIV` = 4, others = 2). |
| Any restraint with a residue class/number suffix (`_CCF3`, `_2`, wildcard `_*`) | residue-scoped (orthogonal flag) | Only removed once **every** residue instance of that class is empty — not just the one the deleted atom belonged to. |

Atom tokens themselves can also be:
- **Ranges**: `C1 > C5` → expands to sequential pairs/list.
- **Wildcards**: `C?` (C1..C9), `C*` (all atoms starting with "C").
- **Residue suffixes**: `_201` (residue number), `_CCF3` (residue class).

This range/wildcard expansion is presently a `TODO` in
`Restraint._parse_line` (`# TODO: resolve ranges like SADI_CCF3 O1 > F9`)
and is not yet implemented anywhere in the codebase. This plan
resolves that TODO as a side effect (see `AtomTokenResolver` below).

## New components

### 1. `AtomTokenResolver` (new module)

Expands a restraint's raw atom tokens (ranges, wildcards, residue
suffixes) into a concrete, ordered list of resolved atom fullnames.
Shared by:
- `Shelxfile._assign_atoms_to_restraints` (refactored to use it,
  behavior-preserving),
- `AtomRestraintGraph.rebuild()` (new use).

### 2. `AtomRestraintGraph` (new module, e.g. `shelx/atom_restraint_graph.py`)

- `rebuild()` — re-scans `shx.restraints`, resolves atom tokens via
  `AtomTokenResolver`, and populates atom→restraint / restraint→atom
  maps.
- `on_atom_deleted(atom)` — for every restraint referencing the atom:
  - if `pairs`: drop the whole pair; remove restraint if below
    `MIN_ATOMS`.
  - if `same_fragments`: drop the atom and its positional counterpart;
    remove restraint if either fragment becomes too short.
  - if `flat`: drop the atom; remove restraint if below `MIN_ATOMS`.
  - if residue-scoped: only remove the restraint once all residues of
    that class are empty (checked via `Residues.all_residues`).

  Also triggers bracket cleanup (see below) and returns a summary
  dict, e.g. `{'restraints': [...], 'brackets': [...]}`.

### 3. Bracket-scope cleanup for `RESI` / `AFIX` / `PART`

- Add `RESI.members` (atoms between this `RESI` card and the next
  `RESI` card) and `RESI.is_empty()`.
- Add a shared `scope_is_empty(shx)` helper usable by `AFIX`, `PART`,
  and `RESI` (all three "bracket" a run of atoms in `_reslist`).
- On atom deletion, if `atom.afix`, `atom.part`, or `atom.resi` scope
  becomes empty, remove the opening card and its paired closer from
  `_reslist` (and from `Residues.all_residues` for `RESI`).
- **Never remove** `AFIX 0`, `PART 0`, or `RESI 0` — these are
  resets/defaults, not scoped groups.

### 4. Supporting collection methods

- `Restraints.remove(restraint)` — mirrors existing `append()`;
  removes from `_restraints` and from `shx._reslist`.
- `Residues.remove(resi)` — symmetric removal helper.

### 5. Restraint class metadata

Add class attributes to each `Restraint` subclass in `cards.py`:

```python
class Restraint(Residue):
    ATOM_GROUPING: str = 'flat'   # 'flat' | 'pairs' | 'same_fragments'
    MIN_ATOMS: int = 2