"""
Generate a grown (complete molecule) .res file from tests/resources/p-31c.res
and write it to ./test_grow.res.

The script:
  1. Reads p-31c.res with ShelXFile.
  2. Calls shx.grow() to complete the molecular fragments using crystal symmetry.
  3. Removes the original asymmetric-unit atoms and the SYMM cards from the
     _reslist (they are no longer needed once the structure is grown).
  4. Inserts the grown atoms immediately before the HKLF card.
  5. Writes the result to ./test_grow.res.
"""

from pathlib import Path
from shelxfile import Shelxfile
from shelxfile.atoms.atom import Atom
from shelxfile.shelx.cards import HKLF, SYMM
from shelxfile.misc.misc import wrap_line

# ── 1. load the file ──────────────────────────────────────────────────────────
shx = Shelxfile(verbose=True)
shx.read_file('tests/resources/p-31c.res')

print(f"Asymmetric unit: {shx.atoms.number} atoms "
      f"({len(shx.atoms.q_peaks)} Q-peaks)")

# ── 2. grow ───────────────────────────────────────────────────────────────────
grown_atoms = shx.grow()
real_grown  = [a for a in grown_atoms if not a.qpeak]
print(f"Grown structure: {len(grown_atoms)} atoms total, "
      f"{len(real_grown)} non-Q-peak atoms")

# ── 3. mark original atoms and SYMM cards for deletion ───────────────────────
for idx, item in enumerate(shx._reslist):
    if isinstance(item, Atom):
        shx.delete_on_write.add(idx)
    elif isinstance(item, SYMM):
        shx.delete_on_write.add(idx)

# ── 4. insert grown atoms just before HKLF ───────────────────────────────────
hklf_idx = next(
    (i for i, x in enumerate(shx._reslist) if isinstance(x, HKLF)),
    len(shx._reslist),
)

# Insert in reverse so that the order is preserved after repeated inserts
for atom in reversed(real_grown):
    shx._reslist.insert(hklf_idx, atom)

print(f"Inserted {len(real_grown)} grown atoms at _reslist position {hklf_idx}")

# ── 5. write ──────────────────────────────────────────────────────────────────
out = Path('./test_grow.res')
shx.write_shelx_file(out)

# Quick sanity check: count atom lines in the output
atom_lines = [
    line for line in out.read_text().splitlines()
    if Shelxfile.is_atom(line)
]
print(f"\nOutput written to: {out.resolve()}")
print(f"Atom lines in output: {len(atom_lines)}")
print(f"\nFirst 5 atom lines:")
for line in atom_lines[:5]:
    print(" ", line.rstrip())
print(f"\nLast 5 atom lines:")
for line in atom_lines[-5:]:
    print(" ", line.rstrip())

