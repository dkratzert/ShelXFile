"""
Generate a grown (complete molecule) .res file from tests/resources/p-31c.res
and write it to ./test_grow.res.

The script calls :meth:`Shelxfile.write_grown_file` which:
  1. Reads p-31c.res with ShelXFile.
  2. Calls shx.grow() to complete the molecular fragments using crystal symmetry.
  3. Changes the space group to P1 (LATT -1, no SYMM cards) so that SHELXL /
     fastmolwidget does not re-apply symmetry to the already-grown atoms.
  4. Preserves disorder parts (PART cards) so the bond graph has no wrong bonds.
  5. Adds a REM warning that the file is not suited for refinement.
  6. Writes the result to ./test_grow.res.
"""

from pathlib import Path
from shelxfile import Shelxfile

# ── 1. load the file ──────────────────────────────────────────────────────────
shx = Shelxfile(verbose=True)
shx.read_file('tests/resources/p-31c.res')

print(f"Asymmetric unit: {shx.atoms.number} atoms "
      f"({len(shx.atoms.q_peaks)} Q-peaks)")

# ── 2. grow and write ─────────────────────────────────────────────────────────
out = Path('./test_grow.res')
shx.write_grown_file(out)

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
