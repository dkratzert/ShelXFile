"""Generate and compare corpus baseline snapshots.

The snapshot is the safety net for the atom/restraint-linking refactor
(plan §E6).  For every structure in the corpus it records:

* the number of parse warnings (``Shelxfile.restraint_errors``),
* a hash of ``Shelxfile.dumps()``,
* whether ``dumps()`` is idempotent (invariant **I2**),
* the exception type, if parsing raised.

Only hashes and relative paths are stored, never structural data, so the
snapshot is safe to commit.  Pass ``--hash-paths`` if even the file names
are sensitive; paths are then replaced by a stable digest.

Usage::

    python -m tests.tools.corpus_snapshot write D:/frames/Workordner
    python -m tests.tools.corpus_snapshot check D:/frames/Workordner
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tempfile
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT / 'src') not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / 'src'))

from shelxfile import Shelxfile  # noqa: E402

DEFAULT_SNAPSHOT = Path(__file__).resolve().parent / 'corpus_baseline.json'


def iter_structures(root: Path) -> list[Path]:
    return sorted(
        [p for p in root.rglob('*.res') if p.is_file()]
        + [p for p in root.rglob('*.ins') if p.is_file()]
    )


def _digest(text: str) -> str:
    return hashlib.sha256(text.encode('utf-8', 'replace')).hexdigest()[:16]


def probe(path: Path) -> dict[str, Any]:
    """Parse *path* and return its baseline record."""
    try:
        first = Shelxfile()
        first.read_file(str(path))
        dumped = first.dumps()
    except Exception as exc:  # noqa: BLE001 - recording, not handling
        return {'error': type(exc).__name__}

    record: dict[str, Any] = {
        'warnings': len(first.restraint_errors or []),
        'dump': _digest(dumped),
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        roundtrip = Path(tmpdir) / path.name
        roundtrip.write_text(dumped, encoding='utf-8', errors='replace')
        try:
            second = Shelxfile()
            second.read_file(str(roundtrip))
            record['idempotent'] = second.dumps() == dumped
        except Exception as exc:  # noqa: BLE001
            record['idempotent'] = False
            record['reparse_error'] = type(exc).__name__
    return record


def build(root: Path, hash_paths: bool = False) -> dict[str, Any]:
    entries: dict[str, Any] = {}
    for path in iter_structures(root):
        key = path.relative_to(root).as_posix()
        if hash_paths:
            key = _digest(key)
        entries[key] = probe(path)
    non_idempotent = sorted(
        k for k, v in entries.items() if v.get('idempotent') is False
    )
    return {
        'version': 1,
        'hashed_paths': hash_paths,
        'count': len(entries),
        'non_idempotent': non_idempotent,
        'errors': sorted(k for k, v in entries.items() if 'error' in v),
        'entries': entries,
    }


def _summarise(snapshot: dict[str, Any]) -> str:
    total = snapshot['count']
    bad = len(snapshot['non_idempotent'])
    err = len(snapshot['errors'])
    return (f'{total} structures, {total - bad - err} idempotent, '
            f'{bad} non-idempotent, {err} parse errors')


def compare(baseline: dict[str, Any], current: dict[str, Any]) -> list[str]:
    """Return human-readable regressions of *current* against *baseline*."""
    problems: list[str] = []
    base_entries = baseline['entries']
    for key, now in current['entries'].items():
        before = base_entries.get(key)
        if before is None:
            continue  # new file: not a regression
        if 'error' in now and 'error' not in before:
            problems.append(f'{key}: now fails to parse ({now["error"]})')
            continue
        if now.get('warnings', 0) > before.get('warnings', 0):
            problems.append(
                f'{key}: parse warnings {before["warnings"]} -> {now["warnings"]}'
            )
        if before.get('idempotent') and not now.get('idempotent', True):
            problems.append(f'{key}: dumps() no longer idempotent')
    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['write', 'check'])
    parser.add_argument('corpus', type=Path)
    parser.add_argument('--snapshot', type=Path, default=DEFAULT_SNAPSHOT)
    parser.add_argument('--hash-paths', action='store_true')
    args = parser.parse_args(argv)

    if not args.corpus.is_dir():
        parser.error(f'corpus directory not found: {args.corpus}')

    current = build(args.corpus, hash_paths=args.hash_paths)

    if args.action == 'write':
        args.snapshot.write_text(json.dumps(current, indent=1), encoding='utf-8')
        print(f'wrote {args.snapshot}: {_summarise(current)}')
        return 0

    if not args.snapshot.is_file():
        print(f'no baseline at {args.snapshot}; run "write" first', file=sys.stderr)
        return 2
    baseline = json.loads(args.snapshot.read_text(encoding='utf-8'))
    problems = compare(baseline, current)
    print(f'current: {_summarise(current)}')
    if problems:
        print(f'{len(problems)} regression(s):', file=sys.stderr)
        for line in problems[:50]:
            print(f'  {line}', file=sys.stderr)
        return 1
    print('no regressions against baseline')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
