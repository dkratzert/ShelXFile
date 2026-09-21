"""Measure dumps() idempotence over the corpus: dumps(parse(dumps(parse(f)))) == dumps(parse(f))."""
import sys
import tempfile
from collections import Counter
from pathlib import Path

from shelxfile import Shelxfile

root = Path(sys.argv[1])
limit = int(sys.argv[2]) if len(sys.argv) > 2 else 0
files = sorted(p for p in root.rglob('*')
               if p.is_file() and p.suffix.lower() in ('.res', '.ins'))
if limit:
    files = files[:limit]

ok = failed = errors = 0
causes = Counter()
examples = {}
scratch = Path(tempfile.mkdtemp())

for path in files:
    try:
        first = Shelxfile()
        first.read_file(str(path))
        once = first.dumps()
        # Round-trip through a file, exactly as tests/test_corpus.py does:
        # read_file and read_string do not take the same route.
        target = scratch / 'roundtrip.ins'
        target.write_text(once, encoding='utf-8', errors='replace')
        second = Shelxfile()
        second.read_file(str(target))
        twice = second.dumps()
    except Exception as exc:
        errors += 1
        causes[f'EXC {type(exc).__name__}: {exc}'] += 1
        examples.setdefault(f'EXC {type(exc).__name__}: {exc}', str(path))
        continue
    if once == twice:
        ok += 1
        continue
    failed += 1
    a = once.splitlines()
    b = twice.splitlines()
    tag = 'line count %d -> %d' % (len(a), len(b))
    for x, y in zip(a, b):
        if x != y:
            kw = x.split()[0] if x.split() else '<blank>'
            tag = f'{kw}: {x!r} -> {y!r}'
            break
    key = tag.split(':')[0]
    causes[key] += 1
    examples.setdefault(key, f'{path}\n     {tag}')

print(f'files      : {len(files)}')
print(f'idempotent : {ok}')
print(f'differing  : {failed}')
print(f'errors     : {errors}')
print('\ncauses:')
for cause, count in causes.most_common(25):
    print(f'  {count:5}  {cause}')
    print(f'         e.g. {examples[cause]}')
