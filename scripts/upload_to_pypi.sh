####
source venv/bin/activate
python3 -m twine upload --repository shelxfile dist/*
#VERSION="$(cat ./setup.cfg | grep version | cut -d' ' -f3)"
#git tag v$VERSION
echo "Tag pushed, a release should appear on https://test.pypi.org/project/shelxfile"