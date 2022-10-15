####
source venv/bin/activate
python3 -m twine upload --repository shelxfile dist/*
git tag v"$(cat ./setup.cfg | grep version | cut -d' ' -f3)"
echo "Now push the tag!"