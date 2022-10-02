# run from project root
source venv/bin/activate
rm -Rf dist
rm -Rf shelxfile.egg-info
pip install pip -U
pip install --upgrade setuptools
python3 -m pip install --upgrade build
python3 -m build
# python3 -m twine upload --repository pypi dist/*