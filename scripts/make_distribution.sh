cd ..
rm -Rf dist
rm -Rf shelxfile.egg-info
pip install --upgrade setuptools
python3 -m build
# python3 -m twine upload --repository pypi dist/*