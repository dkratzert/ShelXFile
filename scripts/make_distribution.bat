cd ..
del dist\*
python3 -m build
rem python3 -m twine upload --repository pypi dist\*