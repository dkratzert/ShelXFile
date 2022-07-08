call venv\Scripts\activate.bat
del /Q dist\*
del /q /s shelxfile.egg-info
pip install --upgrade setuptools
venv\Scripts\python -m build
rem python3 -m twine upload --repository pypi dist\*