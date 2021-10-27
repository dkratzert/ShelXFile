call venv\Scripts\activate.bat
del dist\*
python -m build
rem python3 -m twine upload --repository pypi dist\*