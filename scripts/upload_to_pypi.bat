call venv\Scripts\activate.bat
venv\Scripts\pip install twine

venv\Scripts\python -m twine upload --repository pypi dist\*

echo "Now add a Git tag!"