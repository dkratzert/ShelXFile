rem cd ..
call venv\Scripts\activate.bat
venv\Scripts\pip install twine
rem cd ..
venv\Scripts\python -m twine upload --repository pypi dist\*