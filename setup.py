from setuptools import setup, find_packages

setup(
    name='ShelXFile',
    version='1',
    url='https://github.com/dkratzert/ShelXFile',
    license='Beerware License',
    author='Daniel Kratzert',
    author_email='dkratzert@gmx.de',
    description='A parser for SHELXL results files.',
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.5",
)
