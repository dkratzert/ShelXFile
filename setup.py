from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='shelxfile',
    version='1',
    url='https://github.com/dkratzert/ShelXFile',
    license='Beerware License',
    author='Daniel Kratzert',
    author_email='dkratzert@gmx.de',
    description='A parser for SHELXL results files.',
    package_dir={"": "shelxfile"},
    packages=find_packages(where="shelxfile"),
    python_requires=">=3.5",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Environment :: Console",
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "License :: Beerware",
        "Operating System :: OS Independent",
    ],
)
