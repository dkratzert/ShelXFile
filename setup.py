from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='shelxfile',
    version='4',
    url='https://github.com/dkratzert/ShelXFile',
    license='Beerware License',
    author='Daniel Kratzert',
    author_email='dkratzert@gmx.de',
    description='A parser for SHELXL results files.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(where=".", exclude=['tests']),
    python_requires=">=3.5",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Environment :: Console",
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
)
