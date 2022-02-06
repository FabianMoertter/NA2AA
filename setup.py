from setuptools import setup, find_packages

setup(
    name="na2aa",
    version="0.1.0",
    descripton="Convert DNA sequences to longest amino acide sequence for all six reading frames",
    author="Fabian Moertter",
    author_email="fabian.moertter@gmx.net",
    url="https://github.com/FabianMoertter/na2aa",
    packages=find_packages(include=["na2aa"]),
    install_requires=["biopython", "pandas"],
    entry_points={"console_scripts": ["na2aa=na2aa.na2aa:main"]},
)
