#!/usr/bin/env python
from setuptools import setup, find_packages

__version__ = "0.0.0"
exec(open("pysmFISH/_version.py").read())

setup(
    name="pysmFISH",
    version=__version__,
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "Click",
        "dask",
        "dask-jobqueue",
        "python-json-logger",
        "scikit-learn",
        "scipy",
        "scikit-image",
        "zarr",
        "sympy",
        "Cython",
        "nd2reader",  # pims based: https://github.com/rbnvrw/nd2reader
        "pynndescent",
        "snappy",
        "python-snappy",
        "pyarrow",
        "openpyxl",
        "napari",
        "papermill",
        "pandas",
        "bokeh",
        "ipykernel",
        "asyncssh",
        "jupyter-server-proxy",
    ],
    entry_points="""
		[console_scripts]
	""",
    # metadata
    author="Simone Codeluppi",
    author_email="simone.codeluppi@ki.se",
    description="pipeline for spatial rna mapping",
    license="MIT",
    url="",
)

# dask[complete]=2021.07.0
# 'ctoolz',
# 'scikit-image==0.17.2'
# python -m pip install 'fsspec>=0.3.3'
