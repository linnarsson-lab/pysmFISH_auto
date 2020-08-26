#!/usr/bin/env python
from setuptools import setup, find_packages

__version__ = "0.0.0"
exec(open('pysmFISH/_version.py').read())

setup(
	name="pysmFISH",
	version=__version__,
	packages=find_packages(),
	include_package_data=True,
	install_requires=[
		'Click',
		'dask',
        'dask-jobqueue',
        'python-json-logger',
		'numpy',
		'scikit-learn',
		'scipy',
		'scikit-image',
		'pyyaml',
		'zarr',
        'sympy',
        'Cython',
        'nd2reader', #pims based: https://github.com/rbnvrw/nd2reader
		'napari',
		'prefect',
		'xarray'
	],
	
	# pipeline scripts
	entry_points='''
        [console_scripts]
        pipeline=deploy_pipeline:cli
    ''',

	# metadata
	author="Simone Codeluppi",
	author_email="simone.codeluppi@ki.se",
	description="pipeline for spatial rna mapping",
	license="MIT",
	url="",
)