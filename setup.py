#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup(
    name='densitypeakclustering',
    version='1.0.0',
    description='Implementation of the Clustering by fast search and find of density peaks algorithm',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Pieter Goltstein',
    author_email='xpieter@mac.com',
    url='https://github.com/pgoltstein/densitypeakclustering',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scikit-learn',
        'matplotlib',
    ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
    ],
    python_requires='>=3.7',
    include_package_data=True,
)