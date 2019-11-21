#!/usr/bin/env python

from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='ngsscriptlibrary',
      version='1.0.5',
      description='Scripts for NGS data-analysis',
      long_description=readme(),
      long_description_content_type='text/markdown',
      author='Martin Haagmans',
      author_email='mahaagmans@gmail.com',
      packages=['ngsscriptlibrary'],
      url='https://www.github.com/martinhaagmans/ngsscriptlibrary',
      license='MIT',
      include_package_data=True,
      install_requires=['PyMySQL>=0.8.1', 
                        'pybedtools>=0.7.10', 
                        'PyVCF>=0.6.8',
                        'pysam>=0.14.1', 
                        'pandas>=0.23.0', 
                        'seaborn>=0.8.1', 
                        'matplotlib>=2.2.2'],
      )
