from setuptools import setup,find_packages
import sys, os

setup(name="pymecsim",
      description="Python Cyclic Voltammetry Simulator",
      version='1.0',
      author='Kiran Vaddi',
      author_email='kiranvad@buffalo.edu',
      license='MIT',
      python_requires='>=3.6',
      install_requires=['numpy','scipy', 'pandas'],
      extras_require = {},
      packages=find_packages(),
      long_description=open('Readme.md').read()
)