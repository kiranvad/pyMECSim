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
      long_description=open('README.md').read(),
      long_description_content_type="text/markdown",
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
      ],
)