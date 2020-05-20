Installation Guide for macOS
****************************

Installing the flipchem python package requires python (python 3.x is recommended), a C compiler, and a Fortran compiler. Installation also requires the numpy python package to first be installed. Below is a guide for installation.

Install system packages
=======================

In macOS it is recommended that you obtain `gcc` and `python3` using a package manager such as `homebrew` or `macports`. If you choose `homebrew`, then you can install `gcc` and `python` with::

    $ brew install gcc python

This will install the python headers (traditionally supplied by a `python-devel` package in linux) needed to compile the C and Fortran source code in `flipchem`. If you encounter trouble, feel free to post an issue to github: https://github.com/amisr/flipchem

Install numpy
=============

It is generally a good idea to work within a python virtual environment. This isolates the python environment that you work within from the system python enviroment. That means you are free to install packages without worrying about breaking your system. If you don't already have an environment set up, `here's a nice guide <https://realpython.com/python-virtual-environments-a-primer/>`_.

After activating the python environment, you can install numpy using `pip`::

    pip install numpy

Installing flipchem
===================

Now you can install flipchem to the same python environment using `pip`::

    pip install git+https://github.com/amisr/flipchem.git@v2020.2.1


That's it! Easy peasy.