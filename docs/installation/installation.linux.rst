Installation Guide for Linux
****************************

Installing the flipchem python package requires python (python 3.x is recommended), a C compiler, and a Fortran compiler. Installation also requires the numpy python package to first be installed. Below is a guide for installation.

Install system packages
=======================

Use you linux system package manager, for example `dnf` on Fedora or `aptitude` in Ubuntu, to install::

    gcc python3 python3-devel

You can usually find the OS specific package name with some light googling. If not, feel free to post an issue to github: https://github.com/amisr/flipchem

Install numpy
=============

It is generally a good idea to work within a python virtual environment. This isolates the python environment that you work within from the system python enviroment. That means you are free to install packages without worrying about breaking your system. If you don't already have an environment set up, `here's a nice guide <https://realpython.com/python-virtual-environments-a-primer/>`_.

After activating the python environment, you can install numpy using `pip`::

    pip install numpy

Installing flipchem
===================

Now you can install flipchem to the same python environment using `pip`::

    pip install git+https://github.com/amisr/flipchem.git@v2020.2.0


That's it! Easy peasy.