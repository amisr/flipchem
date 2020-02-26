Installation Guide for Windows
******************************

Installing the flipchem python package requires python (python 3.x is recommended), a C compiler, and a Fortran compiler. Installation also requires the numpy python package to first be installed. Below is a guide for installation.

Install system packages
=======================

On Windows 10, it is recommended that you use the `chocolatey` package manager to install compilers and python. Package managers simplify the management of installed software. One could also install the C and Fortran compilers and python without a package manager, but this will likely be harder than using `chocolatey`. First, `install chocolatey <https://chocolatey.org/>`_.

Next, install the compilers and python. Here we will use miniconda to provide python::

    > choco install mingw32

and::

    > choco install miniconda3 --params="'/AddToPath:1'" --params="'/RegisterPython:1'"

Optional
--------

If you would like to work within a bash shell, then you should also install `git bash`::

    > choco install git.install

and then, open a git-bash shell and run the following command to initialize the conda python environment for bash::

    $ /c/tools/miniconda3/Scripts/conda init bash

After this, you will have to restart the bash console.

Install numpy
=============

Brief aside: It is generally a good idea to work within a python virtual environment. This isolates the python environment that you work within from the system python enviroment. That means you are free to install packages without worrying about breaking your system. If you don't already have an environment set up, `here's a nice guide <https://realpython.com/python-virtual-environments-a-primer/>`_. Since we are using miniconda, which is isolated from the Windows system anyway, this isn't as important, but it is a good habit. You can use `conda` to install the `virtualenv` package::

    conda install virtualenv

To install `numpy`, we'll use `conda`::

    conda install numpy

Installing flipchem
===================

Now you can install flipchem to the same python environment using `pip`::

    pip install git+https://github.com/amisr/flipchem.git@v2020.2.0

If `pip` isn't installed, you can install it using `conda`::

    conda install pip


That's it! Easy peasy.