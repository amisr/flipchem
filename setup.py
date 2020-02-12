"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
Based on: https://github.com/pypa/sampleproject/master/setup.py
"""

from __future__ import absolute_import
import os
import re
from codecs import open
from setuptools import find_packages
# -*- coding: utf-8 -*-
# Check if we can even import numpy, if not, provide a more helpful
# exception message to the user than what it typically provided.
try:
    from numpy.distutils.core import setup, Extension
except ImportError as e:
    text = "There was a problem importing numpy. Do you have it installed?"
    text += "\nImportError Exception:\n%s" % str(e)
    raise Exception(text)


here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# Get version number from __init__.py
regex = "(?<=__version__..\s)\S+"
with open(os.path.join(here,'flipchem/__init__.py'),'r', encoding='utf-8') as f:
    text = f.read()
match = re.findall(regex,text)
version = match[0].strip("'")

# Use mingw32 by default on windows
if os.name == 'nt':
    sfn = os.path.join(os.path.dirname(__file__), 'setup.cfg')
    with open(sfn, 'w') as f:
        f.write("\n[build_ext]\ncompiler = mingw32")

# FLIPCHEM EXTENSION
flipchem_sources = ["src/flipchem/flipchem.pyf",
                    "src/flipchem/flipchem.f"]
flipchem_ext = Extension(name = 'flipchem.ext._f_flipchem',
                         sources = flipchem_sources,
                         extra_f90_compile_args=['--std=legacy','-finit-local-zero','-fno-automatic'],
                         extra_f77_compile_args=['--std=legacy','-finit-local-zero','-fno-automatic'],
                         )

# MSIS EXTENSION
msis_sources = ['src/nrlmsise00/_swigmsis_wrap.c',
                'src/nrlmsise00/nrlmsise-00.c',
                'src/nrlmsise00/nrlmsise-00_data.c']
msis_ext = Extension(name = 'flipchem.ext._c_msis',
                     sources = msis_sources)


packages = find_packages(exclude=['contrib', 'docs', 'tests'])


if __name__ == "__main__":

    setup(
        name='flipchem',

        # Versions should comply with PEP440.  For a discussion on single-sourcing
        # the version across setup.py and the project code, see
        # https://packaging.python.org/en/latest/single_source_version.html
        version=version,

        description='The ion composition model from the Field Line Interhemispheric Plasma (FLIP) ionosphere model wrapped in python.',
        long_description=long_description,

        # The project's main homepage.
        url='',

        # Author details
        author='Ashton S. Reimer',
        author_email='ashton.reimer@sri.com',

        # Choose your license
        license='MIT',

        # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
        classifiers=[
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            'Development Status :: 3 - Alpha',

            # Indicate who your project is intended for
            'Intended Audience :: Developers',
            'Topic :: Software Development :: Build Tools',

            # Pick your license as you wish (should match "license" above)
            'License :: OSI Approved :: MIT License',

            # Specify the Python versions you support here. In particular, ensure
            # that you indicate whether you support Python 2, Python 3 or both.
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.7',
        ],

        # What does your project relate to?
        keywords='flip flipchem ionosphere composition',

        # You can just specify the packages manually here if your project is
        # simple. Or you can use find_packages().
        packages=packages,

        # Alternatively, if you want to distribute just a my_module.py, uncomment
        # this:
        #   py_modules=["my_module"],

        # List run-time dependencies here.  These will be installed by pip when
        # your project is installed. For an analysis of "install_requires" vs pip's
        # requirements files see:
        # https://packaging.python.org/en/latest/requirements.html
        install_requires=['numpy'],

        ext_modules = [flipchem_ext,msis_ext],

        # List additional groups of dependencies here (e.g. development
        # dependencies). You can install these using the following syntax,
        # for example:
        # $ pip install -e .[dev,test]
        # extras_require={
        #     'dev': ['check-manifest'],
        #     'test': ['coverage'],
        # },

        # If there are data files included in your packages that need to be
        # installed, specify them here.  If using Python 2.6 or less, then these
        # have to be included in MANIFEST.in as well.
        package_data={'flipchem': ['dat/*'],
                     },
        include_package_data=True,

        # Although 'package_data' is the preferred approach, in some case you may
        # need to place data files outside of your packages. See:
        # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
        # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
        # data_files=[('my_data', ['data/data_file'])],

        # To provide executable scripts, use entry points in preference to the
        # "scripts" keyword. Entry points provide cross-platform support and allow
        # pip to create the appropriate form of executable for the target platform.
        # entry_points={
        #     'console_scripts': [
        #         'sample=sample:main',
        #     ],
        # },
    )
