"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import path

# Read some info from the where package itself
import where

here = path.abspath(path.dirname(__file__))
exe = where.__executable__

# Get the long description from the relevant file
with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name=where.__name__,
    version=where.__version__,
    description=[s.replace("\n", " ") for s in where.__doc__.strip().split("\n\n")][0],
    long_description=long_description,
    # The project's main homepage.
    url="https://github.com/kartverket/where",
    # Author details
    author="Norwegian Mapping Authority",
    author_email="ann-silje.kirkvik@kartverket.no",
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: Microsoft",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
    ],
    # What does your project relate to?
    keywords="referenceframe vlbi slr gnss",
    # Using find_packages to find all subpackages in where. This is done explicitly because find_packages uses a VERY
    # long time to traverse the directory tree of data (and the exclude option to find_packages only excludes after
    # doing a full traverse).
    packages=["config", "where"] + ["where." + p for p in find_packages(where="where")],
    # List run-time dependencies here.  These will be installed by pip when your project is installed. For an analysis
    # of "install_requires" vs pip's requirements files see: https://packaging.python.org/en/latest/requirements.html
    install_requires=[
        "cython",
        "folium",
        "IPython",
        "jplephem",
        "matplotlib>=3.1.1",
        "midgard>=1.2.0",
        "numpy",
        "pandas",
        "pint",
        "pycurl",
        "scipy",
        "colorama",
        "h5py",
        "netCDF4",
        "python-editor",
        "seaborn",
        "cartopy",
    ],
    # List additional groups of dependencies here (e.g. development dependencies). You can install these using the
    # following syntax, for example:
    #   $ pip install -e .[optional,dev_tools]
    extras_require={
        "optional": [],
        "dev_tools": ["black", "bumpversion", "flake8", "line_profiler", "mypy", "pytest"],
    },
    # If there are data files included in your packages that need to be installed, specify them here.  If using Python
    # 2.6 or less, then these have to be included in MANIFEST.in as well.
    # package_data={
    #     'sample': ['package_data.dat'],
    # },
    # Although 'package_data' is the preferred approach, in some case you may need to place data files outside of your
    # packages. See: http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('my_data', ['data/data_file'])],
    # To provide executable scripts, use entry points in preference to the "scripts" keyword. Entry points provide
    # cross-platform support and allow pip to create the appropriate form of executable for the target platform.
    entry_points={
        "console_scripts": [
            f"{exe}=where.__main__:main",
            f"{exe:profiler}=where.tools.profiler:profiler",
            f"{exe:release}=where_release:main",
            f"{exe:runner}=where.runner:main",
            f"{exe:there}=where.there:main",
            f"{exe:tools}=where.tools:main",
        ],
        "gui_scripts": [],
    },
)
