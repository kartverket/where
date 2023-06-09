from distutils.core import setup
from Cython.Build import cythonize
import sys


# Hack sys.argv before sending to setup
versions = sys.argv[1:]
sys.argv[1:] = ["build_ext", "--inplace"]

# Read all versions from directory listing
if not versions or "all" in versions:
    import os

    versions = [f.split(".")[0] for f in os.listdir() if f.endswith(".pyx")]
    print(versions)

# Setup each version
for version in versions:
    cython_module = "{version}.pyx".format(version=version)
    setup(name="cytest", ext_modules=cythonize(cython_module, language_level=3, annotate=True))
