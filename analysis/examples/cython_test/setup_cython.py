from distutils.core import setup
from Cython.Build import cythonize
import sys


# Hack sys.argv before sending to setup
versions = sys.argv[1:]
sys.argv[1:] = ["build_ext", "--inplace"]

# Read all versions from directory listing
if not versions or "all" in versions:
    import os

    versions = [f[11:].split(".")[0] for f in os.listdir() if f.startswith("tax_cython_") and f.endswith(".pyx")]
    print(versions)

# Setup each version
for version in versions:
    cython_module = "tax_cython_{version}.pyx".format(version=version)
    setup(name="taxes", ext_modules=cythonize(cython_module, language_level=3, annotate=True))
