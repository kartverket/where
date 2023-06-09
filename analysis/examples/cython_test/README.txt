# Simple test/demo of Cython

## Build

To build the examples use the included `setup_cython.py`-script. It can be run
as

```
$ python setup_cython.py
```

to build all cython files. It is also possible to build one or several cython
files explicitly. For instance

```
$ python setup_cython.py plain memview
```

will build `tax_cython_plain.pyx` and `tax_cython_memview.pyx`.


## Run examples

To run the examples use the included `run.py`-script. It will run all versions
on the same testdata. It can be run as

```
$ python run.py
```

which will run on all test files. To run on one or more test files explicitly,
use something like

```
$ python run.py python numpy cython_memview
```

This example will run and compare the running times of `tax_python.py`,
`tax_numpy.py` and `tax_cython_memview.pyx`.
