# Makefile for simple installation of the Where python project
#
# Authors:
# --------
#
# * Geir Arne Hjelle <geir.arne.hjelle@kartverket.no>
#
# $Revision: 15244 $
# $Date: 2018-06-01 23:58:09 +0200 (Fri, 01 Jun 2018) $
# $LastChangedBy: hjegei $

# Programs and directories
F2PY = f2py
F2PYEXTENSION = .cpython-36m-x86_64-linux-gnu.so

DOCSDIR = $(CURDIR)/documents/docs
DOCSDIR_WWW = /home/geosat/doc-where
EXTDIR = $(CURDIR)/where/ext
SOFADIR = $(CURDIR)/external/sofa/src
IERSDIR = $(CURDIR)/external/iers/src
GPT2WDIR = $(CURDIR)/external/gpt2w/src

# Define phony targets (targets that are not files)
.PHONY: develop install cython doc test typing format external sofa iers_2010 gpt2w

# Install in developer mode (no need to reinstall after changing source)
develop:
	pip install --user -e .[optional,dev_tools]

# Install on server: Developer mode for ease of switching tags
server:
	pip install -e .[optional]

# Regular install, freezes the code so must reinstall after changing source code
install:
	pip install --user .[optional]

# Compile Cython files
cython:
	python setup_cython.py

# Create documentation
doc:
	( cd $(DOCSDIR) && make html && cp -R build/html $(DOCSDIR_WWW) )

# Run tests
test:
	pytest -m quick --cov=where

test_system:
	pytest -m system_test --durations=0

typing:
	mypy --ignore-missing-imports where

# Reformat code (for now, only see what black would have done)
black:
	black --line-length 119 where/ analysis/ tests/


######################################################################

# External libraries
external:	sofa iers_2010 gpt2w

# SOFA
sofa:	$(EXTDIR)/sofa$(F2PYEXTENSION)

$(EXTDIR)/sofa$(F2PYEXTENSION):	$(shell find $(SOFADIR) -type f)
	( cd $(EXTDIR) && $(F2PY) -c $(SOFADIR)/sofa.pyf $(SOFADIR)/*.for )

# IERS
iers_2010:	$(EXTDIR)/iers_2010$(F2PYEXTENSION)

$(EXTDIR)/iers_2010$(F2PYEXTENSION):	$(IERSDIR)_2010/libiers-dehant/libiers-dehant.a $(IERSDIR)_2010/libiers-hardisp/libiers-hardisp.a $(shell find $(IERSDIR)_2010 -type f)
	( cd $(EXTDIR) && \
          $(F2PY) -c $(IERSDIR)_2010/iers_2010.pyf $(IERSDIR)_2010/*.F \
                  -liers-dehant -L$(IERSDIR)_2010/libiers-dehant -liers-hardisp -L$(IERSDIR)_2010/libiers-hardisp )

$(IERSDIR)_2010/libiers-dehant/libiers-dehant.a:	$(shell find $(IERSDIR)_2010/libiers-dehant -type f -name *.F)
	( cd $(IERSDIR)_2010/libiers-dehant && make && make clean )

$(IERSDIR)_2010/libiers-hardisp/libiers-hardisp.a:	$(shell find $(IERSDIR)_2010/libiers-hardisp -type f -name *.F)
	( cd $(IERSDIR)_2010/libiers-hardisp && make && make clean )

# GPT2w
gpt2w:	$(EXTDIR)/gpt2w$(F2PYEXTENSION)

$(EXTDIR)/gpt2w$(F2PYEXTENSION):	$(shell find $(GPT2WDIR) -type f)
	( cd $(EXTDIR) && $(F2PY) -c $(GPT2WDIR)/gpt2w.pyf $(GPT2WDIR)/*.f )
