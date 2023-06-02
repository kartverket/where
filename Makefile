# Makefile for simple installation of the Where python project
#
# Authors:
# --------
#
# * Geir Arne Hjelle <geir.arne.hjelle@kartverket.no>

# Programs and directories
# Note from 'dahmic': It seems to be that 'f2py' could not be installed at on a Rocky 8.7 operativ system. 'numpy' has 
#                     overtaken the 'f2py' project. Therefore I used 'F2PY = python -m numpy.f2py' instead of 
#                     'F2PY=f2py.
F2PY = f2py
F2PYEXTENSION = .cpython-36m-x86_64-linux-gnu.so

EXTDIR = $(CURDIR)/where/ext
SOFADIR = $(CURDIR)/external/sofa/src
IERSDIR = $(CURDIR)/external/iers/src
HFEOPDIR = $(CURDIR)/external/iers/hfeop
GPT2WDIR = $(CURDIR)/external/gpt2w/src

# Define phony targets (targets that are not files)
.PHONY: all develop install cython doc test typing format external sofa iers_2010 hf_eop gpt2w

# Everything needed for installation
all:	external cython develop local_config

# Install in developer mode (no need to reinstall after changing source)
develop:
	python -m pip install -e .[optional,dev_tools]

# Install on server: Developer mode for ease of switching tags
server:
	python -m pip install -e .[optional]

# Regular install, freezes the code so must reinstall after changing source code
install:
	python -m pip install --user .[optional]

# Set up local configuration files
local_config:
	python config_wizard.py config/where_local.conf --ignore-existing
	python config_wizard.py config/there_local.conf --ignore-existing
	python config_wizard.py config/files_local.conf --ignore-existing
	python config_wizard.py config/trf_local.conf --ignore-existing

# Compile Cython files
cython:
	python setup_cython.py

# Run tests
test:
	pytest -m quick --cov=where

test_system:
	pytest -m system_test --durations=0

test_vlbi:
	pytest -m vlbi --durations=0

test_slr:
	pytest -m slr --durations=0

test_sisre:
	pytest -m sisre --durations=0

typing:
	mypy --ignore-missing-imports where

# Reformat code
black:
	black --line-length 119 where/ analysis/ tests/ *.py


######################################################################

# External libraries
external:	sofa iers_2010 gpt2w hf_eop

# SOFA
sofa:	$(EXTDIR)/sofa$(F2PYEXTENSION)

$(SOFADIR)/sofa.pyf:
	python download.py sofa

$(EXTDIR)/sofa$(F2PYEXTENSION):	$(SOFADIR)/sofa.pyf $(shell find $(SOFADIR) -type f)
	( cd $(EXTDIR) && $(F2PY) -c $(SOFADIR)/sofa.pyf $(SOFADIR)/*.for )

# IERS
iers_2010:	$(EXTDIR)/iers_2010$(F2PYEXTENSION)

$(IERSDIR)_2010/iers_2010.pyf:
	python download.py iers_2010

$(EXTDIR)/iers_2010$(F2PYEXTENSION):	$(IERSDIR)_2010/iers_2010.pyf $(IERSDIR)_2010/libiers-dehant/libiers-dehant.a $(IERSDIR)_2010/libiers-hardisp/libiers-hardisp.a $(shell find $(IERSDIR)_2010 -type f)
	( cd $(EXTDIR) && \
          $(F2PY) -c $(IERSDIR)_2010/iers_2010.pyf $(IERSDIR)_2010/*.F \
                  -liers-dehant -L$(IERSDIR)_2010/libiers-dehant -liers-hardisp -L$(IERSDIR)_2010/libiers-hardisp )

$(IERSDIR)_2010/libiers-dehant/libiers-dehant.a:	$(shell find $(IERSDIR)_2010/libiers-dehant -type f -name *.F)
	( cd $(IERSDIR)_2010/libiers-dehant && make && make clean )

$(IERSDIR)_2010/libiers-hardisp/libiers-hardisp.a:	$(shell find $(IERSDIR)_2010/libiers-hardisp -type f -name *.F)
	( cd $(IERSDIR)_2010/libiers-hardisp && make && make clean )


# HF EOP
hf_eop:         $(EXTDIR)/hf_eop$(F2PYEXTENSION)

$(HFEOPDIR)/hfeop_xyu_mod.pyf:
	python download.py hf_eop

$(EXTDIR)/hf_eop$(F2PYEXTENSION): $(HFEOPDIR)/hfeop_xyu_mod.pyf   $(shell find $(HFEOPDIR) -type f -name *.f90) 
	( cd $(EXTDIR) && \
           $(F2PY) -c $(HFEOPDIR)/hfeop_xyu_mod.pyf $(HFEOPDIR)/hfeop_xyu_mod.f90 ) 


# GPT2w
gpt2w:	$(EXTDIR)/gpt2w$(F2PYEXTENSION)

$(GPT2WDIR)/gpt2w.pyf:
	python download.py gpt2w

$(EXTDIR)/gpt2w$(F2PYEXTENSION):	$(GPT2WDIR)/gpt2w.pyf $(shell find $(GPT2WDIR) -type f -name *.f)
	( cd $(EXTDIR) && $(F2PY) -c $(GPT2WDIR)/gpt2w.pyf $(GPT2WDIR)/*.f )
