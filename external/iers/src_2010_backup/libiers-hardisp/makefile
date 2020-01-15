#-----------------------------------------------------------------------
#
# Description:  make file for the IERS library for HARDISP subroutines.
#
# Usage:
#
#    To make the library, type:
#
#      make
#
#    To delete all object files, type:
#
#      make clean
#
# This revision:  2010 October 20 
#
# Copright (C) 2008 IERS Conventions Center.  All rights reserved.
#-----------------------------------------------------------------------

# This suite of functions is only compilable by ANSI Fortran compilers -
# give the name of your preferred Fortran compiler and compilation flags
# here.

FC = gfortran 
FF = -c -O -fPIC

#----YOU SHOULDN'T HAVE TO MODIFY ANYTHING BELOW THIS LINE---------

# Name the IERS library in its source location.

IERS_LIB = libiers-hardisp.a

# The list of IERS library object files.

IERS_OBS = ADMINT.o \
           ETUTC.o \
           EVAL.o \
           JULDAT.o \
           LEAP.o \
           MDAY.o \
           RECURS.o \
           SHELLS.o \
           SPLINE.o \
           TOYMD.o \
           TDFRPH.o \
#-----------------------------------------------------------------------

default: $(IERS_LIB)

# Make and install the library.

$(IERS_LIB): $(IERS_OBS)
	ar ru $(IERS_LIB) $?

# Clean up.

clean:
	rm -f $(IERS_OBS)

#-----------------------------------------------------------------------
# The list of object file dependencies.

ADMINT.o : ADMINT.F
	$(FC) $(FF) -o $@ ADMINT.F
ETUTC.o : ETUTC.F 
	$(FC) $(FF) -o $@ ETUTC.F
EVAL.o : EVAL.F 
	$(FC) $(FF) -o $@ EVAL.F
JULDAT.o : JULDAT.F 
	$(FC) $(FF) -o $@ JULDAT.F
LEAP.o : LEAP.F 
	$(FC) $(FF) -o $@ LEAP.F
MDAY.o : MDAY.F 
	$(FC) $(FF) -o $@ MDAY.F
RECURS.o : RECURS.F 
	$(FC) $(FF) -o $@ RECURS.F
SHELLS.o : SHELLS.F 
	$(FC) $(FF) -o $@ SHELLS.F
SPLINE.o : SPLINE.F 
	$(FC) $(FF) -o $@ SPLINE.F
TOYMD.o : TOYMD.F 
	$(FC) $(FF) -o $@ TOYMD.F
TDFRPH.o : TDFRPH.F 
	$(FC) $(FF) -o $@ TDFRPH.F
#-----------------------------------------------------------------------
