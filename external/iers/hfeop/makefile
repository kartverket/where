#####################################################################
#                                                                   #
#    Makefile for  test_hf_eop                                      #
#
#    Linux version.
#                                                                   #
#                                                                   #
#####################################################################

PROGRAM = test_hf_eop 

OBJ=\
hfeop_xyu_mod.o \
jd_to_date_1.o \
test_hf_eop.o 


# Compiler stuff.

FC = gfortran
FFLAGS = -c -static -g -ffree-form 

CC:=gcc
CFLAGS:=-c -static -g -m32 

LDFLAGS:= \
-lz -lm -lc 


$(PROGRAM): $(OBJ)	
	$(FC) -o $(PROGRAM) $(OBJ) $(LDFLAGS) -L/usr/lib/gcc-lib/i586-mandrake-linux-gnu/3.3.1 
	@echo done.

#  The following are special dependency lines for compiling.
.f.o: 
	cd $(<D); $(FC) $(FFLAGS) $<
.c.o:
	cd $(<D); $(CC) $(CFLAGS) $<


clean:
	rm -f *.o; rm -f *.a;
