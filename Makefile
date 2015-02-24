# ============================================================================
# Name        : Makefile
# Author      : 
# Version     :
# Copyright   : 
# Description : Makefile for Damian Bloch
# ============================================================================

.PHONY: all clean

.SUFFIXES: .for

FC=gfortran -g3 -ggdb
#LIBS=-lpthread -llapack -lslatec
LIBS=-lpthread -llapack
FLAGS=
FLAGS+=$(INCS) $(LIBS)

SOURCES=ablochgen.for femutil.for gaussin.for matutil.for mtxutil.for preproc.for uelutil.for bloch.for

OBJECTS=$(SOURCES:.for=.o)

all: $(SOURCES) bloch

bloch: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS) $(FLAGS)


.for.o:
	$(FC) $(INCS) -c $< -o $@

clean:
	rm -f bloch *.mod *.o

echo:
	echo $(OBJS)
