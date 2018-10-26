#!/bin/bash
# 
gfortran -Wall -o bloch.mod -c -g bloch.for
gfortran -Wall -o mat.mod -c -g mat.for
gfortran -Wall -o spring.o -g spring.for bloch.mod mat.mod -llapack
rm *.mod
