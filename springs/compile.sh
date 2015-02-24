#!/bin/bash
# 
gfortran -Wall -o bloch.mod -c -g bloch.for
gfortran -Wall -o mat.mod -c -g mat.for
gfortran -Wall -o SPRING.o -g SPRING.FOR bloch.mod mat.mod -llapack
