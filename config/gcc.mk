CC=gcc-mp-4.4
FC=gfortran-mp-4.4
OPTFLAGS=-O3
DBGFLAGS=-O0 -ggdb -Wall

CFLAGS=
FFLAGS=-fdefault-real-8 -fdefault-double-8 -ffree-line-length-0 -ffixed-line-length-0 -fimplicit-none
FDBG:=-fbounds-check 

LAPACK=-lf77blas -lcblas -latlas

ifeq ($(OS),Darwin)
	LAPACK:=-framework vecLib
endif


