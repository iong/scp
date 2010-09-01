CC=gcc-mp-4.5
FC=gfortran-mp-4.5
OPTFLAGS=-O3
DBGFLAGS=-O0 -ggdb -Wall

CFLAGS=
FFLAGS=-fdefault-real-8 -fdefault-double-8 -ffree-line-length-0 -ffixed-line-length-0
FDBG:=-fbounds-check -Wimplicit

LAPACK=-lf77blas -lcblas -latlas

ifeq ($(OS),Darwin)
	LAPACK:=-framework vecLib
endif


