CC=gcc
FC=gfortran
OPTFLAGS=-O3
DBGFLAGS=-O0 -ggdb -Wall

CFLAGS=
FFLAGS=-ffree-line-length-0 -ffixed-line-length-0 -fimplicit-none
FDBG:=-fbounds-check -ffpe-trap=invalid,zero,overflow,denormal

LAPACK=-lf77blas -lcblas -latlas

ifeq ($(OS),Darwin)
	LAPACK:=-framework vecLib
endif


