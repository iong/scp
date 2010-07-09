OPTFLAGS=-O2 -g
DBGFLAGS=-O0 -ggdb -Wall

CFLAGS=-std=c99 
FFLAGS=-fdefault-real-8 -fdefault-double-8 -ffree-line-length-0 -ffixed-line-length-0
FDBG:=-fbounds-check -Wimplicit

ifeq ($(OS),Darwin)
	LAPACK:=-framework vecLib
endif


