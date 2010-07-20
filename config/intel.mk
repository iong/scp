
OPTFLAGS=-fast
DBGFLAGS=-O0 -g
FDBG:=-traceback
FFLAGS:=-r8 -heap-arrays
LDFLAGS:=-nofor-main -r8  $(LDFLAGS) -heap-arrays -static-intel
MKLLIBS:= mkl_intel_lp64 mkl_sequential mkl_core

ifeq ($(OS),Linux)
	OPTFLAGS=-O3 -no-prec-div
	FDBG:=-check all -ftrapuv -traceback
endif

LAPACK:=$(addprefix $(MKLROOT)/lib/em64t/lib,$(addsuffix .a,$(MKLLIBS)))

