
OPTFLAGS=-O3 -mdynamic-no-pic -no-prec-div
DBGFLAGS=-O0 -g
FDBG:=-check all -ftrapuv -traceback
FDBG:=-traceback
FFLAGS:=-r8 -heap-arrays
#LDFLAGS:=-nofor-main -r8  $(LDFLAGS) -heap-arrays -static-intel
LDFLAGS:=$(LDFLAGS) -static-intel

MKLLIBS:=mkl_sequential mkl_core
MKLLIBS:= mkl_intel_lp64 $(MKLLIBS)
ifeq ($(OS),Linux)
	OPTFLAGS=-O3 -no-prec-div
	MKLLIBS:= mkl_intel_lp64 $(MKLLIBS)
	#LDFLAGS:=-nofor-main  $(LDFLAGS) -heap-arrays -static-intel
	LDFLAGS:=$(LDFLAGS) -static-intel
endif
LAPACK:=$(addprefix $(MKLROOT)/lib/em64t/lib,$(addsuffix .a,$(MKLLIBS)))

