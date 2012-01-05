CC:=icc
FC:=ifort
OPTFLAGS=-fast
DBGFLAGS=-O0 -g -warn unused
FDBG:=-fpe0 -traceback -check all -check noarg_temp_created -ftrapuv
FFLAGS:=
MKLLIBS:= mkl_intel_lp64 mkl_sequential mkl_core

LAPACK:=-mkl=sequential

OS_RELEASE=$(shell uname -r | sed 's/\..*$$//')
ifeq ($(OS_RELEASE), 11)
	LDFLAGS += -Wl,-no_pie
endif
