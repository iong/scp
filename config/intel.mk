CC:=icc
FC:=ifort
OPTFLAGS=-fast
DBGFLAGS=-O0 -g -warn unused
FDBG:=-fpe0 #-traceback -check all -ftrapuv
FFLAGS:=
LDFLAGS:=$(LDFLAGS)
MKLLIBS:= mkl_intel_lp64 mkl_sequential mkl_core

LAPACK:=-framework Intel_MKL
LAPACK:=-L$(MKLROOT)/lib/em64t $(addprefix -l,$(MKLLIBS))

