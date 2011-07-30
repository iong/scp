CC:=icc
FC:=ifort
OPTFLAGS=-fast
DBGFLAGS=-O0 -g -warn unused
FDBG:=-fpe0 -traceback -check all -ftrapuv
FFLAGS:=-heap-arrays
LDFLAGS:=$(LDFLAGS) -heap-arrays
MKLLIBS:= mkl_intel_lp64 mkl_sequential mkl_core

LAPACK:=-mkl=sequential

