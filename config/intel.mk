OPTFLAGS=-O2
DBGFLAGS=-O0 -g -warn unused
FDBG:=-traceback -check all -ftrapuv
FFLAGS:=-r8 -heap-arrays
LDFLAGS:=-r8  $(LDFLAGS) -heap-arrays
MKLLIBS:= mkl_intel_lp64 mkl_sequential mkl_core

LAPACK:=-framework Intel_MKL
LAPACK:=-L$(MKLROOT)/lib/em64t $(addprefix -l,$(MKLLIBS))

