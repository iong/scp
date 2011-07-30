CC=pgcc
FC=pgf95
OPTFLAGS=-fast
DBGFLAGS=-O0 -g -Ktrap=fp
FDBG:=-Mbounds  
FFLAGS:=
LDFLAGS:=$(LDFLAGS) $(FFLAGS)
LAPACK=-lblas
