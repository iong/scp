CC=pgcc
FC=pgf95
OPTFLAGS=-fastsse -Minline
DBGFLAGS=-O0 -g
FDBG:=
FFLAGS:=-r8 -mp
LDFLAGS:=$(LDFLAGS) $(FFLAGS)
LAPACK=-lblas
