CC=pgcc
FC=pgf95
OPTFLAGS=-fastsse -Minline
DBGFLAGS=-O0 -g
FDBG:=
FFLAGS:=-pc 64
LDFLAGS:=$(LDFLAGS) $(FFLAGS)
LAPACK=-lblas
