CC=pgcc
FC=pgf95
OPTFLAGS=-fastsse -Minline
DBGFLAGS=-O0 -g -Mbounds
FDBG:=
FFLAGS:=
LDFLAGS:=$(LDFLAGS) $(FFLAGS)
LAPACK=-lblas
