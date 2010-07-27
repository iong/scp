VPATH=$(dir $(shell readlink $(shell pwd)/$(firstword $(MAKEFILE_LIST))))

DBG=1
COMPILER:=intel

CC:=mpicc
FC:=mpif90
CPPFLAGS:=$(shell pkg-config --cflags gsl) -std=c99
LDFLAGS:=$(shell pkg-config --libs-only-L gsl)
LIBS:=-lgsl

OS=$(shell uname -s)
include config/$(COMPILER).mk

ifdef DBG
	CFLAGS:=$(CFLAGS) $(DBGFLAGS)
	FFLAGS:=$(FFLAGS) $(DBGFLAGS) $(FDBG)
	LDFLAGS:=$(LDFLAGS) $(DBGFLAGS) $(FDBG)
else
	CFLAGS:=$(CFLAGS) $(OPTFLAGS)
	FFLAGS:=$(FFLAGS) $(OPTFLAGS)
	LDFLAGS:=$(LDFLAGS) $(OPTFLAGS)
endif

LIBS:= $(LIBS) $(LAPACK) -lm

VGW=$(addsuffix .o,vgw rhss0 vgwrho vgw0)

all: ljmc

%.o : %.f90
	$(FC) $(FFLAGS) -c $^

%.o : %.f
	$(FC) $(FFLAGS) -c $^

ljmc: main.o dlsode.o vgwspb_H2_4G_Rc_Q_tau_SqrtMeff_Mar03.o vgw.o interaction_lists.o potential_energy.o rhss0.o vgw0.o
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	$(RM) *.o *.mod

.PHONY: libpepc
