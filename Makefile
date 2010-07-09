#DBG=1
VPATH=$(dir $(shell readlink $(CURDIR)/$(firstword $(MAKEFILE_LIST))))
VGW=$(addsuffix .o,vgw rhss0 vgwrho vgw0)

ARCH=x86_64
OS=$(shell uname -s)

CC=mpicc
FC=mpif90

CPPFLAGS=$(shell pkg-config --cflags gsl)
LDFLAGS=$(shell pkg-config --libs-only-L gsl)
LIBS=-lgsl

include gcc.mk
#include intel.mk

ifdef DBG
	CFLAGS:=$(CFLAGS) $(CPPFLAGS) $(DBGFLAGS)
	FFLAGS:=$(FFLAGS) $(CPPFLAGS) $(DBGFLAGS) $(FDBG)
	LDFLAGS:=$(LDFLAGS) $(DBGFLAGS) $(FDBG)
else
	CFLAGS:=$(CFLAGS) $(CPPFLAGS) $(OPTFLAGS)
	FFLAGS:=$(FFLAGS) $(CPPFLAGS) $(OPTFLAGS)
	LDFLAGS:=$(LDFLAGS) $(OPTFLAGS)
endif

LIBS:= $(LIBS) $(LAPACK) -lm


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
