VPATH=$(dir $(shell readlink $(shell pwd)/$(firstword $(MAKEFILE_LIST))))

DBG=1
COMPILER:=intel

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

VGW:=vgw propagation rhs vgw0 vgw1 unpackg rhss0 rhss1\
       interaction_lists dlsode vgwspb_H2_4G_Rc_Q_tau_SqrtMeff_Mar03
VGW:=$(addsuffix .o,$(VGW))

LJMC:=setup_ljmc ljmc mc dump_Umin heat_capacity
LJMC:=$(addsuffix .o,$(LJMC))

all: ljmc

%.o : %.f90
	$(FC) $(FFLAGS) -c $^

%.o : %.f
	$(FC) $(FFLAGS) -c $^

ljmc: $(VGW) $(LJMC)
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	$(RM) *.o *.mod

.PHONY: libpepc
