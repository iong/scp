VPATH:=$(dir $(shell readlink $(shell pwd)/Makefile))

#DBG=1
COMPILER:=pgi

OS=$(shell uname -s)
include $(VPATH)/config/$(COMPILER).mk
CC:=mpicc
FC:=mpif90

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

VGW:=utils propagation vgw unpackg\
       interaction_lists dlsode vgwspb_H2_4G_Rc_Q_tau_SqrtMeff_Mar03
VGW:=$(addsuffix .o,$(VGW))

LJMC:=xyz ljmc_mod mc setup_ljmc ljmc heat_capacity load_defaults populate_cube
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
