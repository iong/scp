VPATH:=$(dir $(realpath Makefile))

#DBG=1
COMPILER:=pgi

OS=$(shell uname -s)
include $(VPATH)/config/$(COMPILER).mk

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

GMD:=xyz spine main

all: gmd

include $(VPATH)/deps.mk

deps.mk:
	./f90deps $(addsuffix .f90,$(VGW) $(GMD)) > $@


%.o : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f
	$(FC) $(FFLAGS) -c $<

gmd: $(addsuffix .o,$(VGW) $(GMD))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

diag: diag.o rs.o utils.o
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	$(RM) *.o *.mod deps.mk

.PHONY: libpepc
