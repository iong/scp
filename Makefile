VPATH:=$(dir $(realpath Makefile))
ifeq "$(VPATH)" ""
    VPATH:=.
endif

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

GMD:=xyz spine gmdshort kubo correlations
MERGECVV= utils mergecvv

all: gmdshort

dbg: dbg.gmd

dbg.%:
	$(MAKE) DBG=1 $*

deps: deps.mk

include deps.mk

deps.mk:
	$(VPATH)/f90deps  $(addsuffix .f90,$(VGW) $(GMD) $(MERGECVV) $(LJ)) > $@


%.o : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f
	$(FC) $(FFLAGS) -c $<

gmdshort: $(addsuffix .o,$(VGW) $(GMD))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	$(RM) *.o *.mod

.PHONY: libpepc deps dbg
