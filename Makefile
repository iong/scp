VPATH:=$(dir $(realpath Makefile))
ifeq "$(VPATH)" ""
    VPATH:=.
endif

#DBG=1
COMPILER:=intel

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

GMD:=xyz spine gmd nose_hoover_chain
MERGECVV= utils mergecvv
LJ:=xyz utils  lj

all: gmd mergecvv

include deps.mk

deps.mk:
	$(VPATH)/f90deps  $(addsuffix .f90,$(VGW) $(GMD) $(MERGECVV) $(LJ)) > $@


%.o : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f
	$(FC) $(FFLAGS) -c $<

gmd: $(addsuffix .o,$(VGW) $(GMD))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

mergecvv: $(addsuffix .o, $(MERGECVV))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

lj: $(addsuffix .o, $(LJ))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

diag: diag.o rs.o utils.o
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	$(RM) *.o *.mod deps.mk

.PHONY: libpepc
