VPATH:=$(dir $(realpath Makefile))
ifeq "$(VPATH)" ""
    VPATH:=.
endif

#DBG=1
COMPILER:=gcc

OS=$(shell uname -s)
include $(VPATH)/config/$(COMPILER).mk

ifdef DBG
	OPTFLAGS := $(DBGFLAGS)
	FFLAGS   += $(FDBG) 
	LDFLAGS  += $(FDBG)
endif

FFLAGS += $(OPTFLAGS)
CFLAGS += $(OPTFLAGS)
LDFLAGS += $(OPTFLAGS)


LIBS += $(LAPACK) -lm

VGW:=utils.f90 lsode.f90 vgw.f90 vgwfm.f90 dlsode.f rs.f
OH:=xyz.f90 mainvars.f90 kubo.f90 correlations.f90 main.f90
objects=$(addsuffix .o,$(basename $(1)))

all: OH

ifneq ($(wildcard $(VPATH)/deps.mk),)
include $(VPATH)/deps.mk
endif

deps:
	$(RM) deps.mk
	gfortran-mp-4.6 -MM -cpp $(FFLAGS) $(VGW) $(OH) > deps.mk
	#$(VPATH)/f90deps $(ALL_SRC) > deps.mk

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f
	$(FC) $(FFLAGS) -c $<

debug:
	@echo $(PWD)
	$(MAKE) -f $(THIS_MAKEFILE) DBG=1

OH: $(call objects,$(VGW) $(OH))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	$(RM) *.o *.mod

.PHONY: deps debug clean
