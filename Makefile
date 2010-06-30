DBG=1
VPATH=$(dir $(shell readlink $(CURDIR)/$(firstword $(MAKEFILE_LIST))))
CC=icc
#FC=gfortran-mp-4.4
FC=ifort
VGW=$(addsuffix .o,vgw rhss0 vgwrho vgw0)

CPPFLAGS=-I/opt/local/include
OPTFLAGS=-O2 -g
DBGFLAGS=-O0 -ggdb -Wall

CFLAGS=-std=c99 
FFLAGS=-fdefault-real-8 -fdefault-double-8 -ffree-line-length-0 -ffixed-line-length-0
FDBG:=-fbounds-check -Wimplicit

LDFLAGS=-L/opt/local/lib
LIBS=-lgsl
LAPACK:=-framework vecLib

ARCH=x86_64
OS=$(shell uname -s)

ifeq ($(FC),ifort)
	OPTFLAGS=-O3 -mdynamic-no-pic -no-prec-div
	DBGFLAGS=-O0 -g
	FDBG:=-check all -ftrapuv -traceback
	FDBG:=-traceback
	FFLAGS:=-r8 -heap-arrays
	LDFLAGS:=-nofor-main -r8  $(LDFLAGS) -heap-arrays -static-intel

	MKLLIBS:=mkl_sequential mkl_core
	MKLLIBS:= mkl_intel_lp64 $(MKLLIBS)
ifeq ($(OS),Linux)
	OPTFLAGS=-O3 -no-prec-div
	MKLLIBS:= mkl_intel_lp64 $(MKLLIBS)
	LDFLAGS:=-nofor-main  $(LDFLAGS) -heap-arrays -static-intel
endif
	LAPACK:=$(addprefix $(MKLROOT)/lib/em64t/lib,$(addsuffix .a,$(MKLLIBS)))
endif

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

#ljmc: main.o dlsode.o vgwspb_H2_4G_Rc_Q_tau_SqrtMeff_Mar03.o
ljmc: f77i.o dlsode.o $(VGW) ljmc.o vgwspb_H2_4G_Rc_Q_tau_SqrtMeff_Mar03.o
	$(FC)  $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	$(RM) *.o *.mod

.PHONY: libpepc
