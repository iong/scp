cmake_minimum_required(VERSION 2.8)

PROJECT(nr Fortran)

add_library(nr STATIC dbrent.f90 frprmn.f90 linmin.f90 mnbrak.f90 nr.f90
	nrtype.f90 nrutil.f90)
