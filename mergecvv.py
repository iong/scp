#!/usr/bin/env python2.6

import sys
import glob
import os.path
from numpy import *
import scipy.io

def read_header(fname):
    f=open(fname, "r")
    labels=[l.strip() for l in f.readlines()]
    f.close()
    return labels


def merge_cvv(d):
    files=glob.glob(d+"/dump/*_Cvv.dat")

    Cxxall=array([loadtxt(zf) for zf in files])
    Cxx=mean(Cxxall, 0)
    chunk_size = size(Cxxall, 0) / 10
    Cxx_=zeros((10,)+shape(Cxx))
    print shape(Cxx_), shape(Cxxall)
    for i in range(0,10):
        Cxx_[i,:,:] = mean(Cxxall[i*chunk_size:(i+1)*chunk_size,:,:], 0)
    Cxxstd = std(Cxx_, 0)

    header=read_header(d+"/dump/map")
    mdict=dict([(header[i], Cxx[:,i]) for i in range(0,len(header))])
    mdict.update([(header[i]+"std", Cxxstd[:,i]) for i in range(1,len(header))])
    scipy.io.savemat(d+'/Cxx.mat', mdict)

    mdict={'Cxxall':Cxxall}
    scipy.io.savemat(d+'/Cxxall.mat', mdict)

    print shape(Cxxstd[:,1])
    for i in range(1,len(header)):
        savetxt(d+'/'+header[i]+".dat", vstack((Cxx[:,(0, i)].T, Cxxstd[:,i])).T)

d="/Users/ionut/data/Cvvshort"
if len(sys.argv)>1:
    d=sys.argv[1]
    
merge_cvv(d)
