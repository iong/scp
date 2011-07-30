#!/usr/bin/env python2.6

import os
import sys
import re
import tempfile

m=re.compile(r'^\s*$')

tempfile.tempdir = os.getcwd()

for fn in sys.argv[1:]:
    f=open(fn, "rw")
    outl = []
    for l in f.readlines():
        if not m.match(l):
            outl.append(l)
        else:
            break
    f.close()
    f=open(fn+".tmp", "w")
    f.writelines(outl)
    f.close()
    os.rename(fn+".tmp", fn)
            

