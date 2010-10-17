#!/bin/bash

cfg=pH2.in
for dt in 5d-4; do
    for rtol in 1d-4 1d-5 1d-6; do
        mv $cfg tmp.in
        sed -e "s/rtol.*$/rtol = $rtol/" -e "s/dt =.*$/dt = $dt/" tmp.in > $cfg
        echo dt=${dt} rtol=${rtol}
        time ./gmd r0_0000020000.xyz
        mv fort.31 dt_${dt}_rtol_${rtol}.dat
    done
done
