utils.mod: utils.o
propagation.mod: propagation.o
vgw.mod: vgw.o
xyz.mod: xyz.o
spine.mod: spine.o
spine.o: utils.mod vgw.mod
gmd.o: propagation.mod spine.mod utils.mod vgw.mod xyz.mod
utils.mod: utils.o
mergecvv.o: utils.mod
xyz.mod: xyz.o
utils.mod: utils.o
lj.o: utils.mod xyz.mod
