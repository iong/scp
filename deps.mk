utils.mod: utils.o
propagation.mod: propagation.o
vgw.mod: vgw.o
xyz.mod: xyz.o
spine.mod: spine.o
spine.o: utils.mod vgw.mod
main.o: propagation.mod spine.mod utils.mod vgw.mod xyz.mod
xyz.mod: xyz.o
mergecvv.o: utils.mod vgw.mod xyz.mod
