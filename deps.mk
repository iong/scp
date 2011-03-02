correlations.mod: correlations.o
correlations.o:  spine.mod utils.mod vgw.mod
gmd.o:  propagation.mod spine.mod utils.mod vgw.mod xyz.mod
gmdshort.o:  correlations.mod propagation.mod spine.mod utils.mod vgw.mod xyz.mod
kubo.o:  spine.mod utils.mod vgw.mod
lj.o:  utils.mod xyz.mod
load_defaults.o:  ljmc.mod mc.mod vgw.mod
mergecvv.o:  utils.mod
propagation.mod: propagation.o
rhs.mod: rhs.o
setup_ljmc.o:  ljmc.mod mc.mod utils.mod
spine.mod: spine.o
spine.o:  utils.mod vgw.mod
utils.mod: utils.o
vgw.mod: vgw.o
vgw.o:  propagation.mod potential_energy.f90 rhss0.f90 rhss1.f90 vgw0.f90 vgw1.f90
vgw0.o:  propagation.mod utils.mod
vgw1.o:  propagation.mod
xyz.mod: xyz.o
