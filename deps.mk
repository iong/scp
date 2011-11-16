utils.o utils.mod: utils.f90
lsode.o lsode_mod.mod: lsode.f90
vgw.o vgw_mod.mod: vgw.f90 utils.mod lsode_mod.mod
vgwfm.o vgwfm_mod.mod: vgwfm.f90 vgw_mod.mod
dlsode.o: dlsode.f
rs.o: rs.f
xyz.o xyz.mod: xyz.f90
mainvars.o mainvars.mod: mainvars.f90 kubo.f90 utils.mod vgwfm_mod.mod
correlations.o correlations.mod: correlations.f90 vgwfm_mod.mod \
 mainvars.mod utils.mod utils.mod
main.o: main.f90 correlations.mod mainvars.mod utils.mod vgwfm_mod.mod \
 xyz.mod
