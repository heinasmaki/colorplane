FCC = gfortran
FFL = -O3 -c

colorplane : colorplane.o colorplane_utils.o colorplane_consts.o
	$(FCC) -o colorplane colorplane.o colorplane_utils.o colorplane_consts.o

colorplane.o : colorplane.f90 colorplane_utils.o colorplane_consts.o
	$(FCC) $(FFL) colorplane.f90

colorplane_utils.o : colorplane_utils.f90 colorplane_consts.o
	$(FCC) $(FFL) colorplane_utils.f90

colorplane_consts.o: colorplane_consts.f90
	$(FCC) $(FFL) colorplane_consts.f90
