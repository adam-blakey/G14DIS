FC = gfortran
FFLAGS = -DARPACK -DMETIS_5_64 -DMUMPS -DDebug -DSPARSKIT -fbounds-check

mainsources = src/quadrature.f90

test:
	$(FC) $(FFLAGS) -o build/$@.out src/test.f90 $(mainsources) -J build

cleanall:
	(rm -f build/*.o)
	(rm -f build/*~)
	(rm -f build/*.mod)
	(rm -f build/*.out)