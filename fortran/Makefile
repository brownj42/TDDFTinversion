FFORTRAN=ifort
LLINKER=ifort
FFFLAGS=-O2 -xHost  -fopenmp  -I${MKLROOT}/include  
LLFLAGS= -L. -larpack  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
FORTRAN=gfortran
LINKER=gfortran
FFLAGS= -O2 -fopenmp  -fdefault-real-8 -fdefault-double-8 
LFLAGS= -L/usr/lib64 -larpack -lopenblas_openmp -llapack -fopenmp  -lm 
OBJECTFILES1= nchoosekmod.o mapmod.o matmul_mod.o minresqlpDataModule.o lsqrDataModule.o lsqrblasInterface.o lsqrModule.o minresqlpBlasModule.o minresqlpModule.o derivedtypes.o secondquant_mod.o secondquant_wrap.o keomod.o readin_data_mod.o potmod.o sortmod.o pot_invert_mod.o density_mod.o density_wrap.o assign_phases_mod.o hamiltonian_mod.o hamiltonian_wrap.o initial_states_mod.o outputdata_mod.o KSadvance_mod.o tddftinversion.o


default: tddftinversion 

tddftinversion: $(OBJECTFILES1) $(OBJECTFILESF)
	$(LINKER) $(OBJECTFILES1) $(OBJECTFILESF) -o  $@ $(LFLAGS)

$(OBJECTFILES1):
	$(FORTRAN) $(FFLAGS) $*.f90 -c

secondquant_wrap.o:secondquant_wrap.f90
nchoosekmod.o:nchoosekmod.f90
secondquant_mod.o:secondquant_mod.f90
mapmod.o:mapmod.f90
hamiltonian_wrap.o:hamiltonian_wrap.f90
readin_data_mod.o:readin_data_mod.f90
KSadvance_mod.o:KSadvance_mod.f90
derivedtypes.o:derivedtypes.f90
keomod.o:keomod.f90
potmod.o:potmod.f90
outputdata_mod.o:outputdata_mod.f90
matmul_mod.o:matmul_mod.f90
density_mod.o:density_mod.f90
density_wrap.o:density_wrap.f90
initial_states_mod.o:initial_states_mod.f90
tddftinversion.o:tddftinversion.f90
hamiltonian_mod.o:hamiltonian_mod.f90
sortmod.o:sortmod.f90
assign_phases_mod.o:assign_phases_mod.f90
pot_invert_mod.o:pot_invert_mod.f90
minresqlpDataModule.o:minresqlpDataModule.f90
minresqlpBlasModule.o:minresqlpBlasModule.f90
minresqlpModule.o:minresqlpModule.f90
lsqrDataModule.o:lsqrDataModule.f90
lsqrblasInteface.o:lsqrblasInterface.f90
lsqrModule.o:lsqrModule.f90


clean:
	rm *.o
	rm *.mod
