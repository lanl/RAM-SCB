# RAM-SCB makefile

# Include statements:
include ../Makefile.conf
include ../Makefile.def

include Makefile.DEPEND
include Makefile.RULES

# Name the executable and library target.
EXE = ${BINDIR}/ram_scb.exe
MY_LIB = libRAM_SCB.a

# Other required libraries.
LIBSHARE = ${LIBDIR}/libSHARE.a
LIBNOMPI  = ${LIBDIR}/libNOMPI.a

# Extra search directories for -I option (-p for Absoft!!!)
# Set twice to work with old and new SWMF versions.
SEARCH_EXTRA = -I${PSPLINE_PATH} -I${NETCDF_PATH}
INCL_EXTRA   = -I${PSPLINE_PATH} -I${NETCDF_PATH}

# Keep RAM_SCB in single precision.
PRECISION = 

# List of fortran object files:
RAM_OBJECTS=                \
	ModRamMain.o	    \
	mpi.o               \
	ModRamMpi.o         \
	nrtype.o            \
	nrutil.o            \
	Module1.o           \
	ModRamIO.o          \
	ModRamCouple.o      \
	ModRamFunctions.o   \
	ModRamTiming.o      \
	ModRamIndices.o     \
	ModRamGeomlt.o	    \
	ModPlane.o          \
	ModRamPl_Ne.o       \
	nr.o                \
	adi_fast.o          \
	erf.o               \
	fct-source.o        \
	irisub_2012.o       \
	iritec.o            \
	iri87.o             \
	igrf.o              \
	iriflip.o           \
	iridreg.o           \
	irifun_2012.o       \
	Module_RAM.o        \
	ModRamSats.o        \
	msis86.o            \
	nrlmsise00_sub.o    \
	recipes2.o          \
	IM_set_parameters.o \
	init_ramscb.o       \
	ram_update_electric.o \
	run_ramscb.o \
	finalize_ramscb.o   \
	TM03_DensPres.o

SCB_OBJECTS= \
	ram_all.o \
	alfges.o\
	anisEq.o\
	barycoord_mine_2D.o\
	baryth.o\
	bounextp.o\
	closed_field.o\
	Computational_domain.o\
	computeBandJacob_initial.o\
	computeEquilibrium.o\
	covsrt.o\
	cubic.o\
	directAlpha.o\
	directPsi.o\
	dmat_solve.o\
	dps_general.o\
	d_swap.o\
	energy.o\
	entropy.o\
	extap.o\
	fdriv.o\
	findR.o\
	findx.o\
	gammln.o\
	gammp.o\
	gaussj.o\
	gcf.o\
	Geopack-2008_dp.o\
	gser.o\
	handle_err.o\
	hRAM.o\
	input.o\
	interpolatepsi.o\
	interpolatepsiR.o\
	Interpolation_delaunay_2D.o\
	iterateAlpha.o\
	iteratePsi.o\
	locate.o\
	latitudes.o\
	mapalpha.o\
	mappsi.o\
	maptheta.o\
	mapthe1d.o\
	mapthe1d_down.o\
	matrixInversion.o\
	metrica.o\
	metric.o\
	metrics.o\
	module_Gauss.o\
	Module_inversion.o\
	Module_points.o\
	module_test.o\
	my_flush.o\
	netcdf_density_dump.o\
	netcdf_flux_dump.o\
	netcdf_flux_dump_3D.o\
	newj.o\
	newk.o\
	newlin2.o\
	PARALLEL_3D_EQ.o\
	plane_scb.o\
	plot_physical.o\
	polint.o\
	pressure.o\
	pressureRad.o\
	psiFunctions.o\
	psiges.o\
	qromb.o\
	qtrap.o\
	quad2.o\
	quad.o\
	secant.o\
	Spline_2D_array.o\
	Spline_2D_derivs.o\
	Spline_2D_point.o\
	Spline_2D_periodic.o\
	Spline_coord_2nd_derivs.o\
	Spline_coord_derivs.o\
	Spline_coord_derivs_large.o\
	Spline_derivs_1D.o\
	Spline_derivs_1D_zeta.o\
	Spline_1D_periodic.o\
	spline.o\
	splint.o\
	svbksb.o\
	test_Convergence_anisotropic.o\
	trapzd.o\
	tridag.o\
	Write_BEandJacob.o\
	Write_pressure_anisotropic.o\
	Write_pressure_isotropic.o\
	pRoeRad.o\
	Interpolation_natgrid_2D.o\
	Interpolation_natgrid_ShellBuild.o\
	SavGol.o\
	SavGol7.o\
	SavGol1D.o\
	ionospheric_potential.o\
	Euler_derivatives.o \
	Interpolation_delaunay_Euler_2D.o \
	Write_ionospheric_potential.o \
	w2k_r8.o\
	w05.o

SLATEC_OBJECTS = ../srcSlatec/*.o
${SLATEC_OBJECTS}:
	cd ../srcSlatec; make slatec


OBJECTS = ${RAM_OBJECTS} ${SCB_OBJECTS}

# List of objects required for executable:
OBJECTS_EXE = Main.o

DEPEND:
	@${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

help: 
	echo Help forthcoming.

# Rules for cleaning: 
distclean:
	/bin/rm -f ../LibLocations.txt
	/bin/rm -f core *.o *.mod

# Build component library:
LIB:	DEPEND
	make ${MY_LIB}
	@echo ' '
	@echo ${MY_LIB} has been brought up to date.
	@echo ' '

${MY_LIB}: ${SLATEC_OBJECTS} ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS} ${SLATEC_OBJECTS}

# Rule for creating executable:
RAM_SCB: DEPEND
	@make ${EXE}
	@echo ' '
	@echo "Program ${EXE} has been brought up to date."
	@echo ' '

# Libraries should be compiled first, because modules are used in main
${EXE}: ${OBJECTS_EXE} ${MY_LIB} ${LIBSHARE}
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ../${MY_LIB}; \
		ar -x ${LIBSHARE};
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} Tmp_/*.o ${LBLAS} ${Lflag1}
	rm -rf Tmp_
