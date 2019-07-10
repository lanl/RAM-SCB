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
SEARCH_EXTRA = -I${NETCDF_PATH} -I${GSL_PATH}
INCL_EXTRA   = -I${NETCDF_PATH} -I${GSL_PATH}

# Keep RAM_SCB in single precision.
PRECISION = -frecord-marker=4 -fdefault-real-8 -fdefault-double-8 -fall-intrinsics -finit-local-zero

# List of fortran object files:
RAM_OBJECTS=                  \
	ModRamBoundary.o      \
	ModRamConst.o         \
	ModRamDrift.o         \
	ModRamEField.o        \
	ModRamFunctions.o     \
	ModRamGrids.o         \
	ModRamIndices.o       \
	ModRamInit.o          \
	ModRamIO.o            \
	ModRamLoss.o          \
	ModRamMain.o          \
	ModRamNCDF.o          \
	ModRamParams.o        \
	ModRamRestart.o       \
	ModRamRun.o           \
	ModRamSats.o          \
	ModRamTiming.o        \
        ModRamVariables.o     \
	ModRamWPI.o           \
        ModRamInjection.o     \
	IM_set_parameters.o

PLE_OBJECTS=                  \
        ModPlane.o            \
        ModRamPl_Ne.o         \
        plane_scb.o           \
        fct-source.o          \
        irisub_2012.o         \
        iritec.o              \
        iriflip.o             \
        iridreg.o             \
        irifun_2012.o         \
        nrlmsise00_sub.o      \
        igrf.o

SCB_OBJECTS=                  \
        ModScbEquation.o      \
	ModScbEuler.o         \
	ModScbFunctions.o     \
	ModScbGrids.o         \
	ModScbInit.o          \
	ModScbIO.o            \
	ModScbMain.o          \
	ModScbParams.o        \
	ModScbRun.o           \
        ModScbVariables.o     \
	ModScbCompute.o       \
	Geopack-2008_dp.o     \
        T89.o		      \
	T96.o		      \
	T02.o		      \
	TS04.o		      \
	TS07.o                \
	w2k_r8.o	      \
	w05.o

SCE_OBJECTS=                  \
	ModSceVariables.o     \
	ModSceGrids.o         \
	ModSceInit.o          \
        ModSceRun.o           \
	ModSceIono.o

GSL_OBJECTS=		\
	ModRamGSL.o	\
	RamGSL.o

CPL_OBJECTS=            \
        ModRamScbRun.o  \
        ModRamScb.o     \
        ModRamSce.o     \
	ModRamCouple.o 

EXT_OBJECTS=              \
        gaussian_filter.o \
	nrtype.o          \
	mpi.o             \
	ModRamMpi.o 

OBJECTS = ${RAM_OBJECTS} ${SCB_OBJECTS} ${SCE_OBJECTS} ${EXT_OBJECTS} ${GSL_OBJECTS} ${CPL_OBJECTS} ${PLE_OBJECTS}

# List of objects required for executable:
OBJECTS_EXE = Main.o

DEPEND:
	${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

help: 
	echo Help forthcoming.

# Rules for cleaning:
clean:
	rm -f *~ core *.o *.mod fort.* a.out *.exe *.a *.so *.protex

distclean:
	/bin/rm -f ../LibLocations.txt
	/bin/rm -f core *.o *.mod

# Build component library:
LIB:	DEPEND
	make ${MY_LIB}
	@echo ' '
	@echo ${MY_LIB} has been brought up to date.
	@echo ' '

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

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
