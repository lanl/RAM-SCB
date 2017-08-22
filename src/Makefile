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
	IM_set_parameters.o
        
SCB_OBJECTS=                  \
        ModScbEquation.o      \
	ModScbEuler.o         \
	ModScbFunctions.o     \
	ModScbGrids.o         \
	ModScbInit.o          \
	ModScbInterp.o        \
	ModScbIO.o            \
	ModScbMain.o          \
	ModScbParams.o        \
	ModScbRun.o           \
	ModScbSpline.o        \
        ModScbVariables.o     \
        scb_ionoPot.o         \
        scb_physQuantities.o  \
        scb_testConvergence.o \
	Geopack-2008_dp.o\
	w2k_r8.o\
	w05.o

CPL_OBJECTS=            \
        ModRamCouple.o  \
        ModRamScb.o     \
        ModScbCouple.o

EXT_OBJECTS=       \
	nrmod.o    \
	nrtype.o   \
	nrutil.o   \
	recipes1.o \
	recipes2.o \
	mpi.o      \
	ModRamMpi.o 

SLATEC_OBJECTS = ../srcSlatec/*.o
${SLATEC_OBJECTS}:
	cd ../srcSlatec; make slatec


OBJECTS = ${RAM_OBJECTS} ${SCB_OBJECTS} ${CPL_OBJECTS} ${EXT_OBJECTS}

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
