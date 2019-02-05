# C language related part of Makefile.conf: Makefile.cc
C_COMPILER_NAME=cc

COMPILE.c     = cc
COMPILE.mpicc = mpicc
COMPILE.mpicxx= mpicxx

CPPLIB = -lstdc++ -lmpi_cxx

#DEBUGC = -g

.SUFFIXES: .c .cpp

FLAGC = ${SEARCH_C} ${FLAGC_EXTRA} -c ${OPT3} ${OPENMPFLAG} ${DEBUGC}

.c.o:
	${COMPILE.c} ${FLAGC} $< -o $@

.cpp.o:
	${COMPILE.mpicxx} ${FLAGC} $< -o $@
