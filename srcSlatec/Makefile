# Basic makefile for building a reduced version of the Slatec library
# for use in Sorin Zaharia's 3D Equilibrium Magnetic Field code.
include ../Makefile.conf
include ../Makefile.def

# Keep RAM_SCB in single precision.
PRECISION = 

OBJECTS = \
	dpchic.o \
        dpchfe.o \
        dpchci.o \
        dpchcs.o \
        dpchce.o \
        dpchfd.o \
        dpchdf.o \
        dpchst.o \
        dpchsw.o \
        dchfdv.o \
	dlmach.o \
        d1mach.o \
        dqage.o  \
        dqagse.o \
        dqk15i.o \
        dqk15.o  \
        dqk21.o  \
        dqk31.o  \
        dqk41.o  \
        dqk51.o  \
        dqk61.o  \
        dqpsrt.o \
        dqagie.o \
        dqags.o  \
        dqelg.o  \
        dqag.o   \
        dqagi.o  \
        dqagie.o \
        dchfev.o \
        j4save.o \
        xerprn.o \
        xersve.o \
        xerhlt.o \
        xercnt.o \
        xermsg.o \
        xgetua.o \
        fdump.o  \
        i1mach.o

slatec: ${OBJECTS}

distclean:
	/bin/rm -f core *.o *.mod