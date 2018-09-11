SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND

OBJECTS = \
	IM_wrapper.o

SEARCH_EXTRA = -I${COUPLERDIR} -I${LIBRARYDIR} -I../src
INCL_EXTRA   = -I${COUPLERDIR} -I${LIBRARYDIR} -I../src

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

MY_LIB = ${LIBDIR}/libIM.a

LIB: DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

RAM_LIB = ../src/libRAM_SCB.a

${MY_LIB}: ${RAM_LIB} ${OBJECTS} ${MAKEFILE_COMP_SELECT}
	cp -f ${RAM_LIB} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

clean:
	rm -f *.o

distclean: clean
	rm -f Makefile.DEPEND
