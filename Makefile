# Top level makefile for IM/RAM_SCB
default : RAM_SCB

include Makefile.def

INSTALLFILES = src/Makefile.DEPEND \
	       src/Makefile.RULES  \
	       srcInterface/Makefile.DEPEND 

help:
	@echo ' '
	@echo ' You can "make" the following:'
	@echo ' '
	@echo '    <default> Build the RAM-SCB executable.'
	@echo ' '
	@echo '    help         (makefile option list)'
	@echo '    install      (install RAM_SCB)'
	@echo '    PDF          (Make PDF version of the documentation)'
	@echo '    HTML         (Make HTML version of the documentation)'
	@echo '    test         (run all tests for RAM-SCB)'
	@echo '    test_help    (show all options for running the tests)'
	@echo '    LIB		(Component library libIM for SWMF)'

PDF:
	@cd doc/Tex; make PDF

RAM_SCB:
	@cd ${SHAREDIR}; 	make LIB
	@cd src;   		make LIB
	@cd src;   		make RAM_SCB

install:
	@touch ${INSTALLFILES}


LIB:
	cd src; make LIB
	cd srcInterface; make LIB

clean:
	@touch ${INSTALLFILES}
	@cd src;          make clean
	@cd srcSlatec;    make clean
	@cd srcInterface; make clean
	@(if [ -d util ];  then cd util;  make clean; fi);
	@(if [ -d share ]; then cd share; make clean; fi);

distclean:
	./Config.pl -uninstall

allclean:
	@touch ${INSTALLFILES}
	@(if [ -d srcNetcdf ];  then rm -rf srcNetcdf;  fi);
	@(if [ -d srcPspline ]; then rm -rf srcPspline; fi);
	@cd src; make distclean
	@cd srcInterface; make distclean
	@cd srcSlatec; make distclean
	rm -f *~

rundir: 
	mkdir -p ${RUNDIR}/IM/output
	cp input/RamIndices.txt ${RUNDIR}/
#	cp input/apf107.dat ${RUNDIR}/
#	cp input/ig_rz.dat ${RUNDIR}/
	cp input/newtau.dat ${RUNDIR}/
	cp input/ne_full.dat ${RUNDIR}/
#	cp input/dgrf*.dat ${RUNDIR}/
#	cp input/igrf*.dat ${RUNDIR}/
#	cp input/ccir*.asc ${RUNDIR}/
#	cp input/ursi*.asc ${RUNDIR}/
	cd ${RUNDIR}; \
	ln -s ../input/bav_diffcoef_chorus_rpa_Kp*.PAonly.dat .
	cd ${RUNDIR}/IM/output; \
		mkdir -p ram/Dsbnd scb/Day00
	cd ${RUNDIR}/IM; \
		mkdir input_ram input_scb output_swmf;    \
		mkdir restartIN restartOUT;               \
		tar xzf ${IMDIR}/input/ramscb_inputs.tgz; \
		mv Input_git/w2k.dat ../;			  \
		mv Input_git/W05_coeff.dat ../;			\
		mv Input_git/omni.txt ../;			  \
		mv Input_git/f2ini* Input_git/*geomlt*.txt input_ram/;	  \
		cp Input_git/hI_output_0000.dat output/scb/hI_output_d20130317_t000000.dat;	  \
		mv Input_git/hI_output_0000.dat output/scb/;         \
		mv Input_git/hI_dipole.dat output/scb/;	          \
		mv Input_git/dipole_config.cdf Input_git/t89*.cdf input_scb/;
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR} ; \
		cp ${IMDIR}/Param/PARAM.in.default ./PARAM.in; \
		ln -s ${BINDIR}/ram_scb.exe .; \
		ln -s IM/* .; \
		rm -f output; \
		ln -s IM/output/ram output_ram; \
		ln -s IM/output/scb output_scb; \
	fi)


#---------
# TESTS
#---------
TESTDIR = run_test

test:
	-@(make test1)
	-@(make test2)
	-@(make test3)

test_help:
	@echo "Preceed all commands with 'make'..."
	@echo "   test				Run all tests in serial."
	@echo "   test1                        Test RAM only (dipole+VS, no SCB.)"
	@echo "   test2                        Test1 + Restart to test restart files."
	@echo "   test3                        Test RAM-SCB using W2k, T89, and electrons."
	@echo "   test1 MPIRUN='mpiexec -n 2'	Run test1 in parallel (on two cpus.)"
	@echo "Each test has commands to perform individual steps, e.g.:"
	@echo "   test1_compile		Compile code for test1"
	@echo "   test1_rundir			Create run directory for test."
	@echo "   test1_run			Execute the test simulation."
	@echo "   test1_check			Compare results to reference solution."


#TEST 1----------------------------------
test1:
	@echo "starting..." > test1.diff
	@echo "test1_compile..." >> test1.diff
	make test1_compile
	@echo "test1_rundir..." >> test1.diff
	make test1_rundir PARAMFILE=PARAM.in.test1
	@echo "test1_run..." >> test1.diff
	make test1_run MPIRUN=
	@echo "test1_check..." >> test1.diff
	make test1_check

test1_compile:
	make

test1_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE="YES"
	cp Param/${PARAMFILE} ${TESTDIR}/PARAM.in
	cp input/sat*.dat ${TESTDIR}/

test1_run:
	cd ${TESTDIR}; ${MPIRUN} ./ram_scb.exe > runlog

test1_check:
	${SCRIPTDIR}/DiffNum.pl -b 			 \
		${TESTDIR}/output_ram/Dsbnd/ds_0108_h.dat \
		${IMDIR}/output/test1/dsbnd.ref          \
		> test1.diff
	${SCRIPTDIR}/DiffNum.pl -b 		        \
		${TESTDIR}/output_ram/pressure_0001.in  \
		${IMDIR}/output/test1/pressure.ref      \
		>> test1.diff
	${SCRIPTDIR}/DiffNum.pl -b -a=0.00001		     \
		${TESTDIR}/output_ram/efield_000.in  \
		${IMDIR}/output/test1/weq01.ref      \
		>> test1.diff
	${SCRIPTDIR}/DiffNum.pl -b 		      \
		${TESTDIR}/output_ram/log_n000000.log \
		${IMDIR}/output/test1/log.ref         \
		>> test1.diff
	#${SCRIPTDIR}/DiffNum.pl -b 	        \
	#	${TESTDIR}/output_ram/sat1.nc   \
	#	${IMDIR}/output/test1/sat1.ref  \
	#	>> test1.diff
	#${SCRIPTDIR}/DiffNum.pl -b 	        \
	#	${TESTDIR}/output_ram/sat2.nc   \
	#	${IMDIR}/output/test1/sat2.ref  \
	#	>> test1.diff
	${SCRIPTDIR}/DiffNum.pl -b		  \
		${TESTDIR}/output_ram/ram000_o.l  \
		${IMDIR}/output/test1/ram_o.l.ref \
		>> test1.diff	
	${SCRIPTDIR}/DiffNum.pl -b		  \
		${TESTDIR}/output_ram/ram000_o.t  \
		${IMDIR}/output/test1/ram_o.t.ref \
		>> test1.diff	
	@echo "Test Successful!"

#TEST 2----------------------------------
test2:
	@echo "starting..." > test2.diff
	@echo "test2_compile..." >> test2.diff
	make test2_compile
	@echo "test2_rundir..." >> test2.diff
	make test2_rundir PARAMFILE=PARAM.in.test2
	@echo "test2_run..." >> test2.diff
	make test2_run MPIRUN=
	@echo "test2_check..." >> test2.diff
	make test2_check

test2_compile:
	make

test2_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE="YES"
	cp Param/${PARAMFILE}.* ${TESTDIR}/

test2_run:
	cd ${TESTDIR}; \
	rm PARAM.in; ln -s PARAM.in.test2.1st PARAM.in; \
	${MPIRUN} ./ram_scb.exe > runlog1; \
	mv IM/restartOUT/* IM/restartIN/;  \
	rm PARAM.in; ln -s PARAM.in.test2.2nd PARAM.in; \
	${MPIRUN} ./ram_scb.exe > runlog2;	

test2_check:
	${SCRIPTDIR}/DiffNum.pl -b -a=0.1	        \
		${TESTDIR}/output_ram/pressure_0001.in  \
		${IMDIR}/output/test1/pressure.ref      \
		> test2.diff
	${IMDIR}/Scripts/CatLog.py ${TESTDIR}/output_ram/log_n*.log
	${SCRIPTDIR}/DiffNum.pl -b -a=0.01	      \
		${TESTDIR}/output_ram/log_n000000.log \
		${IMDIR}/output/test1/log.ref         \
		>> test2.diff
	@echo "Test Successful!"

#TEST 3----------------------------------
test3:
	@echo "starting..." > test3.diff
	@echo "test3_compile..." >> test3.diff
	make test3_compile
	@echo "test3_rundir..." >> test3.diff
	make test3_rundir PARAMFILE=PARAM.in.test3
	@echo "test3_run..." >> test3.diff
	make test3_run MPIRUN=
#	@echo "test3_check..." >> test3.diff
#	make test3_check

test3_compile:
	make

test3_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE="YES"
	cp Param/${PARAMFILE} ${TESTDIR}/PARAM.in

test3_run:
	cd ${TESTDIR}; ${MPIRUN} ./ram_scb.exe > runlog;

test3_check:
	${SCRIPTDIR}/DiffNum.pl -b 		        	    \
		${TESTDIR}/output_ram/pressure_d20050831_t091000.in \
		${IMDIR}/output/test3/pressure.ref      	    \
		> test3.diff
	${SCRIPTDIR}/DiffNum.pl -b		  		  \
		${TESTDIR}/output_ram/efield_d20050831_t090000.in \
		${IMDIR}/output/test3/efield.ref 		  \
		>> test3.diff
	${SCRIPTDIR}/DiffNum.pl -b 		      \
		${TESTDIR}/output_ram/log_n000000.log \
		${IMDIR}/output/test3/log.ref         \
		>> test3.diff
	${SCRIPTDIR}/DiffNum.pl -b 				      \
		${TESTDIR}/output_scb/hI_output_d20050831_t091000.dat \
		${IMDIR}/output/test3/hI.ref 		 	      \
		>> test3.diff
	${SCRIPTDIR}/DiffNum.pl -b 		   \
		${TESTDIR}/output_scb/currents.nc  \
		${IMDIR}/output/test3/currents.ref \
		>> test3.diff
	@echo "Test Successful!"
