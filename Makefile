# Top level makefile for IM/RAM_SCB
default : RAM_SCB

include Makefile.def

srcDir = src
INSTALLFILES = ${srcDir}/Makefile.DEPEND \
	       ${srcDir}/Makefile.RULES  \
	       srcInterface/Makefile.DEPEND \
	       srcExternal/Makefile.DEPEND

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
	@cd ${SHAREDIR}; make LIB
	@cd srcExternal; make LIB
	@cd ${srcDir};   make LIB
	@cd ${srcDir};   make RAM_SCB

install:
	@touch ${INSTALLFILES}


LIB:
	cd ${srcDir};    make LIB
	cd srcInterface; make LIB
	cd srcExternal;  make LIB

clean:
	@touch ${INSTALLFILES}
	@cd ${srcDir};          make clean
	@cd srcInterface; make clean
	@cd srcExternal; make clean
	@(if [ -d util ];  then cd util;  make clean; fi);
	@(if [ -d share ]; then cd share; make clean; fi);

distclean:
	./Config.pl -uninstall

allclean:
	@touch ${INSTALLFILES}
	@(if [ -d srcNetcdf ];  then rm -rf srcNetcdf;  fi);
	@(if [ -d srcPspline ]; then rm -rf srcPspline; fi);
	@cd ${srcDir}; make distclean
	@cd srcInterface; make distclean
	rm -f *~

rundir: 
	mkdir -p ${RUNDIR}/IM/output
	cp input/RamIndices.txt ${RUNDIR}/
	cp input/apf107.dat ${RUNDIR}/
	cp input/ig_rz.dat ${RUNDIR}/
	cp input/newtau.dat ${RUNDIR}/
	cp input/ne_full.dat ${RUNDIR}/
	cp input/initialization.nc ${RUNDIR}/IM/
	cp input/QinDenton_20130317_1min.txt ${RUNDIR}/IM/
	cp input/NitrogenCrossSections.dat ${RUNDIR}/IM/
	cp input/AEindex.txt ${RUNDIR}/
	cd ${RUNDIR}; \
	ln -s ../input/bav_diffcoef_chorus_rpa_Kp*.PAonly.dat .
	cd ${RUNDIR}/IM/output; \
		mkdir -p ram/Dsbnd scb/Day00
	cd ${RUNDIR}/IM; \
		mkdir input_ram input_scb output_swmf;    \
		mkdir restartIN restartOUT;               \
		tar xzf ${IMDIR}/input/ramscb_inputs.tgz; \
		mv Input_git/EMIC_model input_ram/; \
		mv Input_git/HBand input_ram/;  \
		mv Input_git/HeBand input_ram/; \
		mv Input_git/w2k.dat ../;		  \
		mv Input_git/W05_coeff.dat ../;		    \
		mv Input_git/omni.txt ../;		    \
		mv Input_git/*geomlt*.txt input_ram/; \
                mv initialization.nc input_ram/;            \
		mv QinDenton_20130317_1min.txt input_scb/;  \
		mv NitrogenCrossSections.dat input_ram/;
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
TESTDIR1 = run_test1
TESTDIR2 = run_test2
TESTDIR3 = run_test3
TESTDIR4 = run_test4
TESTDIRC = run_test

test:
	@(make test1)
	@(make test2)
	@(make test3)
	@(make test4)

testTravis:
	@(make test3)

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


#TEST TRAVIS
testTravisLong:
	@echo "starting..." > test1.diff
	@echo "test1_compile..." >> test1.diff
	make test1_compile
	@echo "test1_rundir..." >> test1.diff
	make test1_rundir PARAMFILE=PARAM.in.test1
	@echo "test1_run..." >> test1.diff
	make test1_run MPIRUN=
	@echo "test1_check..." >> test1.diff
	make test1_check
	@echo "starting..." > test3.diff
	@echo "test3_compile..." >> test3.diff
	make test3_compile
	@echo "test3_rundir..." >> test3.diff
	make test3_rundir PARAMFILE=PARAM.in.test3
	@echo "test3_run..." >> test3.diff
	make test3_run MPIRUN=
	@echo "starting..." > test4.diff
	@echo "test4_compile..." >> test4.diff
	make test4_compile
	@echo "test4_rundir..." >> test4.diff
	make test4_rundir PARAMFILE=PARAM.in.test4
	@echo "test4_run..." >> test4.diff
	make test4_run MPIRUN=
	@echo "test4_check..." >> test4.diff
	make testTravis_check

#TRAVIS Test
testTravis_check:
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9 \
                ${TESTDIR3}/output_ram/pressure_d20130317_t001500.dat \
                ${TESTDIR4}/output_ram/pressure_d20130317_t001500.dat \
                > testTravis.diff
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9                             \
                ${TESTDIR3}/output_scb/hI_output_d20130317_t001500.dat \
                ${TESTDIR4}/output_scb/hI_output_d20130317_t001500.dat \
                >> testTravis.diff
	ncdump -v "Flux_H","B_xyz"                                     \
               ${TESTDIR3}/output_ram/sat1_d20130317_t000000.nc        \
               | sed -e '1,/data:/d' >                                 \
               ${TESTDIR3}/output_ram/sat1.test
	ncrcat ${TESTDIR4}/output_ram/sat1_d20130317_t000000.nc       \
               ${TESTDIR4}/output_ram/sat1_d20130317_t001000.nc       \
               ${TESTDIR4}/output_ram/sat1.nc
	ncdump -v "Flux_H","B_xyz" ${TESTDIR4}/output_ram/sat1.nc     \
               | sed -e '1,/data:/d' >                                \
               ${TESTDIR4}/output_ram/sat1.test        
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9                            \
                ${TESTDIR3}/output_ram/sat1.test                      \
                ${TESTDIR4}/output_ram/sat1.test                      \
                >> testTravis.diff
	@echo "Test Successful!"

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
	rm -rf ${TESTDIR1}
	make rundir RUNDIR=${TESTDIR1} STANDALONE="YES"
	cp Param/${PARAMFILE} ${TESTDIR1}/PARAM.in
	cp input/sat*.dat ${TESTDIR1}/

test1_run:
	cd ${TESTDIR1}; ${MPIRUN} ./ram_scb.exe | tee runlog

test1_check:
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9				\
		${TESTDIR1}/output_ram/log_d20130317_t000000.log	\
		${IMDIR}/output/test1/log.ref				\
		> test1.diff
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9	                        \
		${TESTDIR1}/output_ram/pressure_d20130317_t001500.dat   \
		${IMDIR}/output/test1/pressure.ref                      \
		>> test1.diff			        
	ncdump -v "Flux_H","B_xyz"                              	\
               ${TESTDIR1}/output_ram/sat1_d20130317_t000000.nc 	\
               | sed -e '1,/data:/d' >                          	\
               ${TESTDIR1}/output_ram/sat1.test
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9 	                        \
		${TESTDIR1}/output_ram/sat1.test                	\
		${IMDIR}/output/test1/sat1.ref                  	\
		>> test1.diff
	ncdump -v "Flux_H","B_xyz"                              	\
               ${TESTDIR1}/output_ram/sat2_d20130317_t000000.nc 	\
               | sed -e '1,/data:/d' >                          	\
               ${TESTDIR1}/output_ram/sat2.test
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9 	                        \
		${TESTDIR1}/output_ram/sat2.test                	\
		${IMDIR}/output/test1/sat2.ref                  	\
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
	rm -rf ${TESTDIR2}
	make rundir RUNDIR=${TESTDIR2} STANDALONE="YES"
	cp Param/${PARAMFILE}.* ${TESTDIR2}/
	cp input/sat*.dat ${TESTDIR2}/

test2_run:
	cd ${TESTDIR2};                                 \
	rm PARAM.in; ln -s PARAM.in.test2.1st PARAM.in; \
	${MPIRUN} ./ram_scb.exe | tee runlog1;              \
	rm PARAM.in; ln -s PARAM.in.test2.2nd PARAM.in; \
	mv restartOUT/*.nc restartIN/restart.nc; \
	mv restartOUT/*.txt restartIN/restart_info.txt; \
	${MPIRUN} ./ram_scb.exe | tee runlog2;	

test2_check:
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9	                      \
		${TESTDIR2}/output_ram/pressure_d20130317_t001500.dat \
		${IMDIR}/output/test1/pressure.ref                    \
		> test2.diff
	ncrcat ${TESTDIR2}/output_ram/sat1_d20130317_t000000.nc       \
	       ${TESTDIR2}/output_ram/sat1_d20130317_t001000.nc       \
	       ${TESTDIR2}/output_ram/sat1.nc
	ncdump -v "Flux_H","B_xyz" ${TESTDIR2}/output_ram/sat1.nc     \
               | sed -e '1,/data:/d' >                                \
               ${TESTDIR2}/output_ram/sat1.test        
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9                            \
                ${TESTDIR2}/output_ram/sat1.test                      \
                ${IMDIR}/output/test1/sat1.ref                        \
                >> test2.diff
	ncrcat ${TESTDIR2}/output_ram/sat2_d20130317_t000000.nc       \
               ${TESTDIR2}/output_ram/sat2_d20130317_t001000.nc       \
               ${TESTDIR2}/output_ram/sat2.nc
	ncdump -v "Flux_H","B_xyz" ${TESTDIR2}/output_ram/sat2.nc     \
               | sed -e '1,/data:/d' > \
               ${TESTDIR2}/output_ram/sat2.test
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9                            \
                ${TESTDIR2}/output_ram/sat2.test                      \
                ${IMDIR}/output/test1/sat2.ref                        \
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
	@echo "test3_check..." >> test3.diff
	make test3_check

test3_compile:
	make

test3_rundir:
	rm -rf ${TESTDIR3}
	make rundir RUNDIR=${TESTDIR3} STANDALONE="YES"
	cp Param/${PARAMFILE} ${TESTDIR3}/PARAM.in
	cp input/sat*.dat ${TESTDIR3}/

test3_run:
	cd ${TESTDIR3}; ${MPIRUN} ./ram_scb.exe | tee runlog;

test3_check:
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9                             \
		${TESTDIR3}/output_ram/pressure_d20130317_t001500.dat  \
		${IMDIR}/output/test3/pressure.ref      	       \
		> test3.diff                                           
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9			       \
		${TESTDIR3}/output_scb/hI_output_d20130317_t001500.dat \
		${IMDIR}/output/test3/hI.ref 		 	       \
		>> test3.diff
	@echo "Test Successful!"

#TEST 4----------------------------------
test4:
	@echo "starting..." > test4.diff
	@echo "test4_compile..." >> test4.diff
	make test4_compile
	@echo "test4_rundir..." >> test4.diff
	make test4_rundir PARAMFILE=PARAM.in.test4
	@echo "test4_run..." >> test4.diff
	make test4_run MPIRUN=
	@echo "test4_check..." >> test4.diff
	make test4_check

test4_compile:
	make

test4_rundir:
	rm -rf ${TESTDIR4}
	make rundir RUNDIR=${TESTDIR4} STANDALONE="YES"
	cp Param/${PARAMFILE}.* ${TESTDIR4}/
	cp input/sat*.dat ${TESTDIR4}/

test4_run:
	cd ${TESTDIR4};                                  \
	rm PARAM.in; ln -s PARAM.in.test4.1st PARAM.in;  \
	${MPIRUN} ./ram_scb.exe | tee runlog1;               \
	rm PARAM.in; ln -s PARAM.in.test4.2nd PARAM.in;  \
	mv restartOUT/*.nc restartIN/restart.nc; \
	mv restartOUT/*.txt restartIN/restart_info.txt; \
	${MPIRUN} ./ram_scb.exe | tee runlog2;      

test4_check:
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9                            \
                ${TESTDIR4}/output_ram/pressure_d20130317_t001500.dat \
                ${IMDIR}/output/test3/pressure.ref                    \
                > test4.diff
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9                             \
                ${TESTDIR4}/output_scb/hI_output_d20130317_t001500.dat \
                ${IMDIR}/output/test3/hI.ref                           \
                >> test4.diff
	ncrcat ${TESTDIR4}/output_ram/sat1_d20130317_t000000.nc       \
               ${TESTDIR4}/output_ram/sat1_d20130317_t001000.nc       \
               ${TESTDIR4}/output_ram/sat1.nc
	ncdump -v "Flux_H","B_xyz" ${TESTDIR4}/output_ram/sat1.nc     \
               | sed -e '1,/data:/d' >                                \
               ${TESTDIR4}/output_ram/sat1.test        
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9                            \
                ${TESTDIR4}/output_ram/sat1.test                      \
                ${IMDIR}/output/test3/sat1.ref                        \
                >> test4.diff
	ncrcat ${TESTDIR4}/output_ram/sat2_d20130317_t000000.nc       \
               ${TESTDIR4}/output_ram/sat2_d20130317_t001000.nc       \
               ${TESTDIR4}/output_ram/sat2.nc
	ncdump -v "Flux_H","B_xyz" ${TESTDIR4}/output_ram/sat2.nc     \
               | sed -e '1,/data:/d' >                                \
               ${TESTDIR4}/output_ram/sat2.test
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9                            \
                ${TESTDIR4}/output_ram/sat2.test                      \
                ${IMDIR}/output/test3/sat2.ref                        \
                >> test4.diff
	@echo "Test Successful!"

#TEST EMIC----------------------------------
testEMIC:
	@echo "starting..." > testEMIC.diff
	@echo "testEMIC_compile..." >> testEMIC.diff
	make testEMIC_compile
	@echo "testEMIC_rundir..." >> testEMIC.diff
	make testEMIC_rundir PARAMFILE=PARAM.in.testEMIC
	@echo "testEMIC_run..." >> testEMIC.diff
	make testEMIC_run MPIRUN=
	@echo "testEMIC_check..." >> testEMIC.diff
	make testEMIC_check

testEMIC_compile:
	make

testEMIC_rundir:
	rm -rf ${TESTDIRC}
	make rundir RUNDIR=${TESTDIRC} STANDALONE="YES"
	cp Param/${PARAMFILE} ${TESTDIRC}/PARAM.in
	cp input/sat*.dat ${TESTDIRC}/

testEMIC_run:
	cd ${TESTDIRC}; ${MPIRUN} ./ram_scb.exe | tee runlog;

testEMIC_check:
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9                              \
	        ${TESTDIRC}/output_ram/log_d20130317_t000000.log        \
	        ${IMDIR}/output/testEMIC/log.ref                           \
	        > testEMIC.diff
	${SCRIPTDIR}/DiffNum.pl -b -a=1e-9                              \
	        ${TESTDIRC}/output_ram/pressure_d20130317_t001500.dat   \
	        ${IMDIR}/output/testEMIC/pressure.ref                      \
	        >> testEMIC.diff
	@echo "Test Successful!"
