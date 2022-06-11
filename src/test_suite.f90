program run_test_suite
    use ModRamMain,      ONLY: nTestPassed, nTestRun
    use ModRamIndices,   ONLY: test_update_indices
    use ModRamFunctions, ONLY: test_Gcoul, test_NewtFit
    
    implicit none
    character(25) :: arg1
    logical :: verbose
  
    !First, make sure the right number of inputs have been provided
    if (command_argument_count().gt.1)THEN
      write(*,*) 'Error: Only one argument is currently supported'
      STOP
    end if
    
    ! check for -v, -verbose, etc. as command line argument
    verbose = .FALSE.  ! default is not verbose
    call get_command_argument(1, arg1)   ! read variable
    if ((index(arg1, '-v') == 1) .or. (index(arg1, '-V') == 1)) verbose = .TRUE.

    ! Tests go here...
    call test_update_indices(verbose)
    call test_Gcoul(verbose)
    call test_NewtFit(verbose)

    ! And after all the tests print the test status
    write(*,*) nTestPassed, "out of ", nTestRun, " tests passed"

    if (nTestPassed .NE. nTestRun) then
      STOP -1  ! if any tests fail, return a negative code so CI will report an error
    end if
  
end program run_test_suite
