# Contributing to RAM-SCBE

Bug reports, suggested enhancements, and other contributions
are always appreciated! RAM-SCBE is free, open-source software
for the space physics community and welcomes both use and
engagement.

## Bug Reports and contributions

### Filing issues

When filing an issue, please include:
- Your operating system name/version
- Which compilers and dependent packages you're using
- Detailed steps, preferably on a minimal example, to reproduce your issue
- Log files and/or full error messages

### Contributing code

Contributions to RAM-SCBE are managed through pull requests
on [our GitHub repository](https://github.com/lanl/RAM-SCB).
Further details on the workflow are given in <WORKFLOW.md>.

If you have any questions about how to write the code, either
comment on the open issue that your code addresses or submit a
draft pull request and ask for help.

## Coding style and standards

RAM-SCBE has developed and evolved over a long time period.
The code has been modernized during updates, but at this time
different coding styles and standards are seen across the code.
As further developments are made, we are encouraging the Fortran
2008 standard.

- RAM-SCBE uses free-format Fortran code
- When using a module, please use `use ModNameHere, ONLY: var, sub` to list the members being used
- Functions are preferred where code is not dependent on state parameters
  - "New-style" functions should be used, such that the result variable is declared on the same line as the function and there is no `return` statement at the end of the function.
- Subroutines should minimize the number of used state parameters
- RAM-SCBE uses a two-space indent, and spaces should be used in-code for readability
- Please provide copious comments. They should be descptive (what is this doing?) as well as explanatory (why is it doing this?)
- When contributing new code, please extend the test suite (see next section)
- When contributing code, please update the documentation as appropriate.
- *Do not* use Fortran 77 features, including:
  - Common blocks
  - GOTO statements
  - Implicit variable typing
  - Statement labels (i.e., numbered statements for flow control)

## Testing in RAM-SCBE

RAM-SCB has two primary approaches to testing. The first is integration
testing, where a short model run is defined and the output is compared
to a known result. This serves partly as a regression test, and partly
to ensure that the simulation runs end-to-end and gives a sensible result.
These tests are run via the Makefile by, e.g., `make test1`

The second approach to testing is unit testing. Rather than define a run and
use the computational resources required for even short simulations, unit
testing concentrates on small pieces of code and checks them, in isolation,
for correctness. This complements the integration approach, and the unit test
suite can be run via the Makefile with `make unittest`

### Adding a new unit test

Ideally, new code will have unit test coverage. Depending on your preferred style,
the unit tests can be written first (so the code will have a pre-defined calling
syntax when you go to write it), or after the fact. Adding tests to the unit test
suite involves two steps:
1. Write the test subroutine at the end of the module where the code being tested is.
2. Add a line calling the test subroutine in `src/test_suite.f90`

The test subroutine should update the module variables `nTestPassed` and `nTestRun`,
from `ModRamMain`. A sample stub for a test is given below:
```
subroutine test_npifunction(verbose)
  ! Comments describing your test(s) go here
  ! 1. Test that npifunction returns a correct approximation to n*pi
  use ModRamMain,    ONLY: Real8_, test_neq_abs, nTestPassed, nTestRun

  implicit none
  logical, intent(in) :: verbose
  logical :: failure, all_pass
  integer :: idx
  real(Real8_), dimension(3) :: expect = (/3.14, 6.28, 9.42/)
  real(Real8_), dimension(3) :: answer

  ! Call the new code and test that everything passes
  all_pass = .TRUE.  ! Start as true and flip to false if anything fails
  do idx=1,3
    answer(idx) = npifunction(real(idx))
    call test_neq_abs(answer(idx), expect(idx), failure)
    ! The failure variable will be .TRUE. if answer and expect are not equal
    if (failure) all_pass = .FALSE.
  end do
  ! set verbose output, so we can inspect what failed if required
  if (verbose) then
    write(*,"(A25,3F12.8))") "test_npifunction: got = ", answer
    write(*,"(A25,3F12.8))") "             expected = ", expect
  end if
  ! and now update the test counters
  nTestRun = nTestRun + 1
  if (all_pass) then
    nTestPassed = nTestPassed + 1
  end if
  ! if required, additional tests for the same function can be
  ! put here, just make sure the counters are updated after every test
end subroutine test_npifunction
```

### Adding a new integration test
Major new functionality should have a short integration test that exercises it.
Configurations for the integration tests are given in <Param> and are named for
the test. E.g., <Param/PARAM.in.test1> or <Param/PARAM.in.testEMIC>. Expected
test results are stored in the appropriately named folder in <output>, e.g.,
<output/test1> or <output/testEMIC> for the examples given previously.