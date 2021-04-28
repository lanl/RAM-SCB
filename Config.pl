#!/usr/bin/perl -i
#Config.pl - The configuration script for RAM-SCB.
#This requires the main Config.pl script, located in: 
#   share/Scripts/        (Standalone mode)
#   ../../share/Scripts   (Framework Module mode)
#
#Note that, different from other module's Config.pls,
#much work is dedicated to the installation of 
#standalone libraries used by RAM-SCB.

use strict;
use CPAN::Version;

# Set identifier information; collect arguments.
our $Component       = "IM";
our $Code            = "RAM-SCB";
our $MakefileDefOrig = 'src/Makefile.def';
our @Arguments       = @ARGV;
my  $IsStandalone = 'True';
my  $DoSetLibs;
my  $DoAutoGSL;
my  $DoAutoNCF;
sub check_exists;
sub gfort_ge_10;

# Find and execute main Config.pl.
# Simultaneously determine if in stand-alone mode.
my $config = "../../share/Scripts/Config.pl";
if(-f $config){
    require $config;
    $IsStandalone = 'False';
}else{
    require "share/Scripts/Config.pl";
}

# Other variables inherited from share/Scripts/Config.pl
our %Remaining;   # Arguments not handled by share/Scripts/Config.pl
our $Help;        # Show help.
our $Verbose;     # Verbose info is printed if true.
our $Install;
our $Show;

# Declare local variables.
my %libs;
my $compiler;
my $extra_flags = "";

# Jump right to help if called:
&print_help if $Help;

# Handle RAM-SCB installation arguments.
foreach(@Arguments){
    if(/^-compiler/) {
        if(/^-compiler=(.+)/) {
            $compiler = $1;
            my $prefix = "gfortran";
            if (substr($compiler, 0, length $prefix) eq $prefix) {
                # If compiler name starts with gfortran.
                if (gfort_ge_10()) {
                    $extra_flags = $extra_flags . " -fallow-argument-mismatch";
                    # Remove following line when share library is updated to remove octal mask setting
                    $extra_flags = $extra_flags . " -fallow-invalid-boz";
                }
            }
        }
        next;
    }
    if(/^-ncdf/) {
        if(/^-ncdf=(.+)/) {
            # Value supplied
            $libs{'netcdf'} =$1;
        } else {
            # No value supplied. try to autodetect
            # Requires nf-config
            check_exists 'nf-config' or die "$0 requires nf-config";
            $DoAutoNCF = 1;
        }
        $DoSetLibs = 1;
        next;};
    if(/^-gsl/) {
        if(/^-gsl=(.+)/) {
            $libs{'gsl'} = $1;
        } else {
            # No value supplied. try to autodetect
            # Requires gsl-config
            check_exists 'gsl-config' or die "$0 requires gsl-config";
            $DoAutoGSL = 1;
        }
        $DoSetLibs = 1;
        next;};
    if(/^-setlibs/)      {$DoSetLibs = 1;  next;};
}

&get_settings;
&show_settings if $Show;

# LIBRARY LOCATIONS::
$DoSetLibs = 1 if($Install);
#die "RAM-SCB ERROR: Change library locations ONLY on installation!\n"
#    if($DoSetLibs and ! $Install);
# Set up library locations.
if($DoSetLibs){
    # Set default library locations unless set by flags.
    if ($DoAutoNCF) {
        print "Auto-detecting NetCDF-Fortran location\n";
        my $npre = `nf-config --prefix`;
        chomp($npre);
        $libs{'netcdf'}=$npre;
    } else {
        $libs{'netcdf'}=$ENV{NETCDFDIR} unless exists $libs{'netcdf'};
    }
    if ($DoAutoGSL) {
        print "Auto-detecting GSL location\n";
        my $gpre = `gsl-config --prefix`;
        chomp($gpre);
        $libs{'gsl'}=$gpre;
    } else {
        $libs{'gsl'}=$ENV{GSLDIR} unless exists $libs{'gsl'};
    }
    &set_libs;
}

exit 0;

#=============================================================================
sub set_libs
    # Set library locations and edit Makefiles appropriately.
    # This subroutine should ONLY be called on an install/reinstall.
{
    print "Setting library locations for RAM_SCB\n";
    
    # Echo libraries used if $Verbose flag set.
    if($Verbose){print "$_ is found at $libs{$_}\n" foreach(keys(%libs))};

    # Build Libflags declaration.  Write locations to file.
    my $add_flg = ' ${Libflags} ';
    my $lib_cmd = "Libflags = ";
    if($IsStandalone eq 'True'){
	open(FILE, '>', 'LibLocations.txt');
    }else{
	open(FILE, '>', '../../LibLocations.txt');
    }
    foreach (keys(%libs)) {
	    if($libs{$_}){
            if (($_ eq 'netcdf') and $DoAutoNCF) {
                my $nclibs = `nf-config --flibs`;
                chomp($nclibs);
             $lib_cmd = $lib_cmd . $nclibs . " ";
            } elsif ($_ eq 'netcdf') {
                $lib_cmd=$lib_cmd . "\\\n\t-L$libs{$_}/lib -l$_ ";
             $lib_cmd=$lib_cmd . '-lnetcdff ';
            }
            if (($_ eq 'gsl') and $DoAutoGSL) {
                my $gsllibs = `gsl-config --libs`;
                chomp($gsllibs);
                $lib_cmd = $lib_cmd . $gsllibs. " ";
            } elsif ($_ eq 'gsl') {
                $lib_cmd=$lib_cmd . "\\\n\t-L$libs{$_}/lib -l$_ ";
                $lib_cmd=$lib_cmd . '-lgslcblas -lm ';
            }
            print FILE "$_  $libs{$_}\n";
        } else {
            $lib_cmd=$lib_cmd . "\\\n\t-l$_ ";
            $lib_cmd=$lib_cmd . '-lnetcdff ' if($_ eq 'netcdf');
            print FILE "$_  (default lib path)\n";
     }
    }

    close(FILE);

    # Choose which makefile.conf to edit based on IsStandalone
    my $MakefileConfEdit = 'Makefile.conf';
    my $MakefileDefEdit  = 'Makefile.def';
    $MakefileConfEdit = '../../Makefile.conf' if($IsStandalone eq 'False');
    $MakefileDefEdit  = '../../Makefile.def'  if($IsStandalone eq 'False');

    print "\nEditing $MakefileConfEdit\n";

    # Open Makefile.conf, insert libflags before Lflag1.
    @ARGV = ($MakefileConfEdit);
    while(<>){
	die "ERROR: Library locations are already set.  Reinstall code first\n"
	    if(/^\s*LibFlags/);
	s/^(\s*CFLAG\s*=.*)\n/$1$extra_flags\n/;
	s/^(\s*Lflag1\s*=.*)\n/$lib_cmd\n$1$add_flg\n/;
	s/^(\s*Lflag2\s*=.*)\n/$1$add_flg\n/;
	print;
    }

    # Some libraries have modules that need to be included during compile.
    # Add the locations of these modules to the local Makefile.def
    my $modpath = '';
    my $gsl_include = '';
    `echo  >> $MakefileDefEdit`;   # Add some spaces
    `echo  >> $MakefileDefEdit`;   # to end of file...
    foreach (keys(%libs)) {
	    if (/netcdf/) {
            if ($DoAutoNCF) {
	            $modpath = `nf-config --includedir`;
                chomp($modpath);
            } else {
	            $modpath = $libs{$_} ? "$libs{$_}/include" : `nc-config --includedir`;
            }
            `echo NETCDF_PATH = $modpath >> $MakefileDefEdit`;
	    }
        if(/gsl/){
            if ($DoAutoGSL) {
                $gsl_include = `gsl-config --prefix`;
                chomp($gsl_include);
                $gsl_include = ${gsl_include}."/include";
                $modpath = $gsl_include;
            } else {
                $modpath = $libs{$_} ? "$libs{$_}/include" : "/usr/lib/gsl/include";
                $gsl_include = $modpath;
            }
            `echo GSL_PATH = $modpath >> $MakefileDefEdit`;
        }
    }

    # The C code using GSL needs the include location too. The path is seemingly not
    # always provided, so we add it explicitly here in Makefile.conf.
    @ARGV = ($MakefileConfEdit);
    while(<>){
	    s/^(\s*FLAGC\s*=.*)\n/FLAGC_EXTRA = -I$gsl_include\n$1$\\n/;
	    print;
    }
}

sub check_exists { 
    my $check = `sh -c 'command -v $_[0]'`; 
    return $check;
}

sub gfort_ge_10 {
    # Test gfortran version. Return 1 if version GE 10
    my $verout = `$compiler -dumpversion`;
    return CPAN::Version->vge("$verout","10") ? 1 : 0;
}

#=============================================================================
sub print_help
    # Print RAM-SCB help.
{
    print "
Additional options for RAM-SCB/Config.pl:

-ncdf=(path)      Set installation path for NetCDF libraries.
-gsl=(path)   Set installation path for GSL libraries. 

These MUST be set if environment variables GSLDIR and
NETCDFDIR are not set.  If both flag and env variable are
set, Config.pl will use the flag value.  This allows for
multiple installations of the libraries. If the flags are set
no path is given (e.g., '-ncdf') then Config.pl will use the
libraries own utility to look up required information.
    
Example installation for RAM-SCB standalone:

     Config.pl -install -ncdf=~/NetCDF/ -mpi=mpich2 -single -compiler=pgf90

";
    exit 0;
}
#=============================================================================
sub get_settings
    # Get current RAM-SCB specific settings.
{
}
#=============================================================================
sub show_settings
    # Display RAM-SCB specific settings.
{
}
#=============================================================================
