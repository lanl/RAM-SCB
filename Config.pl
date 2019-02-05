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

# Set identifier information; collect arguments.
our $Component       = "IM";
our $Code            = "RAM-SCB";
our $MakefileDefOrig = 'src/Makefile.def';
our @Arguments       = @ARGV;
my  $IsStandalone = 'True';
my  $DoSetLibs;
my  $SchemeUser = 'False';

foreach(@Arguments){
    if(/^-scheme/) {$SchemeUser = 'True';
                    push(@Arguments,"-compiler=gfortran");
                    push(@Arguments,"-mpi=openmpi");
                    push(@Arguments,"-openmp");  
                    #system("bash","-c","source schemeSetup.sh");
                    next;
    }
}
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

# Jump right to help if called:
&print_help if $Help;

# Handle RAM-SCB installation arguments.
foreach(@Arguments){
    if(/^-scheme/) {$libs{'gsl'} = "/packages2/.packages2/x86_64-pc-linux-gnu-rhel6/gsl/2.3";
                    $libs{'netcdf'} = "/projects/lanl/Carrington/netcdf";
                    $DoSetLibs = 1;
                    next;
    };
    if(/^-ncdf=(.*)/)    {$libs{'netcdf'} =$1; $DoSetLibs = 1; next;};
    if(/^-gsl=(.*)/)     {$libs{'gsl'}    =$1; $DoSetLibs = 1; next;};
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
    $libs{'netcdf'}=$ENV{NETCDFDIR} unless exists $libs{'netcdf'};
    $libs{'gsl'}=$ENV{GSLDIR} unless exists $libs{'gsl'};
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
    if($SchemeUser eq 'True'){
        $lib_cmd=$lib_cmd . "\\\n\t-L$libs{'netcdf'}/lib -lnetcdff";
        $lib_cmd=$lib_cmd . "\\\n\t-L$libs{'gsl'}/lib -lgsl -lgslcblas -lm";
    }else{
        foreach(keys(%libs)){
	    if($libs{$_}){
	        $lib_cmd=$lib_cmd . "\\\n\t-L$libs{$_}/lib -l$_ ";
    	        $lib_cmd=$lib_cmd . '-lnetcdff ' if($_ eq 'netcdf');
                $lib_cmd=$lib_cmd . '-lgslcblas -lm ' if($_ eq 'gsl');
    	        print FILE "$_  $libs{$_}\n";
    	    }else{
    	        $lib_cmd=$lib_cmd . "\\\n\t-l$_ ";
                $lib_cmd=$lib_cmd . '-lnetcdff ' if($_ eq 'netcdf');
    	        print FILE "$_  (default lib path)\n";
    	    }
        }
    }

    close(FILE);

    # Choose which makefile.conf to edit based on IsStandalone
    my $MakefileConfEdit = 'Makefile.conf';
    my $MakefileDefEdit  = 'Makefile.def';
    $MakefileConfEdit = '../../Makefile.conf' if($IsStandalone eq 'False');
    $MakefileDefEdit  = '../../Makefile.def'  if($IsStandalone eq 'False');

    print "Editing $MakefileConfEdit\n";

    # Open Makefile.conf, insert libflags before Lflag1.
    @ARGV = ($MakefileConfEdit);
    while(<>){
	die "ERROR: Library locations are already set.  Reinstall code first\n"
	    if(/^\s*LibFlags/);
	s/^(\s*Lflag1\s*=.*)\n/$lib_cmd\n$1$add_flg\n/;
	s/^(\s*Lflag2\s*=.*)\n/$1$add_flg\n/;
	print;
    }

    # Some libraries have modules that need to be included during compile.
    # Add the locations of these modules to the local Makefile.def
    my $modpath = '';
    `echo  >> $MakefileDefEdit`;   # Add some spaces
    `echo  >> $MakefileDefEdit`;   # to end of file...
    if($SchemeUser eq 'True') {
        $modpath = "$libs{'netcdf'}/include";
       `echo NETCDF_PATH = $modpath >> $MakefileDefEdit`;
        $modpath = "$libs{'gsl'}/include";
       `echo GSL_PATH = $modpath >> $MakefileDefEdit`;
    }else{
        foreach(keys(%libs)){
	    if(/netcdf/){
	        $modpath = $libs{$_} ? "$libs{$_}/include" : "/usr/lib/netcdf/include";
           `    echo NETCDF_PATH = $modpath >> $MakefileDefEdit`;
	    }
            if(/gsl/){
                $modpath = $libs{$_} ? "$libs{$_}/include" : "/usr/lib/gsl/include";
           `    echo GSL_PATH = $modpath >> $MakefileDefEdit`;
            }
        }
    }

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
multiple installations of the libraries.
    
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
