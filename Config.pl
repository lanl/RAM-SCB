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
    if(/^-ncdf=(.*)/)    {$libs{'netcdf'} =$1; $DoSetLibs = 1; next;};
    if(/^-pspline=(.*)/) {$libs{'pspline'}=$1; $DoSetLibs = 1; next;};
    if(/^-ncarg=(.*)/)   {$libs{'ngmath'} =$1; $DoSetLibs = 1; next;};
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
    $libs{'pspline'}=$ENV{PSPLINEDIR} unless exists $libs{'pspline'};
    $libs{'ngmath'}=$ENV{NCARGDIR} unless exists $libs{'ngmath'};
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
    foreach(keys(%libs)){
	if($libs{$_}){
	    $lib_cmd=$lib_cmd . "\\\n\t-L$libs{$_}/lib -l$_ ";
	    $lib_cmd=$lib_cmd . '-lnetcdff ' if($_ eq 'netcdf');
	    $lib_cmd=$lib_cmd . '-lezcdf '   if($_ eq 'pspline');
	    print FILE "$_  $libs{$_}\n";
	}else{
	    $lib_cmd=$lib_cmd . "\\\n\t-l$_ ";
	    $lib_cmd=$lib_cmd . '-lezcdf ' if($_ eq 'pspline');
	    print FILE "$_  (default lib path)\n";
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
    foreach(keys(%libs)){
	if(/pspline/){
	    $modpath = $libs{$_} ? "$libs{$_}/mod" : "/usr/lib/pspline/mod";
	    `echo PSPLINE_PATH = $modpath >> $MakefileDefEdit`;
	}
	if(/netcdf/){
	    $modpath = $libs{$_} ? "$libs{$_}/include" : "/usr/lib/netcdf/include";
	    `echo NETCDF_PATH = $modpath >> $MakefileDefEdit`;
	}
    }

}
#=============================================================================
#sub get_netcdf
#    # Grab the required NetCDF library from the given location.
#    # If ncdf_path is not given as an argument, the proper 
#    # ENV variable is used.
#{
#    # Check for NetCDF installation path.
#    $ncdf_path = $ENV{NETCDFDIR} unless $ncdf_path;
#    die "ERROR: Installation path for NetCDF not set!  See Help."
#	unless $ncdf_path;
#    
#    # Trim trailing '/' from path as necessary.
#    if($ncdf_path=~/(.*)\/$/){$ncdf_path=$1};
#    print "Using NetCDF install path: $ncdf_path\n";
#
#    # Create directory for NetCDF objects.
#    unless(-d 'srcNetcdf'){mkdir('srcNetcdf',0751) or die $!};
#    `cd srcNetcdf; cp $ncdf_path/lib/libnetcdf.a .; ar -x libnetcdf.a`;
#    `echo $ncdf_path > srcNetcdf/netcdf_install_path.txt`;
#}
##=============================================================================
#sub get_pspline
#    # Grab the required NetCDF library from the given location.
#    # If ncdf_path is not given as an argument, the proper 
#    # ENV variable is used.
#{
#    # Check for NetCDF installation path.
#    $pspl_path = $ENV{PSPLINEDIR} unless $pspl_path;
#    die "ERROR: Installation path for Pspline not set!  See Help."
#	unless $pspl_path;
#    
#    # Trim trailing '/' from path as necessary.
#    if($pspl_path=~/(.*)\/$/){$pspl_path=$1};
#    print "Using Pspline install path: $pspl_path\n";
#
#    # Create directory for Pspline objects.
#    unless(-d 'srcPspline'){mkdir('srcPspline',0751) or die $!};
#    `cd srcPspline; cp $pspl_path/lib/libpspline.a .; ar -x libpspline.a`;
#    `cd srcPspline; cp $pspl_path/lib/libezcdf.a .; ar -x libezcdf.a`;
#    `cd srcPspline; rm -f r8pspltsub.o pspltsub.o`; # Unneeded objects.
#    `cd srcPspline; cp $pspl_path/mod/*.mod .`;
#    `echo $pspl_path > srcPspline/pspline_install_path.txt`;
#}
#=============================================================================
sub print_help
    # Print RAM-SCB help.
{
    print "
Additional options for RAM-SCB/Config.pl:

-ncdf=(path)      Set installation path for NetCDF libraries.
-pspline=(path)   Set installation path for Pspline libraries. 

These MUST be set if environment variables PSPLINEDIR and
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
