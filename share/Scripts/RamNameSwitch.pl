#!/usr/bin/perl -w
#RamNameSwitch.pl

#===========================================================================
#RamNameSwitch                                                      DTW 2009
#
#This script switches all output file names from the old format to the new
#format (or visa versa).
#
# RAM-SCB FILE NAME FORMAT STANDARDS:
#  (prefix)_(extra descriptor)_dYYYYMMDD_tHHMMSS.(extension)
#
# Use: Run in output_ram directory.  
#
#===========================================================================

use strict;

# Declare some variables.
my $yy; my $mm; my $dd; my $hh; my $mn;
my $date_str;


#-------------------------------
# Get lists of files to rename.
#-------------------------------
my @files_prs = (<pressure_*.in>); # Pressure files.
my @files_ram = (<ram*.*>);        # ram*.t and ram*.l files.


#-------------------------------
# Get start time info from user
#-------------------------------
($yy,$mm,$dd,$hh) = &get_epoch;
$mn = 00;
print "Simulation start: $yy/$mm/$dd $hh:$mn UT\n";

#-------------------------------
# Change pressure file names.
#-------------------------------
foreach(@files_prs){
    next unless /pressure_(\d{4})\.in/;  # Extract number from filename.
    die "ERROR!  Could not interpret $_\n" unless($1);
    $date_str = &make_date(0,5*$1);

    # Debug:
    #print "File $_ has date-time of $date_str\n";

    rename($_, "pressure_$date_str.in");
}

foreach(@files_ram){
    next unless /^ram(\d{3})_*(\w{1,2})\.(\w)/; # Extract species and time.
    die "ERROR!  Could not interpret $_\n" unless($1 and $2 and $3);
    $date_str = &make_date($1,0);
    
    # Debug:
    #print "Replacing $_ with ram_$2\_$date_str\.$3\n";

    rename($_, "ram_$2\_$date_str\.$3");
}

#===========================================================================
sub get_epoch
    # Get start date/time from user.
{
    my $yy_in; my $mm_in; my $dd_in; my $hh_in;

    print "Please enter start date and time for this simulation.\n";
    while(1){
	print "Year (YYYY):";
	chomp($yy_in=(<STDIN>));
	last if ($yy_in=~/\d{4}/);
	print "Wrong format! Please try again. \n";
    }
    while(1){
	print "Month (MM):";
	chomp($mm_in=(<STDIN>));
	if($mm_in > 12){print "$mm_in is too big!\n";  next;}
	last if ($mm_in=~/\d{2}/);
	print "Wrong format! Please try again.\n";
    }
    while(1){
	print "Day (DD):";
	chomp($dd_in=(<STDIN>));
	if($dd_in > 31){print "$dd_in is too big!\n";  next;}
	last if ($dd_in=~/\d{2}/);
	print "Wrong format! Please try again.\n";
    }
    while(1){
	print "Hour (HH):";
	chomp($hh_in=(<STDIN>));
	if($hh_in > 24){print "$hh_in is too big!\n";  next;}
	last if ($hh_in=~/\d{2}/);
	print "Wrong format! Please try again.\n";
    }	

    return ($yy_in,$mm_in,$dd_in,$hh_in);
}
    
#===========================================================================
sub make_date
    # Use start time + minutes elapsed + hours elapsed to build time string.
{
    my $hours= shift; 
    my $mins = shift;

    my $IsLeap;
    $IsLeap=1 unless($yy % 4.0);

    my %DaysInMonth  = (
	'01' => 31,
	'02' => 28,
	'03' => 31,
	'04' => 30,
	'05' => 31,
	'06' => 30,
	'07' => 31,
	'08' => 31,
	'09' => 30,
	'10' => 31,
	'11' => 30,
	'12' => 31,
	);
    $DaysInMonth{'02'}=29 if($IsLeap);

    # Time of day:
    my $m_total = $mins + $mn;
    my $m1 = $m_total % 60.0;
    my $h_total = $hours + $hh + int(($m_total - $m1)/60.0);
    my $h1 = $h_total % 24.0;

    # Date:
    my $y1 = $yy;
    my $M1 = $mm;
    my $d1 = $dd + int(($h_total - $h1)/24.0);
    if($d1 > $DaysInMonth{$mm}){
	$d1 = $d1 - $DaysInMonth{$mm};
	$M1++;
    }

    if($M1 > 12.0){
	$M1 = $M1 - 12.0;
	$y1++;
    }

    # Build string and return.
#    my $date_str = 
    return sprintf 
	"d%04d%02d%02d_t%02d%02d%02d", $y1, $M1, $d1, $h1, $m1, 00;
}

#===========================================================================
