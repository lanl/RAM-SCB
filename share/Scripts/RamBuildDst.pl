#!/usr/bin/perl -w
#RamBuildDst.pl

#===========================================================================
#RamBuildDst                                                        DTW 2009
#
# Quickly extracts and sums all Dst values from ram###ss.in files.
#
# Use: Run in output_ram directory.  
#
#===========================================================================

use strict;


# Get lists of files that have the Dst values.
my @h_list = (<ram_h_d*_t*.l>);
my @o_list = (<ram_o_d*_t*.l>);
my @helist = (<ram_he_d*_t*.l>);

# Ensure there are the same number of files for each.
my $nFiles = $#h_list;
die "ERROR: Number of files for each species does not match!\n"
    unless (($nFiles == $#o_list) and ($nFiles == $#helist));

my $dst_sum = 0;
my $dst_h   = 0;
my $dst_he  = 0;
my $dst_o   = 0;
my $h_time;
my $hetime;
my $o_time;

# Create output file and open; write header.
$h_list[0] =~ /^ram_h_(d\d{8}_t\d{6})\.l/;
die "ERROR: Could not interpret file name $h_list[0]\n" 
    unless($1);
open(DST,">dst_$1.txt");
print DST "Header TBD\n";

foreach my $i (0..$nFiles){
    # Check that file times match.
    $h_list[$i]=~/ram_h_(d\d{8}_t\d{6})\.l/;
    $h_time = $1;
    $helist[$i]=~/ram_he_(d\d{8}_t\d{6})\.l/;
    $hetime = $1;
    $o_list[$i]=~/ram_o_(d\d{8}_t\d{6})\.l/;
    $o_time = $1;
    die "ERROR: Mismatching times for files stamped $h_time\n"
	if( ($h_time ne $hetime) or ($h_time ne $o_time) );
    
    # Get dst from each species file.
    $dst_h = &extract_dst($h_list[$i]);
    $dst_he= &extract_dst($helist[$i]);
    $dst_o = &extract_dst($o_list[$i]);
    
    $dst_sum = $dst_h + $dst_he + $dst_o;

    # Write date, time, and dst to file.
    printf DST "%s %E %E %E %E\n",$h_time, $dst_h, $dst_he, $dst_o, $dst_sum;
}

close(DST);
print "Finished!\n";

#===========================================================================
sub extract_dst
    # Extract dst value from input file.
{
    my $file = shift;
    my $dst = 0.0;
    open(INFILE, '<', $file);
    while(<INFILE>){
	if(/Dst\s*=\s*-*\s*(\d+\.\d+E[+-]\d+)/){
	    $dst = $1;
	    $dst = $dst * -1.0 if /=\s-\s/;
	}
    }
    close(INFILE);
    return $dst;
}
#===========================================================================
