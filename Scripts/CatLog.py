#!/usr/bin/env python
'''
Concatenate multiple logfiles with the generic SWMF format (flat, column-
organized ascii) into a single log file that does not contain any over-
lapping points.  It assumes that the first column is either time elapsed or
total run iteration.  If any file does not have the same header as the first,
it is disregarded.

Disclaimer: Always check your appended files.  Incomplete files, extra 
newlines, and mismatched files can yield unexpected results.

Usage: CatLog.py [options] log1 log2 [log3] [log4]...[logN]

The contents of log files 2 through N will be appended to log1 in order.
If a linux-like wild card syntax is provided, a file-glob will be used and
files sorted by name.

Options:
     -h or -help: Display this help info.
     -nocheck:    Do not check for overlapping points, append as-is.
     -rm:         Remove all but the first file given.
     -debug:      Print debug information.

Examples:
    1) Append log_0002.log to log_0001.log, do not check for duplicate lines.
       >CatLog.py -nocheck log_0001.log log_0002.log

    2) Combine all files that begin with 'sat_cluster_n' and remove all but one:
       >CatLog.py -rm sat_cluster*
'''

import sys
from glob import glob
from os import unlink

# Declare important variables.
check=True
debug=False
remove=False
files=[]

#TODO: replace with argparse option parser
for option in sys.argv[1:]:
    # Handle options.
    if option[0]=='-':
        if option == '-nocheck':
            check = False
        elif option == '-rm':
            remove = 'True'
        elif option == '-debug':
            debug = True
        elif option == '-h' or option == '-help':
            print(__doc__)
            exit()
        else:
            print('Unrecognized option: ', option)
            print(__doc__)
            exit()
    else:
        files = files + glob(option)

# Open output file in append mode:
out = open(files.pop(0), 'a+')

# Some machines don't do 'a+' correctly.  Rewind file as necessary.
if out.tell()>0: out.seek(0,0)

# Load and store header:
out.readline() #garbage
head = out.readline() # Header.
nbytes = len(out.readline())
if debug:
    print("DEBUG:\tOpened file %s" % (out.name))
    print("\tEach line is %i characters long." % nbytes)
    print("\tHeader has %i entries." % (len(head.split())) )

# Load last line of file.
out.seek(-1*nbytes, 2) #seek last line.
lasttime = float((out.readline().split())[0])

if debug:
    print("\tLast time = %f." % lasttime)

# Open rest of files, append.
for f in files:
    # No files that end with special characters.
    if f[-1]=='~': continue
    # Open file, slurp lines.
    if debug: print("Processing %s:" % f)
    nextfile = open(f, 'r')
    lines = nextfile.readlines()
    nextfile.close()
    # Read header; skip this file if header is different.
    lines.pop(0)
    nowhead = lines.pop(0)
    if nowhead != head: 
        if debug: 
            print(head)
            print(nowhead)
            print("\tHeader does not match, discarding.")
        continue
    # Jump over overlapping lines:
    if check:
        nSkip=0
        nLines=len(lines)
        while nSkip<nLines:
            if float( (lines[0].split())[0] ) > lasttime:
                break
            else:
                lines.pop(0)
                nSkip += 1
        if debug:
            print("\tFound %i overlapping lines." % nSkip)
    # Append data to first log.
    if len(lines)<1:
        continue
    for l in lines:
        out.write(l)
    # Save "last time".
    lasttime=float( (lines[-1].split())[0] )
    # Delete file when done.
    if remove: 
        if debug:
            print("\tRemoving file.")
        unlink(f)

# Finalize:
out.close()
