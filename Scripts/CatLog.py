#!/usr/bin/env python3
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
     -h or --help:      Display help info.
     -no or --nocheck:  Do not check for overlapping points, append as-is.
     -rm or --remove:   Remove all but the first file given.
     -o or --outfile:   Create new file with concatenated output.
     --debug:           Print debug information.

Examples:
    1) Append log_0002.log to log_0001.log, do not check for duplicate lines.
       > python CatLog.py --nocheck log_0001.log log_0002.log

    2) Combine all log[*].log files into a new, single file:
       > python CatLog.py -o=log_combined.log log*.log
'''

import re
from os import unlink
from shutil import copy
from argparse import ArgumentParser, RawDescriptionHelpFormatter


def get_time(line, index, debug=False):
    '''
    From a string entry, return the "time" the file was written.  This is
    done by splitting the line and taking all items corresponding to the
    indexes within the input list "index", concatenating them, and
    converting the resulting line into a float.  Many preceeding zeros are
    added to each entry to compensate for frequenty chances in within
    log files.
    '''

    parts = line.split()
    keep = ''

    for i in index:
        keep += '{:0>2}'.format(parts[i])

    if debug:
        print('TIME CHECK DEBUG:')
        print('Input Line="{}"'.format(line))
        print('Reducing to {}'.format(keep))

    return int(keep)


# Create argument parser & set up arguments:
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-o", "--outfile", default=None,
                    help="Rather than append to first file, create new file " +
                    "to store output.")
parser.add_argument("files", nargs='+', help="Files to convert.  Can be " +
                    "explicit files or a unix wildcard.")
parser.add_argument("-rm", "--remove", action="store_true",
                    help="Remove all but the first file given.")
parser.add_argument("-no", "--nocheck", action="store_true",
                    help="Do not check for overalpping points; append as-is.")
parser.add_argument("--debug", action="store_true",
                    help="Print debug information.")

# Handle arguments, noting that argparse expands linux wildcards.
args = parser.parse_args()

# Re-order files as the standard unix order will not reflect the
# true order for very large iteration runs.
pattern = '[_\-ned]{0,2}(\d+)[_\-ned]{0,2}(\d+)?\.[(log)|(sat)|(mag)]'
if len(args.files) > 1:
    # Get all file iterations.
    iters = {}
    for f in args.files:
        grab = re.search(pattern, f)
        if grab:
            if None in grab.groups():
                iters[f] = int(grab.groups()[0])
            else:
                iters[f] = int(''.join(grab.groups()))
        else:
            raise ValueError(
                'Could not find iteration number in filename {}'.format(f))
    # Order files by iteration numbers.
    args.files.sort(key=lambda f: iters[f])

if args.debug:
    print("File order = \n")
    for f in args.files:
        print("\t{}\n".format(f))
    if input('Is this correct? [y/N]') != 'y':
        raise Exception('Ordering error.')

# Create new file if "outname" is set:
if args.outfile:
    # Copy first file to new location:
    copy(args.files.pop(0), args.outfile)
else:
    # Append results to first file:
    args.outfile = args.files.pop(0)

# Open output file in append mode:
out = open(args.outfile, 'a+')

# Some machines don't do 'a+' correctly.  Rewind file as necessary.
if out.tell() > 0:
    out.seek(0, 0)

# Load and store header:
out.readline()  # garbage
head = out.readline()  # Header.
nbytes = len(out.readline())

# Using header info, create list of indexes corresponding to time.
time_locs = []  # List of indices
time_vars = []  # List of variable names (for debugging).
# Desired time variable names in order:
search_names = ['year', 'yyyy', 'yy', 'yr', 'doy', 'mo', 'month', 'mm', 'day',
                'dy', 'hour', 'hr', 'hh', 'mm', 'mn', 'min', 'ss', 'sec', 'sc']

iter_names = ['iter', 'it', 'nstep']
for i, part in enumerate(head.split()):
    for s in search_names:
        if s == part.lower():
            time_locs.append(i)
            time_vars.append(s)
            break

# If no time tags are found, try iterations:
if not time_locs:
    for i, part in enumerate(head.split()):
        for s in iter_names:
            if s == part.lower():
                time_locs.append(i)
                time_vars.append(s)
                break
        if time_locs:
            break  # Only want a single iteration tag.

# If nothing was found still, default to first column:
if not time_locs:
    time_locs.append(0)
    time_vars.append('Default (none found)')

if args.debug:
    print("DEBUG:\tOpened file {}" .format(out.name))
    print("\tEach line is {} characters long." .format(nbytes))
    print("\tHeader has {} entries." .format(len(head.split())))
    print("\tUsing the following columns in order for time calculation:")
    for i, s in zip(time_locs, time_vars):
        print("\t[{:02d}] {}".format(i, s))


# Load last line of file.
last_line = out.readlines()[-1]
# Get last time entry:
lasttime = get_time(last_line, time_locs, args.debug)

if args.debug:
    print("\tLast time = {}.".format(lasttime))

# Open rest of files, append.
for f in args.files:
    # No files that end with special characters.
    if f[-1] == '~':
        continue
    # Open file, slurp lines.
    if args.debug:
        print("Processing {}:".format(f))
    nextfile = open(f, 'r')
    lines = nextfile.readlines()
    nextfile.close()
    # Read header; skip this file if header is different.
    lines.pop(0)
    nowhead = lines.pop(0)
    if nowhead != head:
        if args.debug:
            print(head)
            print(nowhead)
            print("\tHeader does not match, discarding.")
        continue
    # Jump over overlapping lines:
    if not args.nocheck:
        nSkip = 0
        nLines = len(lines)
        while nSkip < nLines:
            if get_time(lines[0], time_locs) > lasttime:
                break
            else:
                lines.pop(0)
                nSkip += 1
        if args.debug:
            print("\tFound {} overlapping lines.".format(nSkip))
    # Append data to first log.
    if len(lines) < 1:
        continue
    for line in lines:
        out.write(line)
    # Save "last time".
    if not args.nocheck:
        lasttime = get_time(lines[-1], time_locs)
    # Delete file when done.
    if args.remove:
        if args.debug:
            print("\tRemoving file.")
        unlink(f)

# Finalize:
out.close()
