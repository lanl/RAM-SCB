#!/usr/bin/env python

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from contextlib import contextmanager
import datetime as dt
import os
import sys
import re
import shutil
import subprocess


@contextmanager
def cd(newdir):
    '''Context-managed chdir

    changes back to original directory on exit or failure
    '''
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def setup_parser():
    # Create argument parser & set up arguments:
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-d", "--destdir", default=None, required=True,
                        help="Destination path for run directory")
    parser.add_argument("configfile", nargs=1, help="Run config")
    parser.add_argument("-o", "--overwrite", action="store_true",
                        help="Clobber an existing run directory" +
                        "with the same name/path.")
    parser.add_argument("--debug", action="store_true",
                        help="Print debug information.")
    # Handle arguments, noting that argparse expands linux wildcards.
    args = parser.parse_args()

    return args


def get_ramindices_end(fname='input/RamIndices.txt'):
    # read file and get last date
    with open(fname, 'r') as fh:
        inds = fh.readlines()
    lastline = inds[-1].strip()
    lastdatestr = lastline.split()[0]
    lastdate = dt.datetime.strptime(lastdatestr, '%Y%m%d')

    return lastdate


def parse_block(config_file, srch='#STARTTIME'):
    expr = re.compile(srch)
    with open(config_file, 'r') as fh:
        lines = fh.readlines()
        lines = [line.rstrip() for line in lines]
        for idx, line in enumerate(lines):
            searchyn = re.match(expr, line)
            if searchyn:  # target line
                useidx = idx
    store = [ll.split()[0] for ll in lines[useidx:useidx+7]]
    year = int(store[1])
    month = int(store[2])
    day = int(store[3])
    hour = int(store[4])
    minute = int(store[5])
    second = int(store[6])
    return dt.datetime(year, month, day, hour, minute, second)


def parse_config(config_file):
    '''Placeholder, pending web interface
    '''
    st_date = parse_block(config_file, srch='#STARTTIME')
    end_date = parse_block(config_file, srch='#STOPTIME')
    print('LAUNCH_RUN: Requested {} to {}'.format(st_date, end_date))
    lastramind = get_ramindices_end()
    print('LAUNCH_RUN: RamIndices ends at {}'.format(lastramind))
    sys.stdout.flush()
    if end_date >= lastramind:
        # if run ends after the last date in the Ramindices file,
        # update it
        with cd("Scripts"):
            subprocess.run(['python', 'updateRamIndices.py'])
    return st_date, end_date


def setup_rundir(args):
    '''Make, populate, and move RAM-SCB run directory
    '''
    # first check whether Kp/F10.7 need updating
    st_date, end_date = parse_config(args.configfile[0])
    # Now make rundir and copy in everything we need
    stYYMMDD = '{:04d}-{:02d}-{:02d}'.format(st_date.year,
                                             st_date.month,
                                             st_date.day)
    enYYMMDD = '{:04d}-{:02d}-{:02d}'.format(end_date.year,
                                             end_date.month,
                                             end_date.day)
    if args.overwrite:
        shutil.rmtree('run_ram_ror', ignore_errors=True)
    compl = subprocess.run(['make', 'rundir', 'RUNDIR=run_ram_ror'],
                           check=True, capture_output=True)
    # then make flux boundary files
    cmdline = ' '.join(['python', '/SHIELDS/flux-model/makeGEOboundary.py',
                        f'-s {stYYMMDD}', f'-e {enYYMMDD} -m 0',
                        '-r input -o run_ram_ror/input_ram'])
    compl = subprocess.run(cmdline, shell=True,
                           check=True, stdout=subprocess.PIPE)
    sys.stdout.flush()
    # add supplied PARAM file
    shutil.copyfile(args.configfile[0], 'run_ram_ror/PARAM.in')
    # and move rundir to final location
    if args.overwrite:
        shutil.rmtree(args.destdir, ignore_errors=True)
    shutil.move('run_ram_ror', args.destdir)


def run_model(args):
    '''Launch RAM-SCB and wait so the launch script won't exit
    '''
    with cd(args.destdir):
        print(os.getcwd())
        if not os.path.isfile('ram_scb.exe'):
            raise RuntimeError(' '.join(['RAM-SCB not found in specified',
                                         'directory, or directory not',
                                         'created properly']))
        # launch RAM-SCB, requires all relevant data in run dir
        # and appropriate subdirectories
        pram = subprocess.Popen(['./ram_scb.exe'])
        # entrypoint process in docker must be kept in foreground
        # to prevent container stopping
        pram.wait()


if __name__ == '__main__':
    args = setup_parser()
    setup_rundir(args)
    run_model(args)
