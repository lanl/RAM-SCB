#!/usr/bin/env python

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from contextlib import contextmanager
import datetime as dt
import os
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


def parse_config(config_file):
    '''Placeholder, pending web interface
    '''
    end_date = dt.datetime(2015, 3, 18)  # get from config file
    lastramind = get_ramindices_end()
    if end_date >= lastramind:
        # if run ends after the last date in the Ramindices file,
        # update it
        with cd("Scripts"):
            subprocess.run(['python', 'updateRamIndices.py'])


def setup_rundir(args):
    '''Make, populate, and move RAM-SCB run directory
    '''
    # first check whether Kp/F10.7 need updating
    parse_config(args.configfile)
    # Now copy in everything we need
    # and move rundir to final location
    compl = subprocess.run(['make', 'rundir', 'RUNDIR=run_ram_ror'],
                           check=True, capture_output=True)
    shutil.copyfile(args.configfile[0], 'run_ram_ror/PARAM.in')
    if args.overwrite:
        shutil.rmtree(args.destdir)
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
