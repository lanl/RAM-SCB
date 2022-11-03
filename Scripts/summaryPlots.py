from contextlib import contextmanager
import glob
from argparse import ArgumentParser
import os
import subprocess
import sys
import warnings

import matplotlib.pyplot as plt
from spacepy.pybats import ram
import spacepy.plot as splot
import spacepy.time as spt


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


class RamParser(ArgumentParser):
    def error(self, message):
        """Overload error handler to write help

        Useful in case of any error in usage"""
        sys.stderr.write('error: {}\n'.format(message))
        self.print_help()
        sys.exit(2)


def parserSetup():
    """
    Define and set up argument parser

    Returns
    -------
    parser : argparse.ArgumentParser object
    """
    # Define a command-line option parser and add the options we need
    parser = RamParser()
    parser.add_argument('rundir', metavar='rundir',
                        help='Run directory')
    parser.add_argument("-s", "--startTime", dest="startTime",
                        help="Start time for plots (YYYY-MM-DDThh:mm:ss)")
    parser.add_argument("-e", "--endTime", dest="endTime",
                        help="End time for plots (YYYY-MM-DDThh:mm:ss)")
    parser.add_argument("-o", "--outdir", dest="outdir",
                        help="Folder to save summary plots")

    return parser


def load_logs(options):
    logcands = sorted(glob.glob(os.path.join(options.rundir,
                                             'output_ram',
                                             'log*log')))
    outfile = os.path.join(options.rundir, 'output_ram',
                           'log_combined.log')
    if outfile in logcands:
        oidx = logcands.index(outfile)
        _ = logcands.pop(oidx)
    if len(logcands) > 1:
        subprocess.run(["python", "CatLog.py",
                        f"-o={outfile}",
                        *logcands])
        uselog = outfile
    else:
        uselog = logcands[0]
    log = ram.LogFile(uselog)
    return log


def plotDst(log, options):
    """Make default plot of RAM simulated Dst
    """
    fig = plt.figure(tight_layout=True, figsize=(8, 4.5))
    fig, ax = log.add_dst_quicklook(target=fig, showBiot=False)
    tstart = log['time'][0]
    tstart_ISO = spt.Ticktock(tstart).ISO[0]
    tend = log['time'][-1]
    tend_ISO = spt.Ticktock(tend).ISO[0]
    if options.startTime:
        opt_tst = spt.Ticktock(options.startTime).UTC[0]
        # start needs to be in valid range
        if (opt_tst > tstart) and (opt_tst < tend):
            tstart = opt_tst
        else:
            warnings.warn(f'Supplied start time ({options.startTime}) ' +
                          f'is out of valid range. Using {tstart_ISO}.')
    if options.endTime:
        opt_ten = spt.Ticktock(options.endTime).UTC[0]
        # end needs to be in valid range
        if (opt_ten > tstart) and (opt_ten < tend):
            tend = opt_ten
        else:
            warnings.warn(f'Supplied end time ({options.endTime}) ' +
                          f'is out of valid range. Using {tend_ISO}.')

    ax.set_xlim([tstart, tend])
    # output figure to either rundir/plots or custom
    runname = os.path.split(options.rundir)[-1]
    figname = 'ram_dst_{}.png'.format(runname)
    if options.outdir:
        if not os.path.isdir(options.outdir):
            os.makedirs(options.outdir)
        outfile = os.path.join(options.outdir, figname)
    else:
        plotdir = os.path.join(options.rundir, "plots")
        if not os.path.isdir(plotdir):
            os.makedirs(plotdir)
        outfile = os.path.join(plotdir, figname)
    plt.savefig(outfile)


if __name__ == "__main__":
    parser = parserSetup()
    in_args = parser.parse_args()

    bailflag = False
    if os.path.isdir(in_args.rundir):
        # make sure we don't have a trailing separator
        if (in_args.rundir[-1] == os.path.sep) and \
           (len(in_args.rundir) > 1):
            in_args.rundir = in_args.rundir[:-1]
        # run directory exists, check for presence of logs
        logcand = os.path.join(in_args.rundir, 'output_ram', 'log*log')
        if len(glob.glob(logcand)) > 0:
            splot.style('default')
            uselog = load_logs(in_args)
            plotDst(uselog, in_args)
    else:
        bailflag = True

    if bailflag:
        errmsg = 'Run directory {} does not exist '.format(in_args.rundir)
        errmsg += 'or has no log files'
        parser.error(errmsg)
