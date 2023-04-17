from contextlib import contextmanager
import datetime as dt
import glob
from argparse import ArgumentParser
import os
import re
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
    # move tstart/tend into options struct
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


def parse_ptime(in_str):
    '''Parse date/time in pressure file names and return datetime
    '''
    tparts = re.search('pressure_d(\d{4})(\d{2})(\d{2})_t(\d{2})(\d{2})(\d{2})',
                       in_str).groups()
    return dt.datetime(*[int(tp) for tp in tparts])


def plot_pressure(options):
    presscands = sorted(glob.glob(os.path.join(options.rundir,
                                               'output_ram',
                                               'press*dat')))
    # loop over all files, if in expected time range make plot
    for pcfn in presscands:
        fpath, fmain = os.path.split(pcfn)
        fmain = fmain.split('.')[0]
        ftime = parse_ptime(fmain)
        opt_tst = spt.Ticktock(options.startTime).UTC[0]
        opt_ten = spt.Ticktock(options.endTime).UTC[0]
        if (ftime <= opt_ten) and (ftime >= opt_tst):
            pressf = ram.PressureFile(pcfn)
            plotlist = [key for key in pressf if key.startswith('tot')]
            # Convenience plotting routines in spacepy currently assume a log-transform
            # but this is not useful for anisotropy. When control for that is added to
            # spacepy we can do this more easily...
            # plotlist = [key for key in pressf if (key.startswith('tot') or key.startswith('ani'))]
            for pl in plotlist:
                # TODO: maybe make combined plots, look at axis limits, etc.
                fig, ax, cm, _ = pressf.add_pcol_press(var=pl, add_cbar=True)
                outfn = os.path.join(options.outdir, fmain + f'_{pl}.png')
                plt.savefig(outfn, dpi=200)
                plt.close()


if __name__ == "__main__":
    parser = parserSetup()
    in_args = parser.parse_args()

    flagstatus = [False, False, False]
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
            flagstatus[1] = True
        pfcand = os.path.join(in_args.rundir, 'output_ram', 'press*')
        if len(glob.glob(pfcand)) > 0:
            plot_pressure(in_args)
    else:
        flagstatus[0] = True

    if True in flagstatus:
        errmsg = ''
        if flagstatus[0]:
            errmsg += '\nRun directory {} does not exist\n'.format(in_args.rundir)
        if flagstatus[1]:
            errmsg += '\nNo log files found\n'
        if flagstatus[2]:
            errmsg += '\nNo pressure files found\n'
        parser.error(errmsg)
