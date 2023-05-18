from argparse import ArgumentParser
from contextlib import contextmanager
import datetime as dt
import glob
import itertools as it
import os
import re
import subprocess
import sys
import warnings

import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
import spacepy.datamodel as dm
from spacepy.pybats import PbData
from spacepy.pybats import ram
import spacepy.plot as splot
import spacepy.time as spt
import spacepy.toolbox as tb


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


class Flux2DFile(PbData):
    def __init__(self, filename):
        # Filename:
        super().__init__()
        self.filename = filename
        self._read()

    def _read(self):
        '''
        Load the NCDF4 file.  Should only be called by __init__().
        '''
        with nc4.Dataset(self.filename, "r", format="NETCDF4") as ncf:
            # global attrs
            for attname in ncf.ncattrs():
                self.attrs[attname] = ncf.getncattr(attname)
            # get copies of all variables
            for varname in ncf.variables:
                self[varname] = dm.dmarray(ncf.variables[varname][...].filled())
                # RAM-SCB variables don't have metadata yet, so skip
                # RAM-SCB also doesn't have nested groups, so don't recurse
            # and add netCDF4 "dimensions" as variables
            for dimname in ncf.dimensions:
                self[dimname] = dm.dmarray(np.atleast_1d(ncf.dimensions[dimname].size))
            self.namevars = list(ncf.variables.keys())

        # Get start time, set default if not found.
        try:
            stringtime = self.attrs['start_time'].decode()
        except AttributeError:
            # unicode str already
            stringtime = self.attrs['start_time']
        except KeyError:
            stringtime = '20000101000000'

        # parse start time and use it to create time object.
        self.starttime = dt.datetime.strptime(stringtime, '%Y%m%d%H%M%S')

        # Because timestamp isn't stored in flux file yet...
        # we need to read it from the filename
        fnpart = os.path.split(self.filename)[-1]
        self['Time'] = parse_ftime(fnpart)

    def get_energy_index(self, targ_en):
        '''return index of nearest grid point in energy

        Arguments
        ---------
        targ_en : float
            Target energy (in keV)

        Returns
        -------
        en_idx : integer
            Index of nearest energy in RAM grid
        en_val : float
            Value of nearest energy in RAM grid
        '''
        diffs = np.abs(self['EnergyGrid'] - targ_en)
        en_idx = diffs.argmin()
        en_val = self['EnergyGrid'][en_idx]
        return en_idx, en_val

    def add_flux_plot(self, species='Hydrogen', energy=10, target=None,
                      maxz=1.0, minz=1e-15, add_cbar=True, loc=111,
                      add_tstamp=True, labelsize=15, title='auto',
                      drop_ghost=False):
        '''
        '''
        from matplotlib.colors import LogNorm
        from matplotlib.cm import get_cmap
        from matplotlib.pyplot import colorbar
        from matplotlib.ticker import LogLocator, LogFormatterMathtext

        fig, ax = splot.set_target(target, loc=loc, polar=True)
        var = f'Flux{species}'
        en_idx, en_val = self.get_energy_index(energy)
        alp_idx = 0
        alpha = self['PitchAngleGrid'][alp_idx]
        # Set title.
        if title == 'auto':
            title = f'{species} flux at {en_val:0.2f}keV, '
            title = title + r'$\alpha$='
            title = title + f'{alpha}' + '$^{\circ}$'

        # Set up grid centered on gridpoints.
        T_deg = 15*tb.bin_center_to_edges(self['MLTGrid'][:-1])
        T_rad = np.deg2rad(T_deg) - np.pi/2  # offset is for compat. with adjust_dialplot
        if not drop_ghost:
            R = tb.bin_center_to_edges(self['RadialGrid'])
            p = self[var][alp_idx, en_idx, :-1, :]
        else:
            R = tb.bin_center_to_edges(self['RadialGrid'][1:])
            p = self[var][alp_idx, en_idx, :-1, 1:]
        ax.grid(False)  # as of mpl 3.5 grid must be turned off before calling pcolormesh
        pcol = ax.pcolormesh(T_rad, R, p.T, norm=LogNorm(vmin=minz, vmax=maxz),
                             cmap=get_cmap('inferno'))
        ram._adjust_dialplot(ax, R, title=title, labelsize=15)
        if add_cbar:
            cbar = colorbar(pcol, pad=0.1, ticks=LogLocator(), ax=ax,
                            format=LogFormatterMathtext(), shrink=0.8)
            cbar.set_label('$s^{-1} sr^{-1} cm^{-3} keV^{-1}$')
        else:
            cbar = None
        if add_tstamp:
            tstamp = '{}'.format(self['Time'].isoformat()[:19])
            plt.figtext(0.65, 0.05, tstamp, fontsize=11)

        return fig, ax, pcol, cbar


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
        if (opt_ten >= tstart) and (opt_ten <= tend):
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


def parse_ftime(in_str):
    '''Parse date/time in flux file names and return datetime
    '''
    tparts = re.search('ram_flux_d(\d{4})(\d{2})(\d{2})_t(\d{2})(\d{2})(\d{2})',
                       in_str).groups()
    return dt.datetime(*[int(tp) for tp in tparts])


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
            # get all total pressures for plot
            plotlist = [key for key in pressf if key.startswith('tot')]
            # hold back Helium for summaries
            plotlist = [key for key in plotlist if not key.endswith('He')]
            # Convenience plotting routines in spacepy currently assume a log-transform
            # but this is not useful for anisotropy. When control for that is added to
            # spacepy we can do this more easily...
            # plotlist = [key for key in pressf if (key.startswith('tot')
            #             or key.startswith('ani'))]
            for pl in plotlist:
                if pl == 'tote':
                    # electron pressure is ~10% so shift range 1 order lower
                    minz = 1e-1
                    maxz = 1e2
                else:
                    minz = 1
                    maxz = 1e3
                # TODO: maybe make combined plots, look at axis limits, etc.
                title = '{}'.format(pressf[pl].attrs['label'])
                tstamp = '{}'.format(ftime.isoformat()[:19])
                fig, ax, cm, _ = pressf.add_pcol_press(var=pl, add_cbar=True, title=title,
                                                       minz=minz, maxz=maxz)
                outfn = os.path.join(options.outdir, fmain + f'_{pl}.png')
                plt.figtext(0.05, 0.95, tstamp, fontsize=11)
                plt.savefig(outfn, dpi=200)
                plt.close()


def plot_flux(options):
    fluxcands = sorted(glob.glob(os.path.join(options.rundir,
                                              'output_ram',
                                              'ram_flux*nc')))
    # loop over all files, if in expected time range make plot
    for flfn in fluxcands:
        fpath, fmain = os.path.split(flfn)
        fmain = fmain.split('.')[0]
        ftime = parse_ftime(fmain)
        opt_tst = spt.Ticktock(options.startTime).UTC[0]
        opt_ten = spt.Ticktock(options.endTime).UTC[0]
        if (ftime <= opt_ten) and (ftime >= opt_tst):
            fluxf = Flux2DFile(flfn)
            plotlist = ['Hydrogen', 'OxygenP1', 'Electron']
            splims = {'Hydrogen': [5e2, 1e7],
                      'HeliumP1': [5e1, 1e7],
                      'OxygenP1': [5e0, 5e6],
                      'Electron': [1e2, 1e7]}
            enlist = [1, 10, 30]  # plot energies in keV
            for pl, enval in it.product(plotlist, enlist):
                # TODO: maybe make combined plots, look at axis limits, etc.
                fpdict = {'species': pl, 'add_cbar': True, 'energy': enval,
                          'minz': splims[pl][0], 'maxz': splims[pl][1]}
                if pl == 'Electron':
                    fpdict['drop_ghost'] = True
                    if enval < 5:
                        fpdict['minz'] = 1e3
                        fpdict['maxz'] = 1e8
                fig, ax, cm, _ = fluxf.add_flux_plot(**fpdict)
                outfn = os.path.join(options.outdir, fmain + f'_{pl}_E{enval:0.2f}.png')
                plt.savefig(outfn, dpi=200)
                plt.close()


if __name__ == "__main__":
    parser = parserSetup()
    in_args = parser.parse_args()

    flagstatus = [False, False, False, False]
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
            flagstatus[2] = True
        flcand = os.path.join(in_args.rundir, 'output_ram', 'ram_flux*')
        if len(glob.glob(flcand)) > 0:
            plot_flux(in_args)
        else:
            flagstatus[3] = True
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
        if flagstatus[2]:
            errmsg += '\nNo 2D flux files found\n'
        parser.error(errmsg)
