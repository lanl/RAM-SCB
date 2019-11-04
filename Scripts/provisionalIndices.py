#!/usr/bin/env python

import os, re, copy, bisect, ftplib
from optparse import OptionParser
import datetime as dt
import numpy as np
import spacepy.time as spt
#Py3k/Py2 imports
try:
    import urllib2 as url
except:
    import urllib.request as url
try:
    import cStringIO as io
except:
    import io
try:
    from itertools import izip_longest as zip_longest
except:
    from itertools import zip_longest

defaults = {'rundir': os.path.abspath('../rtRun')}

def parserSetup():
    # Define a command-line option parser and add the options we need
    parser = OptionParser(  usage="%prog [options] PARAMfile",\
                            version="%prog Version 1.00 (June 28, 2017)"  )

    parser.add_option("-d", "--rundir",      dest="rundir",
                            help="Full path to RAM run directory. Default is evaluated as '../rtRun' ")


    return parser


def floatkp(kpstring):
    n = int(kpstring[0])
    dlt = 0
    if len(kpstring)>1:
        kpsym = kpstring[1]
        if kpsym == '-':
            dlt = -0.3
        elif kpsym.lower() == 'o':
            dlt = 0
        elif kpsym == '+':
            dlt = 0.3
        else:
            raise Exception(kpsym + ' is not a recognized symbol')
    return n+dlt

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    # recipe from itertools documentation
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)

def Ap2Kp(ap):
    '''estimate Kp index from Ap value'''
    aparr = np.array([0,2,3,4,5,6,7,9,12,15,18,22,27,32,39,48,56,67,80,94,111,132,154,179,207,236,300,400])
    kparr = np.array([0,0.3,0.7,1,1.3,1.7,2,2.3,2.7,3,3.3,3.7,4,4.3,4.7,5,5.3,5.7,6,6.3,6.7,7,7.3,7.7,8,8.3,8.7,9])
    estKp = kparr[bisect.bisect_left(aparr, ap)]
    return estKp

def validLine(instr):
    out = True if (len(instr.split()[0])>5) else False #first word must be a complete date (yymmdd)
    return out


if __name__=='__main__':
    parser = parserSetup()
    # Parse the args that were (potentially) given to us on the command line
    (options, in_args) = parser.parse_args()

    #roll defaults into options dict and do any checking
    valid_opt = 0
    for key in options.__dict__.keys():
        if options.__dict__[key] != None:
            valid_opt += 1
        try:
            if options.__dict__[key] == None:
                options.__dict__[key] = defaults[key]
        except KeyError:
            pass

    #Now we're done with the option parser setup, start the data managament
    Ap2Kp_vec = np.vectorize(Ap2Kp)
    
    #provisional kp, last month
    kp_lastmm = 'http://www-app3.gfz-potsdam.de/kp_index/pqlyymm.tab'
    #provisional 'quicklook' Kp, current month
    kp_thismm = 'http://www-app3.gfz-potsdam.de/kp_index/qlyymm.tab'
    
    #read data from last, this, (concat), then get provisional 10.7 (somehow)
    #write all out to file
    response = url.urlopen(kp_lastmm)
    indata = response.read()
    try:
        indata = indata.decode('latin-1') #python3 encodes as bytes
    except:
        pass
    
    #now extract the data
    indata = filter(bool, indata.split('\n')) #filtering removes empty strings
    kp_lastmonth = [item for item in indata if validLine(item)]
    data = {}
    data['DateTime'] = np.array([dt.datetime.strptime(dd.split()[0],'%y%m%d') for dd in kp_lastmonth])
    data['Kp'] = np.empty([len(data['DateTime']), 8])
    for idx, line in enumerate(kp_lastmonth):
        data['Kp'][idx,...] = [floatkp(val) for val in line.split()[1:9]]
    offset = copy.copy(idx)+1
    
    response = url.urlopen(kp_thismm)
    indata = response.read()
    try:
        indata = indata.decode('latin-1') #python3 encodes as bytes
    except:
        pass
    
    #now extract the data
    indata = filter(bool, indata.split('\n')) #filtering removes empty strings
    kp_thismonth = [item for item in indata if validLine(item)]
    data['DateTime'] = np.hstack([data['DateTime'], np.asarray([dt.datetime.strptime(dd.split()[0],'%y%m%d') for dd in kp_thismonth])])
    data['Kp'] = np.vstack([data['Kp'], np.empty([len(kp_thismonth), 8])])
    for idx, line in enumerate(kp_thismonth):
        kp_vals = [floatkp(val) for val in line.split()[1:9]]
        if len(kp_vals) == 8:
            data['Kp'][idx+offset,...] = kp_vals
        else:
            data['Kp'][idx+offset, :len(kp_vals)] = kp_vals
            data['Kp'][idx+offset, len(kp_vals):] = kp_vals[-1] #fill rest of day with persistence
    
    #Get the F10.7 data through the present day, for all dates in data['DateTime']
    julians = spt.Ticktock(data['DateTime']).JD + 0.5 #add a half to get noontime flux estimate
    data['F10.7'] = np.empty(len(data['DateTime']))

    #download updated F10.7 from ftp://ftp.geolab.nrcan.gc.ca/data/solar_flux/daily_flux_values/fluxtable.txt
    #only has data from 2004-10-28 (prior data are in different files)
    io_store = io.StringIO()
    ftp = ftplib.FTP('ftp.seismo.nrcan.gc.ca')
    ftp.login()
    ftp.cwd('spaceweather/solar_flux/daily_flux_values')
    ftp.retrlines('RETR {0}'.format('fluxtable.txt'), io_store.write)
    io_store.seek(0)
    conts = io_store.read()
    
    #get first line of data (needs regex as no linebreaks & odd spacing
    #then split up, read dates, interpolate to requested JDs
    char1idx = re.search('\d', conts).start()
    conts = conts[char1idx:].split()
    jds = np.array(conts[2::7]).astype(float)
    AU1flux = np.array(conts[5::7]).astype(float)
    data['F10.7'][...] = np.interp(julians, jds, AU1flux)


    ## FORECAST DATA
    #############################################
    #Now get 45-day Air Force forecast to add...#
    #############################################
    rt_indices = os.path.join(options.rundir, 'RamIndices.txt')
    
    #45 day forecast has daily Ap and F10.7 flux
    latest107file = 'http://services.swpc.noaa.gov/text/45-day-ap-forecast.txt'
    
    response = url.urlopen(latest107file)
    indata = response.read()
    try:
        indata = indata.decode('latin-1') #python3 encodes as bytes
    except:
        pass
    
    #now extract the data
    indata = filter(bool, indata.split('\n')) #filtering removes empty strings
    header, data107, dataAP = [], [], []
    ap, f107 = False, False
    for line in indata:
        if line[0] in (':', "#"):
            header.append(line.strip())
        elif 'AP' in line.upper():
            #AP forecast block, collect data
            ap = True
            continue
        elif '10.7' in line.upper():
            ap = False
            f107 = True
            continue
        elif 'FORECASTER' in line.upper():
            break #last block of file has forecaster info, we're done
        else:
            savedata = dataAP if ap else data107
            savedata.append(line)
    
    data107 = ' '.join(data107)
    data107 = data107.split()
    dataAP  = ' '.join(dataAP)
    dataAP  = dataAP.split()
    
    #collect array of F10.7 values with dates
    data107 = np.array(list(grouper(data107,2))).astype(object)
    data107[:,0] = [dt.datetime.strptime(nn,'%d%b%y') for nn in data107[:,0]]
    data107[:,1] = [float(nn) for nn in data107[:,1]]
    #collect array of AP values with dates
    dataAP = np.array(list(grouper(dataAP,2))).astype(object)[:,1]
    dataAP[...] = [float(nn) for nn in dataAP]
    dataKp = Ap2Kp_vec(dataAP)
    dataKp = np.reshape(np.repeat(dataKp,8),[len(dataKp),8])
    #output
    #get last date in provisional/quicklook array
    lastdate = data['DateTime'][-1]
    keepinds = data107[:,0] > lastdate
    data['DateTime'] = np.hstack([data['DateTime'], data107[keepinds,0]])
    data['Kp'] = np.vstack([data['Kp'], dataKp[keepinds,:]])
    data['F10.7'] = np.hstack([data['F10.7'], data107[keepinds,1]])

    #now write to RamIndices-formatted file
    with open(rt_indices, 'w') as fh:
        fh.write('Provisional RamIndices (Date, Kp, F10.7)\n') #RAM requires 1 header line to read correctly
        for idx, od in enumerate(data['DateTime']):
            of = data['F10.7'][idx]
            usekp = data['Kp'][idx,:].tolist()
            printme = ''.join(['{0: 2.1f}'.format(k) for k in usekp])
            fh.write('{0}{1} {2:05.1f}\n'.format(od.strftime('%Y%m%d'), printme, of))
