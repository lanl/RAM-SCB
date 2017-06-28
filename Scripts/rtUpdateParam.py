#!/usr/bin/env python

import re, os
from optparse import OptionParser
try:
    import StringIO as io
except:
    import io #python 3 compat.

defaults = {'extraTime': 3600, #add one hour
            'outFile': None
            }

def parserSetup():
    # Define a command-line option parser and add the options we need
    parser = OptionParser(  usage="%prog [options] PARAMfile",\
                            version="%prog Version 0.99 (June 28, 2017)"  )

    parser.add_option("-t", "--extraTime",      dest="extraTime",
                            help="Number of extra seconds to run in restarted simulation")

    parser.add_option("-o", "--outFile",      dest="outFile",
                            help="Optional output filename. Default is to overwrite input file.")


    return parser

def replaceSTARTTIME(fn, options):
    '''Scan file for #STARTTIME block and replace entire block with #RESTART
    '''
    expr = re.compile('#STARTTIME')
    startblock = False
    store = io.StringIO()
    with open(fn, 'r') as fh:
        lines = fh.readlines()
        lines = [line.rstrip() for line in lines] #ensure linux newline by replacing whatever is there
        ##replace simulation end time with n + additional time
        for line in lines:
            searchyn = re.match(expr, line)
            if searchyn: #target line
                startblock = True
                store.write('#RESTART\n\n')
            elif line and line[0]!='#' and startblock:
                pass #skip storing these lines
            elif ((not line) or (line[0]=='#')) and startblock:
                startblock = False #exiting STARTTIME block
            else:
                store.write(line+'\n')
    store.seek(0)
    with open(options.outFile,'w') as fh:
        fh.writelines(store.readlines())

def updateSTOP(options):
    '''Scan file for #STOP block and update end time
    '''
    expr = re.compile('(^[0-9]+)*([" "\t]+)(tSimulationMax)')
    store = io.StringIO()
    with open(options.outFile, 'r') as fh:
        lines = fh.readlines()
        ##replace simulation end time with n + additional time
        for line in lines:
            searchyn = re.match(expr, line)
            if searchyn: #target line
                parts = searchyn.groups()
                fmtr = (int(parts[0])+int(options.extraTime), parts[1], parts[2]) 
                line = '{}{}{}\n'.format(*fmtr)
            store.write(line)
    store.seek(0)
    with open(options.outFile,'w') as fh:
        fh.writelines(store.readlines())

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

    #check for any input arguments, if not, print help
    if (not in_args) or (not os.path.isfile(in_args[0])):
        parser.print_help()
        exit()

    paramfile = in_args[0] #name of file to work with
    if options.outFile is None:
        options.outFile = paramfile

    replaceSTARTTIME(paramfile, options)
    updateSTOP(options)
