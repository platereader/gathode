#!/usr/bin/env python
"""
This module implements a CLI for CATHODE.

Chronological life span Analysis Tool for High-throughput Optical
Density Experiments (CATHODE) command line interface (CLI) entry
point.
"""

# CATHODE  Chronological life span Analysis Tool
#          for High-throughput Optical Density Experiments
#
# Copyright (C) 2014 Nils Christian
#
# This file is part of CATHODE.
#
# CATHODE is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# CATHODE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with CATHODE.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
import argparse

from platereader.cls import Cls
from platereader.clsplot import survivalIntegralsToPdf
from platereader._version import __version__

def clsCommandlineInterfaceArgv(argv=['clscli']):
    commandline=' '.join(argv)
    executablename = os.path.basename(argv[0])
    argv.pop(0) # remove first argument (executable name), because we are explicitly passing it to parse_args
    parser = argparse.ArgumentParser(description='CATHODE-CLI (Chronological Life Span Analysis Tool for '+
                                     'High-throughput Optical Density Experiments): '+
                                     'A commandline interface for extracting the chronological life span '+
                                     'from multiple time series of optical density measurements.')
    parser.add_argument('--version', action='store_true', default=None, help='show version and exit')
    parser.add_argument('infiles', metavar='', nargs='*',
                        help='CATHODE file (.cat) or multiple GATHODE files (.gat)')
    parser.add_argument('--day', metavar='float', action='append', type=float, dest='days',
                        help='day of lag/growthrate measurement of the corresponding plate file')
    parser.add_argument('--viabilities', action='store', default=None, help='write viabilities to csv file')
    parser.add_argument('--siplots', action='store', default=None, help='write survival intergal figures to pdf')
    parser.add_argument('--labels', action='store', default=None, help='a file containing sampleid -> label and condition -> label mappings')

    args = parser.parse_args(argv)

    if args.version:
        print(executablename+' '+__version__)
        return 0

    if len(args.infiles) == 0:
        parser.print_help()
        print('\ninput file(s) missing')
        return -1
    elif args.infiles[0].endswith(".cat"):
        cls=Cls(serialisedFilename=args.infiles[0])
    else:
        cls=Cls(args.infiles,args.days)

    if args.viabilities:
        cls.survivalToCsv(args.viabilities)
    if args.siplots:
        labels={}
        if args.labels is not None:
            jsonfile=open(args.labels, 'r')
            labels = json.loads(jsonfile.read())
            jsonfile.close()
        survivalIntegralsToPdf(cls,pdfout=args.siplots,
                               sampleIdToLabel=labels['sampleIdToLabel'] if 'sampleIdToLabel' in labels else None,
                               conditionToLabel=labels['conditionToLabel'] if 'conditionToLabel' in labels else None)

    return 0


def clsCommandlineInterface():
    return clsCommandlineInterfaceArgv(list(sys.argv))

if __name__ == "__main__":
    sys.exit(clsCommandlineInterface())
