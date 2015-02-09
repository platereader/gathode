#! /usr/bin/env python
"""
This module implements a CLI for GATHODE.

Growth Analysis Tool for High-throughput Optical Density Experiments
(GATHODE) command line interface (CLI) entry point.
"""

# GATHODE  Growth Analysis Tool
#          for High-throughput Optical Density Experiments
#
# Copyright (C) 2014 Nils Christian
#
# This file is part of GATHODE.
#
# GATHODE is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# GATHODE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with GATHODE.  If not, see <http://www.gnu.org/licenses/>.


import sys
import os.path
import argparse

from platereader.plate import Plate
from platereader.odplot import plotFullOdPlate
from platereader._version import __version__

def odCommandlineInterfaceArgv(argv=['odcli']):
    commandline=' '.join(argv)
    executablename = os.path.basename(argv[0])
    argv.pop(0) # remove first argument (executable name), because we are explicitly passing it to parse_args
    parser = argparse.ArgumentParser(description='GATHODE-CLI (Growth Analysis Tool for High-throughput '+
                                     'Optical Density Experiments): '+
                                     'A commandline interface for analysing time series of optical density '+
                                     'measurements that were recorded with the help of a plate reader.')
    parser.add_argument('--version', action='store_true', default=None, help='show version and exit')
    parser.add_argument('infile', metavar='file.gat|export.asc', action='store', default=None, nargs='?',
                        help='file to be loaded (gat-file or TECAN ASCII format)')
    parser.add_argument('--pdf', action='store', default=None, help='write figures to pdf')
    parser.add_argument('--csvout', action='store', default=None,   help='write growth data to csv')
    parser.add_argument('--s', '-s', action='store', type=float, default=0.01,
                        help='smoothing factor passed to UnivariateSpline')
    parser.add_argument('--k', '-k', action='store', type=int, default=5,
                        help='degree of smoothing spline')
    parser.add_argument('--hdlin', action='store', type=float, default=1.,
                        help='high density correction linear term')
    parser.add_argument('--hdquad', action='store', type=float, default=0.,
                        help='high density correction quadratic term')
    parser.add_argument('--hdcub', action='store', type=float, default=0.,
                        help='high density correction cubic term')
    parser.add_argument('--logodcutoff', action='store', type=float, default=-5,
                        help='cutoff applied to ln(OD) when calculating various observables')
    parser.add_argument('--maxgrowthlowertimecutoff', action='store', type=float, default=1.,
                        help='lower cutoff applied to t when calculating maximal growth rate')
    parser.add_argument('--maxgrowthuppertimecutoff', action='store', type=float, default=None,
                        help='upper cutoff applied to t when calculating maximal growth rate')
    parser.add_argument('--lagatlogodequals', action='store', type=float, default=-5,
                        help='the lag will be defined as the intersection point of the linear equation '
                        +'of the maximal growth with the given ln(OD) value')
    parser.add_argument('--onlyAveraged', action='store_true', default=True,
                        help='show only the averaged replicates, not each well individually')
    parser.add_argument('--gat', action='store', help='save as gat-file')
    args = parser.parse_args(argv)

    if args.version:
        print(executablename+' '+__version__)
        return 0
    if args.infile is not None:
        plate=Plate(filename=args.infile)
    else:
        parser.print_help()
        print('\ninput file missing')
        return -1
    if plate.readfileformat != 'gat':
        plate.setHighDensityCorrectionLinear(args.hdlin)
        plate.setHighDensityCorrectionQuadratic(args.hdquad)
        plate.setHighDensityCorrectionCubic(args.hdcub)
        plate.setLogOdCutoff(args.logodcutoff)
        plate.setSmoothingS(args.s)
        plate.setSmoothingK(args.k)
        plate.setMaxGrowthLowerTimeCutoff(args.maxgrowthlowertimecutoff)
        plate.setMaxGrowthUpperTimeCutoff(args.maxgrowthuppertimecutoff)
        plate.setLagAtLogOdEquals(args.lagatlogodequals)

    if args.csvout is not None:
        plate.growthParametersToCsv(args.csvout)
    elif args.gat is not None:
        plate.save(args.gat)
    elif args.pdf is not None:
        plotFullOdPlate(plate,pdfout=args.pdf,creator=commandline,showReplicateGroups=args.onlyAveraged)
    else:
        print("don't know what to do")
        return -1
    return 0

def odCommandlineInterface():
    return odCommandlineInterfaceArgv(list(sys.argv))

if __name__ == "__main__":
    sys.exit(odCommandlineInterface())
