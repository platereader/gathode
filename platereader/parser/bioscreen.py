"""
This module implements a parser for files exported by the Bioscreen C platereader, used by GATHODE.

Growth Analysis Tool for High-throughput Optical Density Experiments
(GATHODE) Plate class.
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

import os.path
import re
import numpy
import platereader
from platereader.csvunicode import CsvFileUnicodeReader

def isPlateFormat(filename):
    """
    Calculate a score denoting the likelyhood of this file being a Bioscreen C exported plate.

    :param filename: name of the exported by the plate.
    :type filename: str

    :return: float -- score denoting likelyhood (range: [0,100], 100 means absolutely sure).
    """
    try:
        with CsvFileUnicodeReader(filename, delimiter='\t', quotechar='"') as odreader:
            firstrow=next(odreader)
            if re.search('^READER:\s+Bioscreen',firstrow[0]):
                return 100.
    except:
        pass
    return 0.

def parse(filename):
    """
    Read a Bioscreen C file

    :param odcsvfilename: filename of the Bioscreen export

    :return: numpy.array(float), list( numpy.array(float) ), list(str), str, numpy.array(float)
     -- time (in seconds), optical density readouts, sample ids, plate id, temperature

    :param seperator: split well ids on this seperator to distinguish between sample id and condition
    :type seperator: string
    """
    with CsvFileUnicodeReader(filename, delimiter='\t', quotechar='"') as odreader:
        # first line should contain 'Bioscreen'
        firstrow=next(odreader)
        if not re.search('^READER:\s+Bioscreen\sC',firstrow[0]):
            raise RuntimeError('could not identify Bioscreen file format')
        # next line is the plate id
        plateidrow=next(odreader)
        if not re.search('^TEST NAME:',plateidrow[0]):
            raise RuntimeError('could not identify plate id in Bioscreen file')
        plateId=re.sub("^TEST NAME:\s+","",plateidrow[0])
        # skip lines 3 to 6
        for i in range(3,7):
            row=next(odreader)
        # next column contains the ids
        sampleIds=next(odreader)
        timeId=sampleIds.pop(0) # remove time
        if timeId == 'TenthSec.':
            timemult=.1
        else:
            raise RuntimeError('could not identify time multiplier in Bioscreen file '+odcsvfilename+': '+timeId)
        rawOdList=[]
        for colid in sampleIds:
            # initialise empty list
            rawOdList.append([])
        time=[]
        for row in odreader:
            # divide by 10because these are TenthSec.
            colit=rawOdList.__iter__()
            t=row.pop(0)
            if not len(row):
                # this seems to be the last row that only contains one empty element (already 'pop'ped into t)
                continue
            time.append(timemult*float(t))
            for val in row:
                column=next(colit)
                column.append(val)

    # the time as numpy array
    time=numpy.array(time,dtype=float)

    # create a python list of numpy arrays of raw data
    rawOd=[]
    for l in rawOdList:
        rawOd.append(numpy.array(l,dtype=float))

    # add dummy conditions
    conditions=[]
    for i in range(len(rawOd)):
        conditions.append('')

    temperature=None
    wellids=None
    if len(rawOd) == 100:
        wellids=[str(i) for i in range(1,101)]
    elif len(rawOd) == 200:
        wellids=[str(i) for i in range(1,201)]
    return time, rawOd, sampleIds, conditions, plateId, temperature, wellids
