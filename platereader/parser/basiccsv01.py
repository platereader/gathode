"""
This module implements a parser for files in a simple csv format, used by GATHODE.

Growth Analysis Tool for High-throughput Optical Density Experiments
(GATHODE) csv format parser.
"""

# GATHODE  Growth Analysis Tool
#          for High-throughput Optical Density Experiments
#
# Copyright (C) 2015-2016 Nils Christian
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
import platereader.parser.utils
from platereader.csvunicode import CsvFileUnicodeReader

def isPlateFormat(filename):
    """
    Calculate a score denoting the likelyhood of this file being a GATHODE plate in csv format.

    :param filename: filename of the exported plate.
    :type filename: str

    :return: float -- score denoting likelyhood (range: [0,100], 100 means absolutely sure).
    """
    try:
        dummy = _parseGathodeCsv01AnyCsv(filename)
        return 30
    except:
        pass
    return 0.

def parse(filename):
    """
    Read a GATHODE plate in csv format.

    :param filename: filename of the the csv file

    :return: numpy.array(float), list( numpy.array(float) ), list(str), str, numpy.array(float)
     -- time (in seconds), optical density readouts, sample ids, plate id, temperature
    """
    return _parseGathodeCsv01AnyCsv(filename)

def _parseGathodeCsv01AnyCsv(filename,**csvkwargs):
    csvkwargsCombinations=[
        { 'encoding': 'utf-8', 'dialect': 'excel' },
        { 'encoding': 'utf-8', 'dialect': 'excel-tab' },
        { 'encoding': 'iso-8859-1', 'dialect': 'excel' },
        { 'encoding': 'iso-8859-1', 'dialect': 'excel-tab' },
        { 'encoding': 'utf-16', 'dialect': 'excel' },
        { 'encoding': 'utf-16', 'dialect': 'excel-tab' },
        ]
    for csvkwargs in csvkwargsCombinations:
        try:
            time, rawOd, sampleIds, conditions, plateId, temperature, wellids = _parseGathodeCsv01Helper(filename,**csvkwargs)
            if len(time) > 5 and len(rawOd) > 5:
                return time, rawOd, sampleIds, conditions, plateId, temperature, wellids
        except Exception as err:
            pass
    # if we got here this file is not in the correct format
    raise RuntimeError('could not identify GATHODE csv file format')

def _parseGathodeCsv01Helper(filename,**csvkwargs):
    """ Main parsing function for csv output """
    outbase, outext = os.path.splitext(os.path.basename(filename))
    plateId = outbase

    with CsvFileUnicodeReader(filename,**csvkwargs) as odreader:
        # first line should contain 'time'
        headerrow=next(odreader)
        if headerrow[0].lower() != 'time':
            raise RuntimeError('could not identify GATHODE csv file format, expected "time" in first cell, got "'+str(headerrow[0].lower())+'"')
        rawOdList=[]
        for colid in headerrow:
            # initialise empty list
            rawOdList.append([])
        time=[]
        for row in odreader:
            time.append(row.pop(0))
            platereader.parser.utils.appendElementsToListOfLists(rawOdList,row)

    # remove 'time' from header
    headerrow.pop(0)

    # if last array is empty it should be removed
    if len(rawOdList[-1]) == 0:
        rawOdList.pop()

    # the time as numpy array
    time=numpy.array(time,dtype=float)

    # create a python list of numpy arrays of raw data
    rawOd=[]
    for l in rawOdList:
        rawOd.append(numpy.array(l,dtype=float))

    sampleIds = platereader.parser.utils.replaceIntegerSampleIdsWithSortableStrings(headerrow)

    sampleIds, conditions = platereader.parser.utils.splitSampleCondition(sampleIds,separator='_')

    temperature=None

    return time, rawOd, sampleIds, conditions, plateId, temperature, platereader.plate.Plate.guessWellIds(len(rawOd))
