# -*- coding: utf-8 -*-
# NOTE the coding is needed because we have to strip the ugly °C from the input

"""
This module implements a parser for files exported by the TECAN platereader, used by GATHODE.

Growth Analysis Tool for High-throughput Optical Density Experiments
(GATHODE) TECAN platereader format parser.
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
    Calculate a score denoting the likelyhood of this file being a TECAN exported plate.

    :param filename: name of the exported by the plate.
    :type filename: str

    :return: float -- score denoting likelyhood (range: [0,100], 100 means absolutely sure).
    """
    encodings=["utf-8", "utf-16", "iso-8859-1"]
    for encoding in encodings:
        try:
            return _isTecanFormatHelper(filename,encoding=encoding)
        except:
            pass
    # didn't manage to parse, so this is not a TECAN plate
    return 0.

def _isTecanFormatHelper(filename,encoding="utf-8"):
    with CsvFileUnicodeReader(filename,encoding=encoding,delimiter='\t',quotechar='"') as odreader:
        subfromscore=0
        sampleIdsNCondition=next(odreader)
        # first line may contain ids (first column would be empty)
        if sampleIdsNCondition[0] != '':
            # penalty if first line does not contain ids
            subfromscore=20

        # we expect 96 or 384 wells, plus 2 or 3 extra columns
        if (len(sampleIdsNCondition) != 98 and len(sampleIdsNCondition) != 99
            and len(sampleIdsNCondition) != 386 and len(sampleIdsNCondition) != 387):
            return 0.

        # first entry should be numeric + 's', second numeric + possibly unit of temperature
        # FIXME isn' there a better way to check whether something is numeric?
        secondrow=next(odreader)
        if (re.search('^[+-]?(\d+\.\d+|\d+\.|\.\d+|\d+)([eE][+-]?\d+)?\s*s$',secondrow[0])
            and re.search('^[+-]?(\d+\.\d+|\d+\.|\.\d+|\d+)([eE][+-]?\d+)?',secondrow[1])):
            return 50.-subfromscore
    return 0.

def parse(odcsvfilename,seperator='_'):
    """
    Read a TECAN file.

    :param odcsvfilename: filename of the Tecan export

    :return: numpy.array(float), list( numpy.array(float) ), list(str), str, numpy.array(float), list(str)
     -- time (in seconds), optical density readouts, sample ids, plate id, temperature, wellids

    :param seperator: split well ids on this seperator to distinguish between sample id and condition
    :type seperator: string
    """
    encodings=["utf-8", "utf-16", "iso-8859-1"]
    for encoding in encodings:
        try:
            return _parseTecanCsvExportHelper(odcsvfilename,seperator,encoding=encoding)
        except (UnicodeError, UnicodeDecodeError) as e:
            pass
    raise RuntimeError('Error parsing csv file, tried the encodings '+' '.join(encodings))

def _handleDataRow(rawOdList,row):
    """Helper function for _parseTecanCsvExportHelper"""
    colit=rawOdList.__iter__()
    for val in row:
        column=next(colit)
        column.append(val)      

def _parseTecanCsvExportHelper(odcsvfilename,seperator='_',encoding="utf-8"):
    """Helper function for _parseTecanCsvExport"""

    plateId=os.path.basename(odcsvfilename)

    with CsvFileUnicodeReader(odcsvfilename,encoding=encoding,delimiter='\t',quotechar='"') as odreader:
        sampleIdsNCondition=next(odreader)
        # initialise empty lists according to number of columns
        rawOdList=[[] for idCondition in sampleIdsNCondition]
        if sampleIdsNCondition[0] == '':
            # this was an export including ids in the first line
            # first column is the time (but the column is not marked as such)
            sampleIdsNCondition[0]="time"
            nonSampleIndices=1
            if sampleIdsNCondition[1] == "": # second column is the temperature (but the column is not marked as such)
                sampleIdsNCondition[1]="temperature"
                nonSampleIndices=2
        else:
            # this was an export not containing ids
            if len(sampleIdsNCondition) == 99 or len(sampleIdsNCondition) == 387:
                nonSampleIndices=2
                # this should contain a temperature column
                if not re.search("\s* \xb0C$",sampleIdsNCondition[1]):
                    raise RuntimeError('Error parsing csv file, first line should contain a temperature in second column')
            elif len(sampleIdsNCondition) == 98 or len(sampleIdsNCondition) == 386:
                nonSampleIndices=1
            else:
                raise RuntimeError('Error parsing csv file, does not seem to correspond to 96 or 384 well plate format')

            # handle first row, as it contains data
            _handleDataRow(rawOdList,sampleIdsNCondition)

            # assign numbers to labels
            sampleIdsNCondition=["{:0>3d}".format(i) for i in range(1-nonSampleIndices,len(sampleIdsNCondition)+1-nonSampleIndices)]
            # make sure labels end up as ids, not as condition
            seperator=None

        # read data
        for row in odreader:
            _handleDataRow(rawOdList,row)

    # NOTE remove last item because it is empty
    sampleIdsNCondition.pop()
    rawOdList.pop()

    # sanitise time (remove unit 's')
    newtime=[]
    for t in rawOdList[0]:
        newtime.append(re.sub("s$", "", t))
    rawOdList[0]=newtime

    # the time as numpy array
    time=numpy.array(rawOdList[0],dtype=float)
    # remove time from rawOdList
    rawOdList.pop(0)
    sampleIdsNCondition.pop(0)

    # the temperature
    temperature=None
    if nonSampleIndices >= 2: # second column is the temperature (but the column is not marked as such)
        # remove temperature from rawOdList
        oldtemp=rawOdList.pop(0)
        sampleIdsNCondition.pop(0)
        # sanitise temperature (remove unit '°C')
        newtemp=[]
        for t in oldtemp:
            newtemp.append(re.sub("\s* \xb0C$", "", t))
        temperature=numpy.array(newtemp,dtype=float)

    # create a python list of numpy arrays of raw data
    rawOd=[]
    for l in rawOdList:
        rawOd.append(numpy.array(l,dtype=float))

    # 
    sampleIds=[]
    conditions=[]
    for idCondition in sampleIdsNCondition:
        if seperator is not None:
            (sampleid,sep,condition)=idCondition.rpartition(seperator)
        else:
            sampleid=idCondition
            condition=''
        sampleIds.append(sampleid)
        conditions.append(condition)

    return time, rawOd, sampleIds, conditions, plateId, temperature, platereader.plate.Plate.guessWellIds(len(rawOd))
