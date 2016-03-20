"""
Helper function for parsers of platereader file formats, used by GATHODE.

Growth Analysis Tool for High-throughput Optical Density Experiments
(GATHODE) parser helper functions.
"""

# GATHODE  Growth Analysis Tool
#          for High-throughput Optical Density Experiments
#
# Copyright (C) 2015 Nils Christian
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
import re
import numpy
import platereader
import datetime

def appendElementsToListOfLists(rawOdList,row):
    """Helper function for appending the values of a row to a list of lists"""
    colit=rawOdList.__iter__()
    for val in row:
        column=next(colit)
        column.append(val)

def splitSampleCondition(sampleCondition,separator='_'):
    sampleIds=[]
    conditions=[]
    for idCondition in sampleCondition:
        if separator is not None:
            (sampleid,sep,condition)=idCondition.rpartition(separator)
        else:
            sampleid=idCondition
            condition=''
        sampleIds.append(sampleid)
        conditions.append(condition)

    return sampleIds, conditions
