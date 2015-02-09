"""
This module implements some helper function for numpy.
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

import math
import numpy

def maskedArrayToMeanVar(marr,idcs=None,ddof=0,axis=None):
    """"
    Calls mean and var removing nans.

    Makes sure a None is returned if '--' would be returned by
    masked_array.
    """
    if idcs is None:
        idcs=numpy.isnan(marr)
    if numpy.all(idcs):
        return None, None
    themean=numpy.ma.masked_array(marr,idcs).mean(axis=axis)
    thevar=numpy.ma.masked_array(marr,idcs).var(ddof=ddof,axis=axis)

    # ensure that '--' is replaced with nan/None
    if axis is not None:
        themean[numpy.ma.getmaskarray(themean)]=None
        thevar[numpy.ma.getmaskarray(thevar)]=None            
    else:
        if type(themean) == numpy.ma.core.MaskedConstant:
            themean=None
        if type(thevar) == numpy.ma.core.MaskedConstant:
            thevar=None
        if themean is not None and math.isnan(themean):
            themean=None
        if thevar is not None and math.isnan(thevar):
            thevar=None

    return themean, thevar

def notNanAndGreaterEqual(vec,val):
    if numpy.any(numpy.isnan(vec)):
        # FIXME we make a copy here!
        veccopy = numpy.copy(vec)
        veccopy[numpy.isnan(vec)] = val - 42.
        return veccopy >= val
    return vec >= val

def notNanAndLess(vec,val):
    if numpy.any(numpy.isnan(vec)):
        # FIXME we make a copy here!
        veccopy = numpy.copy(vec)
        veccopy[numpy.isnan(vec)] = val + 42.
        return veccopy < val
    return vec < val

def nonNanSqrt(vec):
    s = numpy.empty(vec.shape)
    s[~numpy.isnan(vec)] = numpy.sqrt(vec[~numpy.isnan(vec)])
    s[numpy.isnan(vec)] = numpy.nan

    return s

def nonNanNonZeroDivide(vecnom,vecden):
    # NOTE this is faster than:
    # s = ma.divide(vecnom,vecden)
    # return s.filled(numpy.nan)
    s = numpy.empty(vecnom.shape)
    idcsok = numpy.logical_and(~numpy.isnan(vecnom), numpy.logical_and(vecden != 0., ~numpy.isnan(vecden)))
    s[idcsok] = numpy.divide(vecnom[idcsok],vecden[idcsok])
    s[~idcsok] = numpy.nan
    return s
