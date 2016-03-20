"""
This module implements the :py:class:`Replicate` class for GATHODE.

Growth Analysis Tool for High-throughput Optical Density Experiments
(GATHODE) Replicate class.
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

import warnings
import math
import numpy
import scipy.stats
from scipy.interpolate import UnivariateSpline
try:
    w = scipy.optimize.OptimizeWarning()
    from scipy.optimize import OptimizeWarning as OptimizeWarning
except AttributeError:
    # this is an old scipy without OptimizeWarning, create a dummy class
    class OptimizeWarning(UserWarning):
        pass

from platereader.numpytools import maskedArrayToMeanVar, notNanAndGreaterEqual, notNanAndLess
from platereader.statusmessage import StatusMessage, Severity

class Replicate(object):
    """
    A time series of optical density measurements for a single well or a group of replicates.
    """

    _isPurePlateParameter={
        'smoothingK': True,
        'smoothingS': True,
        'hdCorrectionLinear': True,
        'hdCorrectionQuadratic': True,
        'hdCorrectionCubic': True,
        'slidingWindowSize': True,
        'lagAtLogOdEquals': True,
        'logOdCutoff': True,
        }

    # for an explanation of this dict see _clearMemoised and _parametersUpdated
    _memoisedDontClear={
        'allowGrowthyieldSlopeNStderrAwayFromZero': set(['derivative',
                                                         'expFitsMu',
                                                         'expFitsOd0Mu',
                                                         'logOd',
                                                         'logOdSmoothed',
                                                         'od',
                                                         'odVar',
                                                         'rawOd',
                                                         'rawOdVar',
                                                         'smoothedOd',
                                                         'smoothedOdDerivative']),
        'allowMaxGrowthrateAtLowerCutoff': set(['derivative',
                                                'expFitsMu',
                                                'expFitsOd0Mu',
                                                'logOd',
                                                'logOdSmoothed',
                                                'od',
                                                'odVar',
                                                'rawOd',
                                                'rawOdVar',
                                                'smoothedOd',
                                                'smoothedOdDerivative']),
        'backgroundIndex': set(['rawOd', 'rawOdVar']),
        'hdCorrectionCubic': set(['rawOd', 'rawOdVar']),
        'hdCorrectionLinear': set(['rawOd', 'rawOdVar']),
        'hdCorrectionQuadratic': set(['rawOd', 'rawOdVar']),
        'lagAtLogOdEquals': set(['derivative',
                                 'expFitsMu',
                                 'expFitsOd0Mu',
                                 'logOd',
                                 'logOdSmoothed',
                                 'od',
                                 'odVar',
                                 'rawOd',
                                 'rawOdVar',
                                 'smoothedOd',
                                 'smoothedOdDerivative']),
        'logOdCutoff': set(['derivative',
                            'expFitsMu',
                            'expFitsOd0Mu',
                            'logOd',
                            'logOdSmoothed',
                            'od',
                            'odVar',
                            'rawOd',
                            'rawOdVar',
                            'smoothedOd',
                            'smoothedOdDerivative']),
        'maxGrowthLowerTimeCutoff': set(['derivative',
                                         'expFitsMu',
                                         'expFitsOd0Mu',
                                         'logOd',
                                         'logOdSmoothed',
                                         'od',
                                         'odVar',
                                         'rawOd',
                                         'rawOdVar',
                                         'smoothedOd',
                                         'smoothedOdDerivative']),
        'maxGrowthUpperTimeCutoff': set(['derivative',
                                         'expFitsMu',
                                         'expFitsOd0Mu',
                                         'logOd',
                                         'logOdSmoothed',
                                         'od',
                                         'odVar',
                                         'rawOd',
                                         'rawOdVar',
                                         'smoothedOd',
                                         'smoothedOdDerivative']),
        'slidingWindowSize': set(['derivative',
                                  'logOd',
                                  'logOdSmoothed',
                                  'od',
                                  'odVar',
                                  'rawOd',
                                  'rawOdVar',
                                  'smoothedOd',
                                  'smoothedOdDerivative']),
        'smoothingK': set(['derivative',
                           'expFitsMu',
                           'expFitsOd0Mu',
                           'logOd',
                           'od',
                           'odVar',
                           'rawOd',
                           'rawOdVar']),
        'smoothingS': set(['derivative',
                           'expFitsMu',
                           'expFitsOd0Mu',
                           'logOd',
                           'od',
                           'odVar',
                           'rawOd',
                           'rawOdVar'])
        }

    def __init__(self,parentPlate=None,wellIndices=None,sampleid=None,condition=None,wellids=None,
                 activeWellIndices=None,isReplicateGroup=False,
                 _unpickled=None,_serialiseFormat=None):
        """
        Constructor.

        :param parentPlate: Object of the plate this well/replicate group is part of.
        :type parentPlate: Plate
        :param wellIndices: Indices into the parent-plate's array of well objects.
        :type wellIndices: list(int)
        :param sampleid: Identifier of this Replicate.
        :type sampleid: str
        :param condition: Growth condition of this Replicate.
        :type condition: str
        :param wellids: Identifiers denoting the location(s) of this replicate's wells within the plate.
        :type wellids: list(str)
        :param activeWellIndices: Indices into the vector of child wells of this replicate group (from interval [0,len(wellIndices)[)
        :type activeWellIndices: list(int)
        :param isReplicateGroup: Whether this is a replicate group (not a single well)
        :param isReplicateGroup: bool

        .. note::
            The following two parameters should only be used when
            implementing a deserialiser.

        :param _unpickled: dictionary of serialised Replicate object; used by plate deserialisation code.
        :param _unpickled: dict
        :param _serialiseFormat: format of serialisation; used by plate deserialisation code.
        :param _serialiseFormat: str
        """
        self.sampleid=None
        self.condition=None
        self.wellids=None
        self._inheritableParameters={}
        self._inheritableParameters['maxGrowthLowerTimeCutoff']=None
        self._inheritableParameters['maxGrowthUpperTimeCutoff']=None
        self._inheritableParameters['allowMaxGrowthrateAtLowerCutoff']=None
        self._inheritableParameters['allowGrowthyieldSlopeNStderrAwayFromZero']=None
        self._wellIndices=None
        self._activeWellIndices=None
        self._backgroundIndex=None
        self._memoised={}
        self.time=None         # timepoints (numpy array, a reference to the parentPlate for quick access)
        self.timeunit=None     # the unit of the time (s, h, ...)
        self.parentPlate=None
        self._replicateGroupParent=None # if and only if this replicate is part of a replicate group _replicateGroupParent is not None
        self._isReplicateGroup=isReplicateGroup # whether this is a replicate group; NOTE checking _replicateGroupParent is not enough
        self.background=None   # background (Replicate object)
        if _unpickled is not None:
            self._deserialise(_unpickled,parentPlate,_serialiseFormat,isReplicateGroup)
        else:
            self.time=parentPlate.time         # timepoints (numpy array, a reference to the parentPlate for quick access)
            self.timeunit=parentPlate.timeunit # the unit of the time (s, h, ...)
            self.parentPlate=parentPlate
            self._setCondition(condition)
            self._setWellIds(wellids)
            self._setSampleId(sampleid,wellIndices)
            self._setWellIndices(wellIndices,activeWellIndices)

    def _deserialise(self,unpickled,parent,serialiseFormat,isReplicateGroup):
        if serialiseFormat is None:
            raise RuntimeError("no serialisation format (version) defined")

        self.time=parent.time
        self.timeunit=parent.timeunit
        self.parentPlate=parent
        self.sampleid=unpickled['sampleId']
        self.condition=unpickled['condition']
        for par in ['maxGrowthLowerTimeCutoff','maxGrowthUpperTimeCutoff',
                    'allowMaxGrowthrateAtLowerCutoff', 'allowGrowthyieldSlopeNStderrAwayFromZero']:
            self._inheritableParameters[par]=unpickled[par]
        self._setWellIds(unpickled['wellIds'])
        self._setWellIndices(unpickled['wellIndices'],unpickled['activeWellIndices'])
        self._tmp_backgroundIndex=unpickled['backgroundIndex']

    def _serialise(self):
        """
        Generates a dictionary of parameters belonging to this well.

        Note that the main payload, i.e. the optical density, is
        stored in the plate's serialisation.

        For internal use only.
        """
        return dict(
            sampleId=self.sampleid,
            condition=self.condition,
            maxGrowthLowerTimeCutoff=self._inheritableParameters['maxGrowthLowerTimeCutoff'],
            maxGrowthUpperTimeCutoff=self._inheritableParameters['maxGrowthUpperTimeCutoff'],
            allowMaxGrowthrateAtLowerCutoff=self._inheritableParameters['allowMaxGrowthrateAtLowerCutoff'],
            allowGrowthyieldSlopeNStderrAwayFromZero=self._inheritableParameters['allowGrowthyieldSlopeNStderrAwayFromZero'],
            wellIds=self.wellids,
            wellIndices=self._wellIndices,
            activeWellIndices=self._activeWellIndices,
            backgroundIndex=self._backgroundIndex,
            )

    def _invalidate(self):
        """
        Invalidate well. The well is not usable anymore afterwards.

        This is used to make sure code that holds a reference to an
        outdated well fails, instead of producing werid results.

        For internal use only.
        """
        self._memoised={}
        self.sampleid=None
        self.condition=None
        self.wellids=None
        self._inheritableParameters=None
        self._wellIndices=None
        self._activeWellIndices=None
        self._backgroundIndex=None
        self._replicateGroupParent=None
        self.background=None

    def _setMetadata(self,metadata):
        if metadata is None:
            metadata={}
        if 'sample' in metadata:
            if metadata['sample'] is not None:
                self.sampleid=metadata['sample']
            else:
                self.sampleid=''
        if 'condition' in metadata:
            if metadata['condition'] is not None:
                self.condition=metadata['condition']
            else:
                self.condition=''

    def _getMetadata(self):
        return {
            'sample': self.sampleid,
            'condition': self.condition,
            }

    def _setSampleId(self,sampleid,wellIndices=None):
        if sampleid is not None:
            self.sampleid=sampleid
        elif wellIndices is not None:
            wellids=[]
            samplecond={}
            for tcidx in wellIndices:
                tc=self.parentPlate.wells[tcidx]
                if tc.wellids is not None:
                    wellids.extend(tc.wellids)
                # gather all sampe condition pairs
                samplecond[tc.sampleid+" "+tc.condition]=1
                sampleid=tc.sampleid
            if len(samplecond) == 1:
                self.sampleid=sampleid
            else:
                self.sampleid=' '.join(list(samplecond.keys()))
            if len(wellids) > 0:
                self._setWellIds(wellids)
        else:
            raise RuntimeError("setting sampleid failed")

    def _setCondition(self,condition):
        self.condition=condition

    def _setWellIds(self,wellids):
        if wellids is None:
            self.wellids=None
            return
        if type(wellids) is not list:
            raise RuntimeError("wellids should be a list")
        self.wellids=wellids

    def fullId(self,withPlateId=False):
        """
        Return string containing enough information to uniquely identify this Replicate.

        :param withPlateId: If True includes the identifier string of the plate.
        :type withPlateId: bool

        :return: str -- A human readable string identifying this well/set of wells.

        The string contains enough information to uniquely
        identify the well (sample id, condition, well location,
        active child wells)
        """
        fullId=self.sampleid+' '+self.condition
        fullId+=' '+self.activeChildWellIdStr()
        if self.isReplicateGroup():
            fullId+=' (replicate group)'
        if withPlateId and self.parentPlate.plateId is not None:
            fullId+=' '+self.parentPlate.plateId
        return fullId


    def _setWellIndices(self,wellIndices,activeWellIndices=None):
        """
        Set the indices into plate's array of well objects and optionally which ones of these are active.

        For internal use only: should only be used in first initialisation

        :param wellIndices: Indices into the parent-plate's array of well objects.
        :type wellIndices: list(int)
        :param activeWellIndices: Indices into the vector of child wells of this replicate group (from interval [0,len(wellIndices)[)
        :type activeWellIndices: list(int)
        """
        self._wellIndices=wellIndices
        if wellIndices is not None:
            self._setActiveChildWellIndices(activeWellIndices)

    def _setReplicateGroupParent(self,replicateGroupParent):
        """
        For internal use only.

        :param replicateGroupParent: The replicate group that this one should be part of.
        :type replicateGroupParent: Replicate
        """
        self._replicateGroupParent=replicateGroupParent

    def replicateGroupParent(self):
        """
        :return: Replicate -- The replicate group this well is a child of.

        A :py:class:`single-well <.Replicate>` usually is
        part of a replicate group. The child objects of a
        replicate group are used to calculate properties such as
        growth rate, lag time etc. This method returns the
        replicate group parent if this is a single-well object.
        """
        return self._replicateGroupParent

    def isReplicateGroup(self):
        """
        :return: bool -- True if this is a replicate group (not a single well).

        For details see the :py:meth:`documentation of
        replicateGroupParent <.replicateGroupParent>`.
        """
        return self._isReplicateGroup


    def childWellIndices(self):
        """
        Return the indices into the plate's single-well object array.

        :return: list(int) -- Absolute indices of all child wells.

        .. note::

            the indices are *global indices* into the array of
            single-well objects of the plate.
        """
        return self._wellIndices

    def childWells(self):
        """
        Return the single-well child objects.

        :return: list(Replicate) -- child wells of this replicate.
        """
        if self.childWellIndices() is None:
            return []
        tcs=[]
        for idx in self.childWellIndices():
            tcs.append(self.parentPlate.wells[idx])

        return tcs

    def _setActiveChildWellIndices(self,activeWellIndices):
        """
        Set active well indices.

        For internal use only.
        """
        if activeWellIndices is None:
            activeWellIndices = []
            for idx in range(0,len(self.childWellIndices())):
                if self.parentPlate._rawOd[self.childWellIndices()[idx]] is not None:
                    activeWellIndices.append(idx)
        else:
            self._checkActiveWellIndices(activeWellIndices)

        self._activeWellIndices=activeWellIndices

    def _checkActiveWellIndices(self,activeWellIndices):
        """
        Check that active child well indices are in range and no duplicates exists.

        For internal use only.
        """
        includedIndices=set()
        for idx in activeWellIndices:
            if idx < 0 or idx >= len(self.childWellIndices()):
                raise IndexError("local index "+str(idx)+" out of range for "+self.sampleid+" "+self.condition)
            if self.parentPlate._rawOd[self.childWellIndices()[idx]] is None:
                raise RuntimeError("a well without data must not be active (local index "+str(idx)+
                                   ", well index "+str(self.childWellIndices()[idx])+") "+self.fullId())
            # make sure each index is only mentioned once
            if idx in includedIndices:
                raise RuntimeError("local index "+str(idx)+" given multiple times: "+self.fullId())
            includedIndices.add(idx)

    def _calculateRawOd(self):
        """
        Calculate the raw optical density and its variance.

        For internal use only.

        Helper method to calculate the raw optical density (and the
        variance) of this Replicate. This method does not call
        callbacks.
        """
        self._memoised['rawOd']=None
        self._memoised['rawOdVar']=None
        self._checkActiveWellIndices(self.activeChildWellIndices())
        if len(self.activeChildWellIndices()) > 0:
            mat=numpy.zeros([len(self.activeChildWellIndices()),len(self.parentPlate._rawOd[0])])
            ridx=0
            for localdataidx in self.activeChildWellIndices():
                # set row of the matrix
                mat[ridx,]=self.parentPlate._rawOd[self.childWellIndices()[localdataidx]]
                ridx+=1
            # calculate mean of the rawOd for the active data indices
            self._memoised['rawOd'] = mat.mean(axis=0)

        if len(self.activeChildWellIndices()) > 1:
            # calculate variance of the rawOd for the active data indices
            self._memoised['rawOdVar'] = mat.var(axis=0,ddof=1)

    def setActiveChildWellIndices(self,activeWellIndices):
        """
        Set the active wells of this Replicate.

        :param activeWellIndices: Indices into the vector of child wells of this replicate group.
        :type activeWellIndices: list(int)

        A replicate group object consists of multiple wells (its
        children). This method changes which of these are used
        when calculating properties such as growth rate, lag time
        etc.

        .. note::

            these indices are *local* indices into the list of
            well indices of this replicate group, i.e. the maximal
            value should be (number of child wells - 1).
        """
        if not self.isReplicateGroup():
            raise RuntimeError(self.fullId()+": cannot change active child well indices, this is not a replicate group.")
        self._setActiveChildWellIndices(activeWellIndices)
        # tell parent that "self" changed (this may be the background of another replicate)
        self._parametersUpdated('activewells')

    def activateChildWellIndex(self,index,value):
        """
        The child with the given index will be (de-)activated.

        :param index: Index into the vector of child wells of this replicate group.
        :type index: int
        :param value: Set to True to activate, False to deactivate.
        :type index: bool

        For details see the :py:meth:`documentation of
        setActiveChildWellIndices <.Replicate.setActiveChildWellIndices>`.
        """
        activeIndices=self.activeChildWellIndices()
        if value and index not in activeIndices:
            activeIndices.append(index)
            activeIndices=sorted(activeIndices)
        elif not value and index in activeIndices:
            activeIndices.remove(index)

        self.setActiveChildWellIndices(activeIndices)

    def activeChildWellIndices(self):
        """
        Return the indices of the active child-wells.

        :return: list(int) -- Local indices of active child wells.

        .. note::

            these indices are *local* indices into the list of
            well indices of this replicate, i.e. the maximal
            value should be (number of child wells - 1).

        For details see the :py:meth:`documentation of
        setActiveChildWellIndices <.Replicate.setActiveChildWellIndices>`.
        """
        return self._activeWellIndices

    def activeChildWells(self):
        """
        Return the objects of the active child-wells.

        :return: list(Replicate) -- active child wells of this Replicate.

        For details see the :py:meth:`documentation of
        setActiveChildWellIndices <.Replicate.setActiveChildWellIndices>`.
        """
        if self.childWellIndices() is None:
            return []
        tcs=[]
        for localdataidx in self.activeChildWellIndices():
            tcs.append(self.parentPlate.wells[self.childWellIndices()[localdataidx]])

        return tcs

    def activeChildWellIds(self):
        """
        Return ids of active child-wells as list of strings.

        If well ids are set, return these, otherwise return the indices of the wells.

        :return: list(str) -- A list of human readable strings identifying this' active well/set of wells.

        For details see the :py:meth:`documentation of
        setActiveChildWellIndices <.Replicate.setActiveChildWellIndices>`.
        """
        lst=[]
        if self.isReplicateGroup():
            for tc in self.activeChildWells():
                lst.extend(tc.activeChildWellIds())
        else:
            if self.wellids is not None and self.wellids[0] is not None:
                wellid=self.wellids[0]
            else:
                wellid=str(self.childWellIndices()[0])
            lst.append(wellid)
        return lst

    def activeChildWellIdStr(self):
        """
        Return ids of active child-wells as string.

        If well ids are set, return these, otherwise return the indices of the wells.

        :return: str -- A human readable string identifying this' active well/set of wells.
        """
        return ' '.join(self.activeChildWellIds())

    def _setBackgroundIndex(self,indexOfParentReplicateGroup):
        """
        Set the index of the background Replicate object.

        For internal use only.

        :param indexOfParentReplicateGroup: Absolute index of backgound replicate group.
        :type indexOfParentReplicateGroup: int
        """
        self._backgroundIndex=indexOfParentReplicateGroup
        if indexOfParentReplicateGroup is None:
            self.background=None
        else:
            self.background=self.parentPlate.replicateGroups[indexOfParentReplicateGroup]
        self._parametersUpdated(par='backgroundIndex')

    def _clearMemoised(self,par=None):
        """
        Clear memoised ("cached") results.

        For internal use only.

        :param par: The name of the parameter that was changed.
        :type par: str

        See :py:meth:`_parametersUpdated
        <.Plate._parametersUpdated>` for more
        information.
        """
        # if no parameter was given or parameter is not listed, play
        # it safe here and clear all memoised results
        if par is None or par not in Replicate._memoisedDontClear:
            self._memoised.clear()
            return
        # clear all values from memoised results except those that are explicitely excluded
        for key in list(self._memoised.keys()):
            if key not in Replicate._memoisedDontClear[par]:
                self._memoised.pop(key)

    def _parametersUpdated(self,par=None,dontRecurse=False):
        """
        Notify replicate(s) that a parameter changed and memoised results should be deleted.

        For internal use only.

        :param par: The name of the parameter that was changed.
        :type par: str

        The object memoises some results that are expensive to
        calculate. When a parameter is updated, the result may not
        be valid anymore and should get removed from the "cache".
        If par is given, this method can decide which results
        should be removed.
        """
        if self.isReplicateGroup() and not dontRecurse:
            for tc in self.childWells():
                tc._clearMemoised(par)
        self._clearMemoised(par)
        if not self.isReplicateGroup() and self._replicateGroupParent is not None and not dontRecurse:
            self._replicateGroupParent._clearMemoised(par)
        # in case self is the background of some replicates notify these
        # (needs only to be done for 'activewells', as background wells don't
        # depend on any other parameter (only rawOd of background wells is used))
        if par is not None and par == 'activewells' and not dontRecurse:
            self.parentPlate._replicateChanged(self)
        self.parentPlate.modified=True

    def _setExplicitParameter(self,par,val):
        """
        Change value of the given (editible) parameter.

        For internal use only.

        :param par: The name of the parameter that will be changed.
        :type par: str
        :param val: The new value.
        """
        if par in Replicate._isPurePlateParameter and Replicate._isPurePlateParameter[par]:
            raise RuntimeError("_setExplicitParameter: parameter "+par+" cannot be set as this is a plate-wide parameter")
        if par not in self._inheritableParameters:
            raise RuntimeError("_setExplicitParameter: unknown parameter "+par)
        self._inheritableParameters[par]=val
        self._parametersUpdated(par)

    def _getExplicitParameter(self,par):
        """
        Get explicit value of parameter (inherited value is not considered).

        For internal use only.

        :param par: The name of the parameter.
        :type par: str
        """
        if par in Replicate._isPurePlateParameter and Replicate._isPurePlateParameter[par]:
            return None
        if par not in self._inheritableParameters:
            raise RuntimeError("_getExplicitParameter: unknown parameter "+par)
        return self._inheritableParameters[par]

    def getParameter(self,par):
        """
        Return the requested parameter.

        :param par: The name of the parameter.
        :type par: str

        If the parameter is explicitly set for this replicate
        this value returned. Otherwise return an inherited value
        from either the replicate group this one is part of,
        or the plate default.

        See chapter :ref:`parameters <gat parameters>` for details of
        parameter handling and available parameters.
        """
        if par in Replicate._isPurePlateParameter and Replicate._isPurePlateParameter[par]:
            # this is a plate parameter (global parameter), get it from plate object
            return self.parentPlate._getDefaultParameter(par)
        if par not in self._inheritableParameters:
            raise RuntimeError("getParameter: unknown parameter "+par)
        if self._inheritableParameters[par] is not None:
            # parameter is explicitly set
            return self._inheritableParameters[par]
        if self.replicateGroupParent() is not None and self.replicateGroupParent().getParameter(par) is not None:
            # take parameter from parental replicate group
            return self.replicateGroupParent().getParameter(par)
        # neither explicitly set nor implicitely via parental replicate group: use plate default
        return self.parentPlate._getDefaultParameter(par)

    def parameterIsEditible(self,par):
        """
        Return True if this is a parameter that is editible for a Replicate.

        :return: bool -- True if parameter can be edited.

        Some parameters can only be changed per Replicate, some only
        per Plate. This method is used to distinguish between them.

        See chapter :ref:`parameters <gat parameters>` for details of
        parameter handling and available parameters.
        """
        if par in Replicate._isPurePlateParameter and Replicate._isPurePlateParameter[par]:
            return False
        if par not in self._inheritableParameters:
            raise RuntimeError("parameterIsEditible: unknown parameter "+par)

        return True

    def parameterIsExplicitlySet(self,par):
        """
        Return True if this is parameter is explicitly set.

        :param par: The name of the parameter.
        :type par: str

        :return: bool -- True if parameter is explicitly set.

        If a parameter is explicitly set for a replicate it overrides
        an inherited value. This method is used to tell whether this
        is the case.

        See chapter :ref:`parameters <gat parameters>` for details of
        parameter handling and available parameters.
        """
        return self._getExplicitParameter(par) is not None

    def activeChildReplicatesHaveExplicitParameter(self,par):
        """
        Return True if this is a ReplicateGroup and for at least one of the active child-wells the given parameter is explicitly set.

        :param par: The name of the parameter.
        :type par: str

        :return: bool -- True if parameter is explicitly set in one of the active child wells.

        See chapter :ref:`parameters <gat parameters>` for details of
        parameter handling and available parameters.
        """
        if not self.isReplicateGroup():
            return False
        activeChilds=self.activeChildWells()
        if len(activeChilds) == 0:
            return False
        for childtc in activeChilds:
            if childtc._getExplicitParameter(par) is not None:
                return True
        return False

    def hdCorrectionLinear(self):
        """
        :return: float -- Coefficient of linear term of high density correction.
        """
        return self.getParameter('hdCorrectionLinear')

    def hdCorrectionQuadratic(self):
        """
        :return: float -- Coefficient of quadratic term of high density correction.
        """
        return self.getParameter('hdCorrectionQuadratic')

    def hdCorrectionCubic(self):
        """
        :return: float -- Coefficient of cubic term of high density correction.
        """
        return self.getParameter('hdCorrectionCubic')

    def smoothingK(self):
        """
        :return: int -- Degree of the smoothing spline.
        """
        return self.getParameter('smoothingK')

    def smoothingS(self):
        """
        :return: float -- Smoothing factor used to choose the number of knots.
        """
        return self.getParameter('smoothingS')

    def maxGrowthLowerTimeCutoff(self):
        """
        :return: float -- Lower limit of interval in which the maximal growth should be searched.
        """
        return self.getParameter('maxGrowthLowerTimeCutoff')

    def maxGrowthUpperTimeCutoff(self):
        """
        :return: float -- Upper limit of interval in which the maximal growth should be searched.
        """
        return self.getParameter('maxGrowthUpperTimeCutoff')

    def logOdCutoff(self):
        """
        :return: float -- Cutoff value of log(OD).
        """
        return self.getParameter('logOdCutoff')

    def slidingWindowSize(self):
        """
        :return: int -- Number of datapoints of sliding windows.

        The value that is used for local exponential fit (growth rate) and linear regression (growth yield).
        """
        return self.getParameter('slidingWindowSize')

    def lagAtLogOdEquals(self):
        """
        :return: float -- Value of log(OD) used to define the lag time.
        """
        return self.getParameter('lagAtLogOdEquals')

    def allowMaxGrowthrateAtLowerCutoff(self):
        """
        :return: bool -- Whether or not the maximal growth rate is allowed at the interval limits.
        """
        return self.getParameter('allowMaxGrowthrateAtLowerCutoff')

    def allowGrowthyieldSlopeNStderrAwayFromZero(self):
        """
        :return: int -- Number of standard errors the slope may be away from zero to detect stationary state of yield.
        """
        return self.getParameter('allowGrowthyieldSlopeNStderrAwayFromZero')

    def setMaxGrowthLowerTimeCutoff(self,t):
        """Set lower limit of interval in which the maximal growth should be searched."""
        self._setExplicitParameter('maxGrowthLowerTimeCutoff',t)

    def setMaxGrowthUpperTimeCutoff(self,t):
        """Set upper limit of interval in which the maximal growth should be searched."""
        self._setExplicitParameter('maxGrowthUpperTimeCutoff',t)

    def setAllowMaxGrowthrateAtLowerCutoff(self,boolval):
        """(Dis)allow the maximal growth rate at the lower interval limit."""
        self._setExplicitParameter('allowMaxGrowthrateAtLowerCutoff',boolval)

    def setAllowGrowthyieldSlopeNStderrAwayFromZero(self,numTimesStderr):
        """Set number of standard errors the slope may be away from zero to detect stationary state of yield."""
        self._setExplicitParameter('allowGrowthyieldSlopeNStderrAwayFromZero',numTimesStderr)

    def rawOd(self):
        """
        Return the raw readout of the optical density.

        The raw optical density is not preprocessed, i.e. the
        background has not been subtracted and it is not corrected for
        non-linearities at high densities.
        """
        if 'rawOd' not in self._memoised:
            self._calculateRawOd()
        return self._memoised['rawOd']

    def rawOdVar(self):
        """
        Return the variance of raw readout of the optical density.
        """
        if 'rawOdVar' not in self._memoised:
            self._calculateRawOd()
        return self._memoised['rawOdVar']

    def _calculateOd(self):
        """
        Set od and odVar.

        For internal use only.
        """
        self._memoised['od']=None
        self._memoised['odVar']=None
        if self.rawOd() is None:
            return
        if self.background is None:
            return
        if self.hdCorrectionLinear() is None or self.hdCorrectionQuadratic() is None or self.hdCorrectionCubic() is None:
            return

        # subtract background from rawOd
        rawdiff=self.rawOd()-self.background.rawOd()
        self._memoised['od']=self.hdCorrectionLinear()*rawdiff + self.hdCorrectionQuadratic()*rawdiff**2 + self.hdCorrectionCubic()*rawdiff**3
        if self.rawOdVar() is not None and self.background.rawOdVar() is not None:
            # simple error propagation of _uncorrelated_ variables
            rawdiffVar=self.rawOdVar() + self.background.rawOdVar()
            self._memoised['odVar']=((self.hdCorrectionLinear()
                                      + 2*self.hdCorrectionQuadratic()*rawdiff
                                      + 3*self.hdCorrectionCubic()*rawdiff**2)**2 * rawdiffVar)

    def od(self):
        """
        Return the background- and high-density corrected optical density.
        """
        if 'od' not in self._memoised:
            self._calculateOd()
        return self._memoised['od']

    def odVar(self):
        """
        Return the variance of the background- and high-density corrected optical density.
        """
        if 'odVar' not in self._memoised:
            self._calculateOd()
        return self._memoised['odVar']

    def derivative(self):
        """
        Return the (left) derivative of the background- and high-density corrected optical density.
        """
        if 'derivative' in self._memoised:
            return self._memoised['derivative']

        self._memoised['derivative']=None
        od=self.od()
        if od is not None:
            # simple "left" derivative
            self._memoised['derivative']=numpy.diff(od)/numpy.diff(self.time)

        return self._memoised['derivative']

    def smoothedOd(self):
        """
        Return smoothing spline of optical density.

        :return: numpy.array(float) -- Smoothed optical density.
        """
        # FIXME instead of smoothing the time series itself it may be worthwhile to try
        # http://dx.doi.org/10.5402/2011/164564

        if 'smoothedOd' in self._memoised:
            return self._memoised['smoothedOd']
        if self.od() is None:
            return None

        self._memoised['smoothedOd']=None
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")

            f=UnivariateSpline(self.time,self.od(),k=self.smoothingK(),s=self.smoothingS())
            if len(w):
                print("smoothing for sample '"+self.sampleid+"' condition '"+self.condition+"' failed")
            else:
                self._memoised['smoothedOd']=f(self.time)

        return self._memoised['smoothedOd']

    def smoothedOdDerivative(self):
        """
        Return derivative of smoothed optical density.

        :return: numpy.array(float) -- Derivative of smoothed optical density.
        """
        if 'smoothedOdDerivative' in self._memoised:
            return self._memoised['smoothedOdDerivative']

        self._memoised['smoothedOdDerivative']=None
        smoothedOd=self.smoothedOd()
        if smoothedOd is not None:
            self._memoised['smoothedOdDerivative']=numpy.diff(smoothedOd)/numpy.diff(self.time)

        return self._memoised['smoothedOdDerivative']

    def logOd(self):
        """
        Return the logarithm of the background- and high-density corrected optical density.
        """
        if 'logOd' in self._memoised:
            return self._memoised['logOd']

        self._memoised['logOd']=None
        od=self.od()
        if od is not None:
            self._memoised['logOd']=numpy.zeros([len(self._memoised['od'])])
            idcs=self._memoised['od'] >= 1e-35 # FIXME an abitrary threshold
            self._memoised['logOd'][idcs]=numpy.log(self._memoised['od'][idcs])
            self._memoised['logOd'][~idcs]=numpy.nan

        return self._memoised['logOd']

    def logOdSmoothed(self):
        """
        Return smoothed logarithm optical density.

        :return: numpy.array(float) -- Smoothed log(OD).
        """
        if 'logOdSmoothed' in self._memoised:
            return self._memoised['logOdSmoothed']

        self._memoised['logOdSmoothed']=None
        try:
            with warnings.catch_warnings(record=True) as w:
                # Cause all warnings to always be triggered.
                warnings.simplefilter("always")

                nonnanidcs=~numpy.isnan(self.logOd()) # only use non-nan values for the interpolation
                #print 'len(nonnanidcs)', nonnanidcs.shape
                f=UnivariateSpline(self.time[nonnanidcs],self.logOd()[nonnanidcs],k=self.smoothingK(),s=self.smoothingS())

                if len(w):
                    print("logsmoothing for sample '"+self.fullId()+"' failed (warning)")
                else:
                    self._memoised['logOdSmoothed']=f(self.time)
        except :
            print("logsmoothing for sample '"+self.fullId()+"' failed (error)")

        return self._memoised['logOdSmoothed']

    def logOdDerivative(self):
        """
        Return derivative of logarithmised optical density.

        :return: numpy.array(float) -- Derivative of logarithmised optical density.
        """
        if self.od() is None:
            return None
        return numpy.diff(self.logOd())/numpy.diff(self.time)

    def logOdDerivativeFromNonLog(self):
        """
        Return derivative of logarithmised optical density (by chain rule).

        :return: numpy.array(float) -- Derivative of logarithmised optical density.
        """
        if self.od() is None:
            return None

        od0=self.od()[:-1]
        nonzeroidcs = od0 != 0
        logodderivative=numpy.zeros([len(self.od())-1])
        logodderivative[nonzeroidcs]=1/od0[nonzeroidcs]*self.derivative()[nonzeroidcs]
        logodderivative[~nonzeroidcs]=numpy.nan

        return logodderivative

    def logOdDerivativeFromNonLogSmoothed(self):
        """
        Return smoothed derivative of logarithmised optical density (by chain rule).

        :return: numpy.array(float) -- Derivative of logarithmised optical density.
        """
        if self.smoothedOd() is None:
            return None

        od0=self.smoothedOd()[:-1]
        nonzeroidcs = od0 != 0
        logodderivative=numpy.zeros([len(self.od())-1])
        logodderivative[nonzeroidcs]=1/od0[nonzeroidcs]*self.smoothedOdDerivative()[nonzeroidcs]
        logodderivative[~nonzeroidcs]=numpy.nan

        return logodderivative

    def expFitsOd0Mu(self):
        """
        Return parameters for fitted exponential functions (two parameters: OD0 and mu).

        :return: numpy.array(float), numpy.array(float), numpy.array(float), numpy.array(float) -- mean(mu), var(mu), mean(od0), var(od0)

        Exponential functions are fitted piecewise to :math:`w` (see
        :ref:`fit window`) data points. Note that the inital value
        :math:`\odz_{t0}` assumes that :math:`t=0` for the first data
        point of the fit window; to convert it to the real time you
        should use the formula :math:`\odz = \odz_{t0} e^{-\mu
        t_0}`.

        """
        if 'expFitsOd0Mu' in self._memoised:
            c=self._memoised['expFitsOd0Mu']
            return c['mu'], c['muvar'], c['od0'], c['od0var']

        c={}
        c['mu'], c['muvar'], c['od0'], c['od0var'] = self._localODexpFit(fitOd0=True,useSmoothed=False)
        self._memoised['expFitsOd0Mu']=c
        return c['mu'], c['muvar'], c['od0'], c['od0var']

    def expFitsMu(self):
        """
        Return parameters for fitted exponential functions (one parameter: mu).

        :return: numpy.array(float), numpy.array(float) -- mean(mu), var(mu)

        Exponential functions are fitted piecewise to :math:`w` (see
        :ref:`fit window`) data points. Only the growth rate
        :math:`\mu` is fitted, the first data point from :math:`w` is
        used as initial value of the exponential function.
        """
        if 'expFitsMu' in self._memoised:
            c=self._memoised['expFitsMu']
            return c['mu'], c['muvar']

        c={}
        c['mu'], c['muvar'], od0Dummy, od0varDummy = self._localODexpFit(fitOd0=False,useSmoothed=False)
        self._memoised['expFitsMu']=c
        return c['mu'], c['muvar']

    def _localODexpFit(self,fitOd0=True,useSmoothed=False):
        """
        Return parameters for fitted exponential functions.

        :return: numpy.array(float), numpy.array(float), numpy.array(float), numpy.array(float) -- mean(mu), var(mu), mean(od0), var(od0)

        For internal use only.
        """
        if self.od() is None:
            return None, None, None, None
        slidingWindowSize=self.slidingWindowSize()
        if slidingWindowSize is None:
            raise RuntimeError("slidingWindowSize for "+self.fullId(withPlateId=True)+" is None")
        if slidingWindowSize < 3:
            raise RuntimeError("slidingWindowSize for "+self.fullId(withPlateId=True)+" is too small: "+str(slidingWindowSize))

        if self.isReplicateGroup():
            # here we average over the underlying wells
            mu=numpy.zeros([len(self.activeChildWellIndices()), len(self.od())-slidingWindowSize])
            od0=numpy.zeros([len(self.activeChildWellIndices()), len(self.od())-slidingWindowSize])
            i=0
            for tc in self.activeChildWells():
                mu[i], muvarDummy, od0[i], od0varDummy = tc._localODexpFit(fitOd0=fitOd0,useSmoothed=useSmoothed)
                i+=1

            mumean, muvar = maskedArrayToMeanVar(mu, ddof=1, axis=0)
            od0mean, od0var = maskedArrayToMeanVar(od0, ddof=1, axis=0)

            return mumean, muvar, od0mean, od0var

        if useSmoothed is True:
            thisod=self.smoothedOd()
        else:
            thisod=self.od()
        if thisod is None:
            return None, None, None, None

        mu=numpy.zeros([len(self.od())-slidingWindowSize])
        od0=numpy.empty([len(self.od())-slidingWindowSize])
        od0.fill(numpy.nan)
        for i in range(0,mu.shape[0]):
            if fitOd0 is True:
                # function to fit: OD_fit(t[i+j]) = OD_i * exp(mu*(t[i+j] - t[i]))
                # i.e. OD_i and mu are fit parameters, and a good guess for OD_i is OD_meas(t[i])
                f = lambda tminti, *p: p[0] * numpy.exp(p[1]*(tminti))
                with warnings.catch_warnings(record=True) as w:
                    # Turn all warnings into errors.
                    warnings.simplefilter("error")
                    try:
                        popt, pcov = scipy.optimize.curve_fit(f,
                                                              xdata=self.time[i:i+slidingWindowSize]-self.time[i],
                                                              ydata=thisod[i:i+slidingWindowSize],
                                                              p0=[thisod[i],1])
                        mu[i]=popt[1]
                        od0[i]=popt[0]
                    except (RuntimeError, OptimizeWarning) as e:
                        mu[i]=numpy.nan
                        od0[i]=numpy.nan
            else:
                # function to fit: OD_fit(t[i+j]) = OD[t[i]] * exp(mu*(t[i+j] - t[i]))
                # i.e. mu is the fit parameter
                f = lambda tminti, *p: thisod[i] * numpy.exp(p[0]*(tminti))

                with warnings.catch_warnings(record=True) as w:
                    # Turn all warnings into errors.
                    warnings.simplefilter("error")
                    try:
                        popt, pcov = scipy.optimize.curve_fit(f,
                                                              xdata=self.time[i:i+slidingWindowSize]-self.time[i],
                                                              ydata=thisod[i:i+slidingWindowSize],
                                                              p0=[1])
                        mu[i]=popt[0]
                    except (RuntimeError, OptimizeWarning) as e:
                        mu[i]=numpy.nan

        return mu, None, od0, None

    def maxGrowthrateFromLogOdDerivative(self):
        """
        Return parameters for exponential function at maximal growth rate (determined from log(OD) derivative).

        :return: mean(mu), var(mu), mean(od0), var(od0), mean(maxt), var(maxt), mean(lag), var(lag), method, statuses
        """
        return self._maxGrowthrate(method='nonlogsmoothed')

    def maxGrowthrate(self):
        """
        Return parameters for exponential function at maximal growth rate (determined from fitted exponential functions).

        :return: mean(mu), var(mu), mean(od0), var(od0), mean(maxt), var(maxt), mean(lag), var(lag), method, statuses
        """
        return self._maxGrowthrate(method='expfit')

    def _maxGrowthrate(self,method='expfit',detailsInMessage=True):
        """
        Return parameters for exponential function at maximal growth rate.

        :param method: The method used to determine the maximal growth.
        :type method: str

        :return: mean(mu), var(mu), mean(od0), var(od0), mean(maxt), var(maxt), mean(lag), var(lag), method, statuses

        For internal use only.

        FIXME enhance documentation of returned values!
        """
        lagAtLogOdEquals=self.lagAtLogOdEquals()
        allowMuMaxAtLowerCutoff=self.allowMaxGrowthrateAtLowerCutoff()

        method2statustext={
            'nonlogsmoothed': 'smoothed',
            'expfit': 'exp. fit',
            }
        methodtxt=method2statustext[method]

        if self.isReplicateGroup():
            # here we average over the underlying wells
            mumax=numpy.zeros([len(self.activeChildWellIndices())])
            od0max=numpy.zeros([len(self.activeChildWellIndices())])
            maxt=numpy.zeros([len(self.activeChildWellIndices())])
            lag=numpy.zeros([len(self.activeChildWellIndices())])
            allstatuses=StatusMessage()
            statuses=StatusMessage()
            alllagstatuses=StatusMessage()
            lagstatuses=StatusMessage()
            i=0
            for tc in self.activeChildWells():
                mumax[i], d1, od0max[i], d2, maxt[i], d3, lag[i], d4, themethod, status = tc._maxGrowthrate(method,detailsInMessage=False)
                if status is not None:
                    allstatuses.addStatus(status)
                    alllagstatuses.addStatus(status.statusesWithKey('lag ('+methodtxt+'):'))
                if status is not None and ~numpy.isnan(mumax[i]):
                    statuses.addStatus(status.statusesWithKey('max. growth rate ('+methodtxt+'):'))
                if status is not None and ~numpy.isnan(lag[i]):
                    lagstatuses.addStatus(status.statusesWithKey('lag ('+methodtxt+'):'))
                i+=1
            idcs=numpy.isnan(mumax)
            lagidcs=numpy.isnan(lag)
            if numpy.all(idcs):
                return None, None, None, None, None, None, None, None, method, allstatuses
            mumaxmean, mumaxvar =   maskedArrayToMeanVar(mumax, idcs, ddof=1)
            od0maxmean, od0maxvar = maskedArrayToMeanVar(od0max, idcs, ddof=1)
            maxtmean, maxtvar =     maskedArrayToMeanVar(maxt, idcs, ddof=1)
            lagmean, lagvar =       maskedArrayToMeanVar(lag, lagidcs, ddof=1)

            if lagmean is None:
                statuses.addStatus(alllagstatuses)
            else:
                statuses.addStatus(lagstatuses)
            return mumaxmean, mumaxvar, od0maxmean, od0maxvar, maxtmean, maxtvar, lagmean, lagvar, method, statuses

        if method == 'nonlogsmoothed':
            mu = self.logOdDerivativeFromNonLogSmoothed()
            t = self.time[:-1]
            logod = None
            if self.logOd() is not None:
                logod = self.logOd()[:-1]
            od0 = None
            if self.smoothedOd() is not None:
                od0 = self.smoothedOd()[:-1]
        elif method == 'expfit':
            mu, muvar, od0, od0var = self.expFitsOd0Mu()
            t = self.time[int(math.floor(self.slidingWindowSize()/2.)):int(-math.ceil(self.slidingWindowSize()/2.))]
            logod = None
            if self.logOd() is not None:
                logod = self.logOd()[int(math.floor(self.slidingWindowSize()/2.)):int(-math.ceil(self.slidingWindowSize()/2.))]

        # only consider non-None growth rates
        fidcs = numpy.isfinite(mu) if mu is not None else numpy.array([False for i in range(0,len(t))])
        idcs=numpy.copy(fidcs)
        # only consider growth rates where log(OD) is above cutoff
        if self.logOdCutoff() is not None and self.logOd() is not None:
            cidcs=notNanAndGreaterEqual(logod,self.logOdCutoff())
            idcs=numpy.logical_and(idcs,cidcs)

        if not numpy.any(idcs):
            return None, None, None, None, None, None, None, None, method, StatusMessage(
                key='max. growth rate ('+methodtxt+'):',shortmsg='growthrate({methodtxt}):noMu',
                longmsg='no growth rate could be determined',
                severity=Severity.failed, methodtxt=methodtxt)

        # only consider growth rates which lie within time-cutoff
        maxGrowthLowerTimeCutoff=self.maxGrowthLowerTimeCutoff()
        maxGrowthUpperTimeCutoff=self.maxGrowthUpperTimeCutoff()
        if maxGrowthLowerTimeCutoff is not None:
            lidcs=t >= maxGrowthLowerTimeCutoff
            idcs=numpy.logical_and(idcs,lidcs)

        if maxGrowthUpperTimeCutoff is not None:
            uidcs=t <= maxGrowthUpperTimeCutoff
            idcs=numpy.logical_and(idcs,uidcs)

        if not numpy.any(idcs):
            return None, None, None, None, None, None, None, None, method, StatusMessage(
                key='max. growth rate ('+methodtxt+'):',shortmsg='growthrate({methodtxt}):noMuWithinLimits',
                longmsg='no growth rate within cutoff limits',
                severity=Severity.failed, methodtxt=methodtxt)

        maxidx=numpy.argmax(mu[idcs])
        mumax=mu[idcs][maxidx]
        maxt=t[idcs][maxidx]
        globalmaxidx=idcs.nonzero()[0][maxidx]
        # convert od0 at time[i] to od0 at t=0
        # mu*(t-t0) + math.log(od0_t0) == mu*t - mu*t0 + log(od0_t0)
        # ==> log(od0)=log(od0_t0)-mu*t0
        # ==> od0 = od0_t0 * exp(-mu*t0)
        if method == 'expfit':
            od0max=od0[idcs][maxidx] * math.exp(-mumax*self.time[:-self.slidingWindowSize()][idcs][maxidx])
        else:
            od0max=od0[idcs][maxidx] * math.exp(-mumax*t[idcs][maxidx])

        # unrealistic location of maximal growth is neglected
        if mumax <= 0:
            return None, None, None, None, None, None, None, None, method, StatusMessage(
                key='max. growth rate ('+methodtxt+'):',shortmsg='growthrate({methodtxt}):mumaxLt0',
                longmsg='maximal growth rate less than zero',
                severity=Severity.failed, methodtxt=methodtxt)
        if od0max <= 0:
            return None, None, None, None, None, None, None, None, method, StatusMessage(
                key='max. growth rate ('+methodtxt+'):',shortmsg='growthrate({methodtxt}):od0maxLt0',
                longmsg='initial OD is less than zero',
                severity=Severity.failed, methodtxt=methodtxt)
        if lagAtLogOdEquals is not None and mumax*maxt+math.log(od0max) < lagAtLogOdEquals:
            details=' ('+str(mumax*maxt+math.log(od0max))+' < '+str(lagAtLogOdEquals)+')' if detailsInMessage else ''
            return None, None, None, None, None, None, None, None, method, StatusMessage(
                key='max. growth rate ('+methodtxt+'):',shortmsg='growthrate({methodtxt}):lagGtValAtMumax',
                longmsg='OD at lag is greater than OD at time of maximal growth'+details,
                severity=Severity.failed, methodtxt=methodtxt)
        status=None

        # there may be three reasons why this max is located at a "lower" or "upper" cutoff:
        # - logOdCutoff forces datapoints next to the maximum to NaN
        # - maxGrowth{Lower,Upper}TimeCutoff is next to the maximum
        # - the growthrate could not be detected next to the maximum
        key=None
        if globalmaxidx == 0 or not fidcs[globalmaxidx-1]:
            details=' (t='+str(maxt)+')' if detailsInMessage else ''
            key='max. growth rate ('+methodtxt+'):'
            shortmsg='growthrate({methodtxt}):maxAtFirstNonNan'
            longmsg='next to maximum there is no growth rate defined'+details
        elif maxGrowthLowerTimeCutoff is not None and not lidcs[globalmaxidx-1]:
            wouldgivemax=', would need a lower cutoff of '+str(t[lidcs.nonzero()[0][0]-1]) if 0 < lidcs.nonzero()[0][0]-1 else ''
            details=' (t='+str(maxt)+wouldgivemax+')' if detailsInMessage else ''
            key='max. growth rate ('+methodtxt+'):'
            shortmsg='growthrate({methodtxt}):maxAtLowerCutoff'
            longmsg='located at lower cutoff'+details
        elif self.logOdCutoff() is not None and not cidcs[globalmaxidx-1]:
            details=' (at t<'+str(maxt)+')' if detailsInMessage else ''
            key='max. growth rate ('+methodtxt+'):'
            shortmsg='growthrate({methodtxt}):maxAtLogOdCutoff'
            longmsg='located at log(OD) cutoff'+details
        if key is not None:
            if not allowMuMaxAtLowerCutoff:
                return None, None, None, None, None, None, None, None, method, StatusMessage(
                    key=key,
                    shortmsg=shortmsg+'Removed',
                    longmsg='maximal growth rate rejected: '+longmsg,
                    severity=Severity.failed, methodtxt=methodtxt)
            else:
                status=StatusMessage(
                    key=key,shortmsg=shortmsg,
                    longmsg='maximal growth rate: '+longmsg,
                    severity=Severity.warning, methodtxt=methodtxt)

        if globalmaxidx == len(t)-1 or not fidcs[globalmaxidx+1]:
            details=' (t='+str(maxt)+')' if detailsInMessage else ''
            return None, None, None, None, None, None, None, None, method, StatusMessage(
                key='max. growth rate ('+methodtxt+'):',shortmsg='growthrate({methodtxt}):maxAtUpperCutoffRemoved',
                longmsg='maximal growth rate rejected: there is no growth rate defined for greater times'+details,
                severity=Severity.failed, methodtxt=methodtxt)
        elif maxGrowthUpperTimeCutoff is not None and not uidcs[globalmaxidx+1]:
            wouldgivemax=', would need an upper cutoff of '+str(t[uidcs.nonzero()[0][-1]+1]) if len(t) > uidcs.nonzero()[0][-1]+1 else ''
            details=' (t='+str(maxt)+wouldgivemax+')' if detailsInMessage else ''
            return None, None, None, None, None, None, None, None, method, StatusMessage(
                key='max. growth rate ('+methodtxt+'):',shortmsg='growthrate({methodtxt}):maxAtUpperCutoffRemoved',
                longmsg='maximal growth rate rejected: located at upper cutoff'+details,
                severity=Severity.failed, methodtxt=methodtxt)
        elif self.logOdCutoff() is not None and not cidcs[globalmaxidx+1]:
            # this will probably never happen, but just in case
            details=' (at t>'+str(maxt)+')' if detailsInMessage else ''
            return None, None, None, None, None, None, None, None, method, StatusMessage(
                key='max. growth rate ('+methodtxt+'):',shortmsg='growthrate({methodtxt}):maxAtLogOdCutoffRemoved',
                longmsg='maximal growth rate rejected: located at log(OD) cutoff'+details,
                severity=Severity.failed, methodtxt=methodtxt)

        lag=None
        lagstatus=None
        if lagAtLogOdEquals is None:
            pass # no message, it is clear that without this the lag cannot be calculated
        else:
            # mumax <= 0 and od0max <= 0 have been exluded above
            lag=(lagAtLogOdEquals-math.log(od0max))/(mumax)
            if lag < 0:
                details=' (at t='+str(lag)+')' if detailsInMessage else ''
                lagstatus = StatusMessage(
                    key='lag ('+methodtxt+'):',shortmsg='lag({methodtxt}):lessThanZero',
                    longmsg='lag rejected: negative'+details,
                    severity=Severity.failed, methodtxt=methodtxt)
                lag = None
        if status is not None and lagstatus is not None:
            listStatus=StatusMessage()
            listStatus.addStatus(status)
            listStatus.addStatus(lagstatus)
            status=listStatus
        elif lagstatus is not None:
            status=lagstatus

        return mumax, None, od0max, None, maxt, None, lag, None, method, status

    def odSlopemaxIntercept(self):
        """
        Return maximal slope, intercept, time of maximal slope and time index of the (linear) OD.

        :return: float, float, float, float, float, float, numpy.array(int), StatusMessage -- mean(slope), var(slope), mean(intercept), var(intercept), mean(timemax), var(timemax), timemaxIdcs, status
        """

        if self.isReplicateGroup():
            # here we average over the underlying wells
            slopemax=numpy.empty([len(self.activeChildWellIndices())])
            slopemax.fill(numpy.nan)
            interceptmax=numpy.empty([len(self.activeChildWellIndices())])
            interceptmax.fill(numpy.nan)
            timemax=numpy.empty([len(self.activeChildWellIndices())])
            timemax.fill(numpy.nan)
            timemaxIdx=numpy.empty([len(self.activeChildWellIndices())])
            timemaxIdx.fill(numpy.nan)
            allstatuses=StatusMessage()
            statuses=StatusMessage()
            i=0
            for tc in self.activeChildWells():
                slopemax[i], slopemaxVarDummy, interceptmax[i], interceptmaxVarDummy, timemax[i], timemaxVarDummy, timemaxIndices, status = tc.odSlopemaxIntercept()
                if timemaxIndices is not None:
                    timemaxIdx[i]=timemaxIndices[0]
                if status is not None:
                    allstatuses.addStatus(status)
                if status is not None and ~numpy.isnan(slopemax[i]):
                    statuses.addStatus(status)
                i+=1
            idcs=numpy.isnan(slopemax)
            if numpy.all(idcs):
                return None, None, None, None, None, None, None, allstatuses
            slopemaxmean, slopemaxvar =   maskedArrayToMeanVar(slopemax, idcs, ddof=1)
            interceptmaxmean, interceptmaxvar = maskedArrayToMeanVar(interceptmax, idcs, ddof=1)
            timemaxmean, timemaxvar =     maskedArrayToMeanVar(timemax, idcs, ddof=1)

            return slopemaxmean, slopemaxvar, interceptmaxmean, interceptmaxvar, timemaxmean, timemaxvar, timemaxIdx, statuses

        if self.smoothedOdDerivative() is None:
            return None, None, None, None, None, None, None, StatusMessage(
                key='odSlopemaxIntercept',
                shortmsg='odSlopemaxIntercept:noSlope',
                longmsg='derivative of smoothed optical density could not be calculated (probably smoothing failed)',
                severity=Severity.failed)

        smoothedod=self.smoothedOd()[:-1]
        smoothedderivative=self.smoothedOdDerivative()
        t=self.time[:-1]

        maxGrowthLowerTimeCutoff=self.maxGrowthLowerTimeCutoff()
        if maxGrowthLowerTimeCutoff is not None:
            # lower cutoff index: first element that is greater equal lower cutoff value
            greaterLowerCutoff = (self.time[:-1] >= maxGrowthLowerTimeCutoff).nonzero()[0]
            # if there is no such index, set the index to a value triggers below status 'noSlopeWithinLimits'
            tlower=greaterLowerCutoff[0] if len(greaterLowerCutoff) else len(self.time)
        else:
            tlower=0
        # upper cutoff index
        tupper=len(self.time)-1
        if tupper-tlower < 1:
            return None, None, None, None, None, None, None, StatusMessage(
                key='odSlopemaxIntercept',
                shortmsg='odSlopemaxIntercept:noSlopeWithinLimits',
                longmsg='no derivative of smoothed optical density for times greater than maxGrowthLowerTimeCutoff',
                severity=Severity.failed)

        if self.logOdCutoff() is not None and self.logOd() is not None:
            idcs=notNanAndGreaterEqual(self.logOd()[tlower:tupper],self.logOdCutoff())
        else:
            idcs=numpy.array([True for i in range(tlower,tupper)],dtype=bool)
        if not idcs.any():
            return None, None, None, None, None, None, None, StatusMessage(
                key='odSlopemaxIntercept',
                shortmsg='odSlopemaxIntercept:noSlopeForLogodGreaterCutoff',
                longmsg='no derivative of smoothed optical density for which log(OD) is greater equal cutoff',
                severity=Severity.failed)

        localMaxIdx = numpy.argmax(smoothedderivative[tlower:tupper][idcs])
        timemaxIdx = idcs.nonzero()[0][localMaxIdx] + tlower
        slopemax = smoothedderivative[timemaxIdx]
        timemax=self.time[timemaxIdx]
        interceptmax=smoothedod[timemaxIdx] - slopemax * timemax

        # unrealistic location of maximal growth is neglected
        if slopemax < 0:
            return None, None, None, None, None, None, None, StatusMessage(
                key='odSlopemaxIntercept',
                shortmsg='odSlopemaxIntercept:slopemaxLt0',
                longmsg='no positive slope could be determined',
                severity=Severity.failed)
        if smoothedod[timemaxIdx] <= 0:
            return None, None, None, None, None, None, None, StatusMessage(
                key='odSlopemaxIntercept',
                shortmsg='odSlopemaxIntercept:odAtMaxLt0',
                longmsg='optical density at maximal slope is less than zero',
                severity=Severity.failed)
        if interceptmax >= 0:
            return None, None, None, None, None, None, None, StatusMessage(
                key='odSlopemaxIntercept',
                shortmsg='odSlopemaxIntercept:interceptGt0',
                longmsg='at intercept (t=0) optical density is greater than zero',
                severity=Severity.failed)

        return slopemax, None, interceptmax, None, timemax, None, numpy.array([timemaxIdx]), None

    def _windowMeanVar(self,fromidx=None,toidx=None,useSmoothed=False):
        """
        Mean and its variance of a sliding window.

        For each index within [fromidx:toidx] the mean and
        variance for a window of slidingWindowSize datapoints is
        rerturned.

        For internal use only.
        """
        slidingWindowSize=self.slidingWindowSize()
        if fromidx is None:
            fromidx=0
        if toidx is None:
            toidx=self.od().size-slidingWindowSize
        if toidx > self.od().size-slidingWindowSize:
            toidx=self.od().size-slidingWindowSize
        if fromidx > toidx:
            return None, None
        if useSmoothed is True:
            thisod=self.smoothedOd()
        else:
            thisod=self.od()

        mean=numpy.zeros([toidx-fromidx])
        sigma2=numpy.zeros([toidx-fromidx])
        for i in range(fromidx,toidx):
            window=thisod[i:i+slidingWindowSize]
            mean[i-fromidx]=window.mean()
            sigma2[i-fromidx]=window.var()

        return mean, sigma2

    def _windowSlope(self,fromidx=None,toidx=None,useSmoothed=False):
        """
        Slope estimated by linear regression in a sliding window.

        For each index within [fromidx:toidx] a linear regression
        is performed for a window of slidingWindowSize datapoints.

        For internal use only.
        """
        slidingWindowSize=self.slidingWindowSize()
        if fromidx is None:
            fromidx=0
        if toidx is None:
            toidx=self.od().size-slidingWindowSize
        if toidx > self.od().size-slidingWindowSize:
            toidx=self.od().size-slidingWindowSize
        if fromidx >= toidx:
            return None, None, None
        if useSmoothed is True:
            thisod=self.smoothedOd()
        else:
            thisod=self.od()

        slope = numpy.zeros([toidx-fromidx])
        slopeStdErr = numpy.zeros([toidx-fromidx])
        intercept = numpy.zeros([toidx-fromidx])
        for i in range(fromidx,toidx):
            window=thisod[i:i+slidingWindowSize]
            slope[i-fromidx], intercept[i-fromidx], thisr_value, thisp_value, slopeStdErr[i-fromidx] = scipy.stats.linregress(
                self.time[i:i+slidingWindowSize],
                window)

        return slope, slopeStdErr, intercept

    def growthyield(self):
        """
        Return estimate of the growth yield by looking for a slope (from linear regression in a time window) compatible with zero.

        :return: float, float, float, float, StatusMessage -- mean(yield), var(yield), mean(timeOfYield), var(timeOfYield), status
        """
        return self._growthyield(useSmoothed=True)

    def _growthyield(self,useSmoothed=False):
        """
        Return estimate of the growth yield by looking for a slope (from linear regression in a time window) compatible with zero.

        :param useSmoothed: Whether to use the smoothed optical density.
        :type useSmoothed: bool

        :return: float, float, float, float, StatusMessage -- mean(yield), var(yield), mean(timeOfYield), var(timeOfYield), status

        For internal use only.
        """

        if self.od() is None:
            return None, None, None, None, StatusMessage(key='growthyield',shortmsg='growthyield:noOd',
                                                                             longmsg='no non-raw optical density',
                                                                             severity=Severity.failed)

        if self.isReplicateGroup():
            # here we average over the underlying wells
            growthyield=numpy.zeros([len(self.activeChildWellIndices())])
            tgrowthyield=numpy.zeros([len(self.activeChildWellIndices())])
            allstatuses=StatusMessage()
            statuses=StatusMessage()
            i=0
            for tc in self.activeChildWells():
                growthyield[i], d1, tgrowthyield[i], d2, status = tc._growthyield(useSmoothed)
                allstatuses.addStatus(status)
                if status is not None and ~numpy.isnan(growthyield[i]):
                    statuses.addStatus(status)
                i+=1
            idcs=numpy.isnan(growthyield)
            if numpy.all(idcs):
                return None, None, None, None, allstatuses
            (growthyieldmean, growthyieldvar)=maskedArrayToMeanVar(growthyield, idcs, ddof=1)
            (tgrowthyieldmean, tgrowthyieldvar)=maskedArrayToMeanVar(tgrowthyield, idcs, ddof=1)

            return growthyieldmean, growthyieldvar, tgrowthyieldmean, tgrowthyieldvar, statuses


        # the timepoint of the yield should be above the maximal growth (linear scale, because this is more stable)
        slidingWindowSize=self.slidingWindowSize()
        slopemax, slopemaxVar, interceptmax, interceptmaxVar, timemaxslope, timemaxslopeVar, timemaxIndices, plainSlopeStatus=self.odSlopemaxIntercept()
        if timemaxslope is None:
            return None, None, None, None, StatusMessage(key='growthyield',shortmsg='growthyield:noSlopemax',
                                                                             longmsg='no maximal slope',
                                                                             severity=Severity.failed)
        timemaxIdx=timemaxIndices[0]
        if timemaxIdx >= len(self.time)-slidingWindowSize:
            return None, None, None, None, StatusMessage(key='growthyield',shortmsg='growthyield:slopemaxUnusable',
                                                                             longmsg='unusable maximal slope',
                                                                             severity=Severity.failed)
        odAtMaxLinSlope = slopemax*timemaxslope+interceptmax

        slope, slopeStdErr, intercept = self._windowSlope(fromidx=timemaxIdx,useSmoothed=useSmoothed)
        mean, sigma2 = self._windowMeanVar(fromidx=timemaxIdx,useSmoothed=useSmoothed)
        tmb=int(math.floor(slidingWindowSize/2.))
        tmt=int(math.ceil(slidingWindowSize/2.))

        # a slice within [skipStart:] :
        # 1. the the OD at the growthyield should be larger than than the OD at maximal linear growth
        # 2. the time of the growthyield should be larger than the time of the maximal growth
        # 3. the fitted slope (linear fit within a window around time) at that time should be zero
        validIndices=(numpy.logical_and(mean > odAtMaxLinSlope,
                                        numpy.logical_and(slope+slopeStdErr >= 0, slope-slopeStdErr < 0)))

        status=None
        # a second try, this time with larger errors (looping with n=[2:allowNStderrDeviations]):
        allowNStderrDeviations=self.allowGrowthyieldSlopeNStderrAwayFromZero()
        n=2
        while (not validIndices.any() and n <= allowNStderrDeviations):
            validIndices=(numpy.logical_and(mean > odAtMaxLinSlope,
                                            numpy.logical_and(slope + n*slopeStdErr >= 0, slope - n*slopeStdErr < 0)))
            status=StatusMessage(key='growthyield',shortmsg='growthyield:within{severity2}Stderr',
                                                     longmsg='slope zero within {severity2} standard errors',
                                                     severity=Severity.warning,severity2=n)
            n+=1

        if not validIndices.any():
            return None, None, None, None, StatusMessage(key='growthyield',shortmsg='growthyield:noValidIndices',
                                                                             longmsg='no window with a slope compatible with zero',
                                                                             severity=Severity.failed)

        # pick the largest value of those that are valid
        yieldidx=numpy.argmax(mean[validIndices])
        growthyield=mean[validIndices][yieldidx]
        tgrowthyield=self.time[tmb+timemaxIdx:-tmt][validIndices][yieldidx]

        if growthyield<0:
            return None, None, None, None, StatusMessage(key='growthyield',shortmsg='growthyield:negativeYield',
                                                                             longmsg='invalid yield (negative)',
                                                                             severity=Severity.failed)

        return growthyield, None, tgrowthyield, None, status

    @staticmethod
    def growthrateToDoublingTime(mu,mu_var):
        """
        Calculate doubling time from growth rate.

        :param mu: The growth rate.
        :type mu: float
        :param mu_var: The variance of the growth rate.
        :type mu_var: float

        :return: float, float -- doubling time, variance of doubling time
        """
        doublingtime=None
        doublingtimevar=None
        if mu is not None and ~numpy.isnan(mu):
            doublingtime=math.log(2)/mu
            if mu_var is not None and ~numpy.isnan(mu_var):
                # error propagation from doub(mu)=ln(2)/mu :
                # doubvar(mu)= (d{doub(mu)}/d{mu})^2 * muvar = (-ln(2)/mu^2)^2 * muvar = ln(2)^2/mu^4 * muvar
                doublingtimevar=math.pow(math.log(2),2)/math.pow(mu,4)*mu_var
        return doublingtime, doublingtimevar
