"""
This module implements the :py:class:`ClsReplicate` class for CATHODE.

Chronological life span Analysis Tool for High-throughput Optical
Density Experiments (CATHODE) Replicate class.
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

import math
import numpy

from platereader.replicate import Replicate
from platereader.numpytools import maskedArrayToMeanVar
from platereader.statusmessage import StatusMessage, Severity

class ClsReplicate(object):
    """
    """

    def __init__(self,parent=None,
                 listOfOdReplicates=None,days=None,
                 wellIndices=None,isReplicateGroup=False,
                 wellids=None,sampleid=None,condition=None,
                 status=StatusMessage()):
        self.sampleid=sampleid
        self.condition=condition
        self.wellids=wellids
        self._isReplicateGroup=isReplicateGroup
        self._wellIndices=wellIndices

        self.days=days
        self._setActiveChildWellIndices(None) # NOTE this activates all children
        self.clsParent=parent
        self._odReplicates=listOfOdReplicates
        self.initDiffStatus=status

    def _serialiseLightweight(self):
        return dict(
            sampleId=self.sampleid,
            condition=self.condition,
            wellIds=self.wellids,
            wellIndices=self._wellIndices,
            activeWellIndices=self._activeWellIndices,
            )

    def _deserialiseLightweight(self,unpickled):
        """
        Apply unpickled data on top of existing plates.

        When re-reading files we would like to preserve the state of
        activated wells. This method does exactly this: an initialised
        Cls object (basically an array of re-read
        Plate objects) is amended with data that was
        previously saved.

        For internal use only.
        """
        if self.sampleid != unpickled['sampleId']:
            raise RuntimeError('sampleids do not match: '+self.sampleid+' '+unpickled['sampleId'])
        if self.condition != unpickled['condition']:
            raise RuntimeError('conditions do not match: '+self.condition+' '+unpickled['condition'])
        if self.wellids != unpickled['wellIds']:
            raise RuntimeError('wellids do not match: '+self.wellids+' '+unpickled['wellIds'])
        if self._wellIndices != unpickled['wellIndices']:
            raise RuntimeError('wellIndices do not match: '+self._wellIndices+' '+unpickled['wellIndices'])
        self._activeWellIndices=unpickled['activeWellIndices']

    def isReplicateGroup(self):
        """
        :return: bool -- True if this is a replicate group.

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
        Return the single-well cls objects.

        :return: list(Replicate) -- child cls wells of this replicate.
        """
        if self.childWellIndices() is None:
            return []
        tcs=[]
        for idx in self.childWellIndices():
            tcs.append(self.clsParent.clsWells[idx])

        return tcs

    def _setActiveChildWellIndices(self,activeWellIndices):
        """
        a helper function that does not call callbacks
        """
        if activeWellIndices is None:
            activeWellIndices=list(range(0,len(self.childWellIndices())))

        if len(activeWellIndices) > 0:
            includedIndices=set()
            for localdataidx in activeWellIndices:
                if localdataidx < 0 or localdataidx >= len(self.childWellIndices()):
                    raise RuntimeError("local index "+str(localdataidx)+" out of range")
                if localdataidx in includedIndices:
                    raise RuntimeError("local index "+str(localdataidx)+" given multiple times")
                includedIndices.add(localdataidx)

        self._activeWellIndices=activeWellIndices

    def setActiveChildWellIndices(self,activeWellIndices):
        """
        Set the active wells indices.

        A replicate group object consists of multiple wells (its
        children). This method changes which of these are used
        when calculating properties such as viability etc.

        .. note::

            these indices are *local* indices into the list of
            well indices of this replicate group, i.e. the maximal
            value should be (number of child wells - 1).
        """
        if not self.isReplicateGroup():
            raise RuntimeError(self.fullId()+": cannot change active child well indices, this is not a replicate group.")
        self._setActiveChildWellIndices(activeWellIndices)
        self.clsParent.modified=True

    def activateChildWellIndex(self,index,value):
        """
        The child with the given index will be (de-)activated.

        :param index: Index into the vector of child ClsReplicate of this replicate group.
        :type index: int
        :param value: Set to True to activate, False to deactivate.
        :type index: bool

        For details see the :py:meth:`documentation of
        setActiveChildWellIndices <.ClsReplicate.setActiveChildWellIndices>`.
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
        Return the indices of the active children.

        :return: list(int) -- Local indices of active child wells.

        .. note::

            these indices are *local* indices into the list of
            well indices of this replicate, i.e. the maximal
            value should be (number of child wells - 1).

        For details see the :py:meth:`documentation of
        setActiveChildWellIndices <.ClsReplicate.setActiveChildWellIndices>`.
        """
        return self._activeWellIndices

    def activeChildWells(self):
        """
        :return: list(Replicate) -- active child wells of this Replicate.
        """
        if self.childWellIndices() is None:
            return []
        tcs=[]
        for localdataidx in self.activeChildWellIndices():
            tcs.append(self.clsParent.clsWells[self.childWellIndices()[localdataidx]])

        return tcs

    def activeChildWellIds(self):
        """
        Return ids of active wells as list of strings.

        If well ids are set, return these, otherwise return the indices of the wells.

        :return: list(str) -- A list of human readable strings identifying this' active well/set of wells.
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
        Return ids of active wells as string.

        If well ids are set, return these, otherwise return the indices of the wells.

        :return: str -- A human readable string identifying this' active well/set of wells.
        """
        return ' '.join(self.activeChildWellIds())

    def odReplicatesNdays(self):
        """
        :return: list(Replicate), list(float) -- replicate objects of the underlying plates, days.
        """
        return self._odReplicates, self.days

    def fullId(self):
        """
        Return string containing enough information to uniquely identify this ClsReplicate.

        :return: str -- A human readable string identifying this well/set of wells.

        The string contains enough information to uniquely
        identify the well (sample id, condition, well location,
        active child wells)
        """
        return self.sampleid+' '+self.condition+' '+self.activeChildWellIdStr()

    def viability(self):
        """
        Calculate the viability.

        :return: list(float), list(float), list(float), StatusMessage -- days, viability, var(viability), status.
        """
        if self.isReplicateGroup():
            # here we average over the underlying replicates
            viability=numpy.zeros([len(self.activeChildWellIndices()), self.days.shape[0]])

            allstatuses=StatusMessage()
            statuses=StatusMessage()
            i=0
            for clstc in self.activeChildWells():
                days, viability[i], viabilityvar, status = clstc.viability()
                if status is not None:
                    allstatuses.addStatus(status)
                if status is not None and days is not None:
                    statuses.addStatus(status)
                i+=1

            idcs=numpy.isnan(viability)
            if numpy.all(idcs):
                allstatuses.addStatus(self.initDiffStatus)
                return None, None, None, allstatuses

            viabilitymean, viabilityvar = maskedArrayToMeanVar(viability, idcs=idcs, ddof=1, axis=0)

            statuses.addStatus(self.initDiffStatus)
            return self.days, viabilitymean, viabilityvar, statuses

        viability=numpy.zeros([len(self._odReplicates)])
        viabilityvar=numpy.zeros([len(self._odReplicates)])
        mu_ef, mu_ef_var, od0_ef, od0_ef_var, maxt_ef, maxt_ef_var, lag_ef__0, lag_ef_var__0, method_ef, status = self._odReplicates[0].maxGrowthrate()
        if lag_ef__0 is None:
            return None, None, None, StatusMessage(
                key='viability:noInitialLag', severity=Severity.failed,
                longmsg='lag could not be extract for first timepoint')

        tcidx=-1
        for tc in self._odReplicates:
            tcidx+=1
            if tc is None:
                viability[tcidx]=None
                viabilityvar[tcidx]=None
                print('WARNING no Replicate defined for day '+str(self.days[tcidx])+' and sample/condition '+self._odReplicates[0].fullId())
                continue

            mu_ef, mu_ef_var, od0_ef, od0_ef_var, maxt_ef, maxt_ef_var, lag_ef, lag_ef_var, method_ef, status = tc.maxGrowthrate()
            doublingtime, doublingtimevar = Replicate.growthrateToDoublingTime(mu_ef,mu_ef_var)
            if lag_ef is not None and doublingtime is not None:
                deltaT=lag_ef-lag_ef__0
                # f(deltaT,doubling) = 2^(- deltaT/doubling )  [ == e^(-ln(2) deltaT/doubling) ]
                viability[tcidx]=math.pow(2,-deltaT/doublingtime)
            elif doublingtime is None:
                # if we could not extract a doubling time we assume that there was no growth at all
                viability[tcidx]=0.
                viabilityvar[tcidx]=0.
            else:
                viability[tcidx]=None
            if lag_ef_var__0 is not None and lag_ef_var is not None and doublingtimevar is not None:
                deltaTvar=lag_ef_var+lag_ef_var__0
                # df/ddeltaT  = - ln(2)/doubling  2^(- deltaT/doubling )
                # df/doubling = + ln(2) deltaT 1/doubling^2   2^(- deltaT/doubling )
                viabilityvar[tcidx]=( math.pow( math.log(2)/doublingtime * math.pow(2,-deltaT/doublingtime),2 ) * doublingtimevar
                                      + math.pow( math.log(2)*deltaT  * math.pow(doublingtime,-2) * math.pow(2,-deltaT/doublingtime),2 ) * deltaTvar)
            else:
                viabilityvar[tcidx]=None

        return self.days, viability, viabilityvar, self.initDiffStatus

    def survivalIntegral(self):
        """
        Calculate the survival integral.

        :return: float, float, StatusMessage -- survival integral, var(survival integral), status.
        """
        if self.isReplicateGroup():
            # here we average over the underlying replicates
            si=numpy.zeros([len(self.activeChildWellIndices())])

            allstatuses=StatusMessage()
            statuses=StatusMessage()
            i=0
            for clstc in self.activeChildWells():
                si[i], sivar, status = clstc.survivalIntegral()
                if status is not None:
                    allstatuses.addStatus(status)
                if status is not None and days is not None:
                    statuses.addStatus(status)
                i+=1

            idcs=numpy.isnan(si)
            if numpy.all(idcs):
                allstatuses.addStatus(self.initDiffStatus)
                return None, None, allstatuses

            simean, sivar = maskedArrayToMeanVar(si, ddof=1)

            statuses.addStatus(self.initDiffStatus)
            return simean, sivar, statuses

        days, viability, viabilityvar, initDiffStatus=self.viability()
        si=None
        if viability is not None and days is not None:
            si=numpy.trapz(viability,x=days)
        return si,None,None
