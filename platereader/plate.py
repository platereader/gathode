"""
This module implements the :py:class:`Plate` class for GATHODE.

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
import math
import numpy
import json
import bz2

import platereader
from platereader.replicate import Replicate
from platereader.statusmessage import StatusMessage, Severity
from platereader.csvunicode import CsvFileUnicodeWriter, CsvFileUnicodeReader
from platereader.parser import tecan, bioscreen

class Plate(object):
    """
    Class containing the wells and holding plate-wide parameters.
    """

    _parser2module = platereader.parser.modulenameToModule(
        list(platereader.parser.getModulesOfNamespace(platereader.parser)),
        replace='platereader.parser.',
        lower=True)

    _isNotPlateParameter={
        'allowMaxGrowthrateAtLowerCutoff': True,
        'allowGrowthyieldSlopeNStderrAwayFromZero': True,
        }

    def __init__(self,filename=None,fileformat=None,
                 time=None,rawOds=None,
                 sampleIds=None,conditions=None,wellids=None,plateId=None):
        """
        Constructor.

        If filename is not None and fileformat is None some heuristics
        are used to identify the file format.

        :param filename: name of serialised Plate or ascii file exported by the plate reader.
        :type filename: str

        :param fileformat: string indicating the format ('gat', 'tecan')
        :type fileformat: str

        :param time: array of timepoints when optical density was measured
        :type time: numpy.array(float)

        :param rawOds: list of optical density arrays
        :type rawOds: list( numpy.array(float) )

        :param sampleIds: list of sample names corresponding to the array of optical densities
        :type sampleIds: list(str)

        :param conditions: list of conditions under which the samples where grown
        :type conditions: list(str)

        :param plateId: name of this plate
        :type plateId: str
        """
        self.plateId=None
        self._rawOd=None
        self.wells=None
        self.time=None
        self.temperature=None
        self.timeunit=None
        self._inheritableParameters={}
        # default parameters
        self._inheritableParameters['maxGrowthLowerTimeCutoff']=None
        self._inheritableParameters['maxGrowthUpperTimeCutoff']=None
        self._inheritableParameters['allowMaxGrowthrateAtLowerCutoff']=False
        self._inheritableParameters['allowGrowthyieldSlopeNStderrAwayFromZero']=1
        # pure plate parameters
        self._inheritableParameters['logOdCutoff']=None
        self._inheritableParameters['lagAtLogOdEquals']=-5
        self._inheritableParameters['slidingWindowSize']=10
        self._inheritableParameters['hdCorrectionLinear']=None
        self._inheritableParameters['hdCorrectionQuadratic']=None
        self._inheritableParameters['hdCorrectionCubic']=None
        self._inheritableParameters['smoothingK']=5
        self._inheritableParameters['smoothingS']=0.01
        self._loadStatus=StatusMessage()
        self._capitaliseBackgroundIds=['blank','background']
        self._clearMetadata()
        if filename is not None:
            if not os.path.exists(filename):
                raise IOError("No such file or directory: '"+filename+"'")
            if fileformat is None:
                if filename.endswith('.gat'):
                    fileformat='gat'
                else:
                    scorefileformat=[]
                    for fileformat in Plate._parser2module:
                        score=Plate._parser2module[fileformat].isPlateFormat(filename)
                        if score > 0.:
                            scorefileformat.append({'score': score, 'fileformat': fileformat})
                    scorefileformat = sorted(scorefileformat, key=lambda k: k['score'],reverse=True)
                    if not len(scorefileformat):
                        raise Plate.UnknownFileFormat(filename,detailedError='Cannot determine file format')
                    fileformat=scorefileformat[0]['fileformat']
            if fileformat == 'gat':
                self._load(filename)
            elif fileformat in Plate._parser2module:
                time, rawOd, sampleIds, conditions, plateId, temperature, wellids=Plate._parser2module[fileformat].parse(filename)
                self._initFromArrays(time,rawOd,sampleIds,conditions,plateId=plateId,temperature=temperature,wellids=wellids)
            else:
                raise Plate.UnknownFileFormat(filename,serFormat=fileformat)
            self.readfileformat=fileformat
        elif rawOds is not None:
            self._initFromArrays(time,rawOds,sampleIds,conditions,plateId=plateId,wellids=wellids)
        else:
            raise RuntimeError('could not construct Plate, neither filename nor arrays given')
        self.modified=False

    def _clearReplicateGroups(self):
        if hasattr(self,'replicateGroups'):
            for tc in self.replicateGroups:
                # NOTE invalidating here so code holding references to these fails
                tc._invalidate()

        self.replicateGroups=None
        self._backgroundGroupIndices=None
        self._sampleConditionToReplicateGroupIdcs=None # an associative array mapping replicate groups by sample ID to a list of Replicate object indices
        self._conditionToReplicateGroupIdx=None  # an associative array mapping condition to a list of replicate group object indices

    def _clearMetadata(self):
        self._clearReplicateGroups()
        self._setBackgroundForAllReplicates(None)
        self._conditionToWellIdx=None # an associative array mapping condition to a list of Replicate objects
        self._sampleConditionToWellIdcs=None # an associative array mapping wells (sample IDs) to a list of Replicate object indices

    def _load(self,filename):
        with bz2.BZ2File(filename, 'r') as rfile:
            pickled=rfile.read().decode("utf-8")

        try:
            unpickled = json.loads(pickled)
        except ValueError as err:
            raise Plate.UnknownFileFormat(filename,detailedError=str(err))
        return self._deserialise(unpickled,filename)

    def _deserialise(self,unpickled,filename):
        if 'format' not in unpickled:
            raise Plate.UnknownFileFormat(filename,detailedError='no "format" keyword found in file')

        serFormatVersion=unpickled['formatversion'] if 'formatversion' in unpickled else 'undefined'
        if unpickled['format'] != 'opticaldensityplate' or serFormatVersion != '1':
            raise Plate.UnknownFileFormat(filename,serFormat=unpickled['format'],serFormatVersion=serFormatVersion)

        parkeys=[
            # default parameters
            'maxGrowthLowerTimeCutoff',
            'maxGrowthUpperTimeCutoff',
            'allowMaxGrowthrateAtLowerCutoff',
            'allowGrowthyieldSlopeNStderrAwayFromZero',
            # pure plate parameters
            'logOdCutoff',
            'lagAtLogOdEquals',
            'slidingWindowSize',
            'hdCorrectionLinear',
            'hdCorrectionQuadratic',
            'hdCorrectionCubic',
            'smoothingK',
            'smoothingS'
            ]

        # reset these to make sure defaults given to constructor are not used for serialised plate
        for par in self._inheritableParameters:
            self._inheritableParameters[par]=None

        self.plateId=unpickled['plateId']
        self.time=numpy.array(unpickled['time'],dtype=float)
        self.timeunit=unpickled['timeunit']
        # defaut parameters, some of which can be overridden by the individual replicates
        for par in parkeys:
            self._inheritableParameters[par]=unpickled[par] 

        if 'temperature' in unpickled:
            self.temperature=numpy.array(unpickled['temperature'],dtype=float)
        self._rawOd=[]
        for lst in unpickled['rawOd']:
            self._rawOd.append(numpy.array(lst,dtype=float))
        self.wells=[]
        for tcup in unpickled['wells']:
            self.wells.append(Replicate(_unpickled=tcup,parentPlate=self,_serialiseFormat=unpickled['format']))
        self.replicateGroups=[]
        for tcup in unpickled['replicateGroup']:
            comptc=Replicate(_unpickled=tcup,parentPlate=self,_serialiseFormat=unpickled['format'],isReplicateGroup=True)
            self.replicateGroups.append(comptc)
            # set parental replicate group of the children
            for childtc in comptc.childWells():
                childtc._setReplicateGroupParent(comptc)

        # deferred to here: set the background index
        for tc in self.wells:
            if tc._tmp_backgroundIndex is not None:
                tc._setBackgroundIndex(tc._tmp_backgroundIndex)
        for tc in self.replicateGroups:
            if tc._tmp_backgroundIndex is not None:
                tc._setBackgroundIndex(tc._tmp_backgroundIndex)
        # reset background indices, as these have been initialised
        # before setting the replicate's backgrounds
        self._backgroundWellIndices=None
        self._backgroundGroupIndices=None

        self._setBackgroundStatus()

    def _serialise(self):
        """
        Generates a dictionary of the plate data and parameters.

        For internal use only.
        """
        parkeys=[
            # default parameters
            'maxGrowthLowerTimeCutoff',
            'maxGrowthUpperTimeCutoff',
            'allowMaxGrowthrateAtLowerCutoff',
            'allowGrowthyieldSlopeNStderrAwayFromZero',
            # pure plate parameters
            'logOdCutoff',
            'lagAtLogOdEquals',
            'slidingWindowSize',
            'hdCorrectionLinear',
            'hdCorrectionQuadratic',
            'hdCorrectionCubic',
            'smoothingK',
            'smoothingS'
            ]
        sr=dict()
        sr["format"]='opticaldensityplate'
        sr["formatversion"]='1' # this is an unsigned integer
        sr['plateId']=self.plateId
        sr['time']=self.time.tolist()
        sr['timeunit']=self.timeunit
        for key in parkeys:
            sr[key]=self._inheritableParameters[key]
        if self.temperature is not None:
            sr['temperature']=self.temperature.tolist()
        sr['rawOd']=[]
        for raw in self._rawOd:
            sr['rawOd'].append(raw.tolist())
        sr['wells']=[]
        for tc in self.wells:
            sr['wells'].append(tc._serialise())
        sr['replicateGroup']=[]
        for tc in self.replicateGroups:
            sr['replicateGroup'].append(tc._serialise())

        return sr

    def save(self,filename):
        """
        Saves the plate content in a file.

        :param filename: Name of the file.
        :type filename: str

        :return: StatusMessage/None -- non-fatal notifications.
        """
        status=None
        if not filename.endswith('.gat'):
            root, ext = os.path.splitext(filename)
            status=StatusMessage(
                key='Saving file',shortmsg='wrongExtension',
                longmsg=('GATHODE uses a file extension that is different from"'+ext+'". '
                         +'This means that a future version of this program will not be able to open this file with the graphical user interface. '
                         +'Please make save the file with the ".gat" extension.'),
                severity=Severity.warning)
        sr=self._serialise()
        pickled = json.dumps(sr)
        with bz2.BZ2File(filename, 'w') as wfile:
            wfile.write(pickled.encode('utf-8'))
        self.modified=False
        return status

    def _explicitlySetParsInChildWells(self,par):
        """
        Explicitly set parameters in wells to their inherited values.

        This can be used when replicate groups get removed
        (e.g. setting new metadata) but the parameters should be
        preserved. You most likely want to call
        _reduceExplicitParameter once new replicate groups have been
        created.

        For internal use only.
        """
        for tc in self.replicateGroups:
            for child in tc.childWells():
                # copy parameters from replicate group to the child
                child._setExplicitParameter(par,child.getParameter(par))

    def _reduceExplicitParameter(self,par):
        """
        Sets the parameter par of self and wells such that it is shared by most of its children.

        For internal use only.

        :param par: the parameter for which a smaller set of values is created
        :type par: string
        """
        # check what could be the plate default for this parameter
        parvals=Plate._getChildParvalOccurrence(self,par)
        platedefault=Plate._chooseDefaultFromOccurrences(parvals)
        # set parameters in replicate groups; if one of a groups's children has the same value
        # as the platedefault use that one, otherwise try find another value for the group
        for tc in self.replicateGroups:
            Plate._reduceExplicitParametersHelper(tc,par,platedefault)
        # now set the plate default (important: this has to be done *after* the Replicates are changed!)
        Plate._reduceExplicitParametersHelper(self,par,platedefault)

        return platedefault

    @staticmethod
    def _reduceExplicitParametersHelper(obj,par,parentdefault):
        """
        Helper function for _reduceExplicitParameter

        For internal use only.

        :param obj: the parameter for that a smaller set of values is created
        :type obj: Plate/Replicates
        :param par: the parameter for that a smaller set of values is created
        :type par: string

        Will be called with plate and Replicate objects.
        """
        # gather occurrence of each value for this parameter in children
        parvals=Plate._getChildParvalOccurrence(obj,par)
        # the value that occurs most often will become the replicate group's value
        newdefaultval=Plate._chooseDefaultFromOccurrences(parvals,parentdefault)
        # only if none of the children got value None we can copy values up
        if newdefaultval is None:
            return
        # delete consensus value from children
        for child in Plate._getChildren(obj):
            if newdefaultval == child.getParameter(par):
                child._setExplicitParameter(par,None)
        # set consensus value for replicate group parent
        obj._setExplicitParameter(par,newdefaultval)

    @staticmethod
    def _getChildParvalOccurrence(obj,par):
        """
        Return count of parameter values of all leaf children (at the lowest level of the hierarchy).

        For internal use only.

        :return: dict -- { value1: countValue1, value2: countValue2, ...}
        """
        if isinstance(obj, Replicate) and not obj.isReplicateGroup():
            # this is a single well
            val=obj.getParameter(par)
            return {val: 1}
        else:
            parvals={}
            for child in Plate._getChildren(obj):
                childparvals=Plate._getChildParvalOccurrence(child,par)
                # assemble childrens' results into the main dictionary
                for val in childparvals:
                    if val not in parvals:
                        parvals[val]=0
                    parvals[val]+=childparvals[val]
            return parvals

    @staticmethod
    def _chooseDefaultFromOccurrences(parvals,parentdefault=None):
        """
        Return the value of a parameter that occurs most often in leaf children.

        For internal use only.

        Can be called both without parentdefault (for the whole plate)
        and with parentdefault (for ReplicateGroups).

        :param parvals: output from _getChildParvalOccurrence
        :type parvals: dict

        :return: float -- the most occuring parameter value for this plate or ReplicateGroup
        """
        if None in parvals:
            return None
        maxcnt=0
        maxval=None
        parvalkeys=list(parvals.keys())
        parvalkeys.sort()
        for val in parvalkeys:
            # if there is a value corresponding to the plate default choose that one
            if parentdefault is not None and val == parentdefault:
                return parentdefault
            # choose maximal occurring as default
            if parvals[val] > maxcnt:
                maxval=val
                maxcnt=parvals[val]

        return maxval

    @staticmethod
    def _getChildren(obj):
        if isinstance(obj, Plate):
            # this is a plate
            return obj.replicateGroups
        else:
            # this is a replicate group
            return obj.childWells()

    @staticmethod
    def capitaliseId(sampleId,capitaliseThese):
        """
        Capitalise id if in given list.

        :param sampleId: sample id; if this matches capitaliseThese it will be capitalised
        :type sampleId: str

        :param capitaliseThese: list of sample ids that correspond to samples that should be capitalised
        :type capitaliseThese: list(str)

        :return: str -- sample id (capitalised if it matches one of capitaliseThese)
        """
        for bgid in capitaliseThese:
            if sampleId.upper() == bgid.upper():
                return bgid.upper()
        return sampleId

    def _initFromArrays(self,time,rawOd,sampleIds,conditions,plateId=None,temperature=None,wellids=None):
        """
        Initialises a plate from numpy arrays.

        For internal use only.

        :param time: array of timepoints when optical density was measured
        :type time: numpy.array(float)

        :param rawOd: list of optical density arrays
        :type rawOd: list( numpy.array(float) )

        :param sampleIds: list of sample names corresponding to the array of optical densities
        :type sampleIds: list(str)

        :param conditions: list of conditions under which the samples where grown
        :type conditions: list(str)

        :param plateId: name of this plate
        :type plateId: str

        :param temperature: array of the temperature
        :type time: numpy.array(float)

        :param wellids: array of ids for the wells (e.g. A1 to P24)
        :type wellids: list(str)
        """

        if len(rawOd) != len(sampleIds):
            raise RuntimeError('number of raw optical density arrays is different from number of sample ids')
        if len(sampleIds) != len(conditions):
            raise RuntimeError('number of sample ids is different from number of conditions')
        if wellids is not None and len(wellids) != len(set(wellids)):
            raise RuntimeError('ids in wellids are not unique')

        self.plateId=plateId

        self.time=time/3600.
        self.timeunit="h"

        self._rawOd=rawOd

        # make sure that background is correctly identified even if case is different
        newSampleIds=[]
        for sampleid in sampleIds:
            newSampleIds.append(Plate.capitaliseId(sampleid,self._capitaliseBackgroundIds))

        # create replicate objects for single wells from data (NOTE ids may exist multiple times, therefore this is not an associative array)
        self.wells=[]
        tcidx=0
        for sampleid in newSampleIds:
            wellid = [wellids[tcidx]] if wellids is not None else None
            self.wells.append(Replicate(self,[tcidx],sampleid,conditions[tcidx],wellid))
            # NOTE that on purpose this index is only increased for samples (not for time, temperature, ...)
            tcidx+=1

        self._createReplicateGroupsFromSampleIdsNConditions()
        # use guessed background sampleIds to set background of single well and replicate groups
        self._setBackgroundForAllReplicates(self._guessBackgroundSampleIds())

    def wellMetadataOk(self,metadata):
        """
        Check that the given metadata (i.e. sample id, growth condition) is valid and can be applied.

        This basically checks that there is the right amount of
        metadata entries and these contain sample ids and conditions.

        :param metadata: array of metadata dictionaries
        :type metadata: list(dict)

        :return: bool, StatusMessage -- True if ok, False otherwise (and a StatusMessage with details)
        """

        if len(metadata) != len(self.wells):
            return False, StatusMessage(
                key='Wrong metadata length:',shortmsg='metadata:wrongLength',
                longmsg=('Number of metadata entries ('+str(len(metadata))+
                         ') is different from number of wells '+str(len(self.wells))),
                severity=Severity.failed)

        idx=0
        for metdat in metadata:
            idx+=1
            if len(metdat.keys()) != 2 or 'sample' not in metdat or 'condition' not in metdat:
                thekeys='"'+('" "'.join(sorted(metdat.keys())))+'"' if len(metdat.keys()) else 'nothing'
                return False, StatusMessage(
                    key='Wrong metadata elements:',shortmsg='metadata:wrongLength',
                    longmsg=('metadata for entry '+str(idx)+' contains '+thekeys+
                             ', but should contain "condition" and "sample"'),
                    severity=Severity.failed)

        return True, StatusMessage()

    def setWellMetadata(self,metadata):
        """
        Set the metadata (e.g. sample id, growth condition) of the wells.

        :param metadata: array of metadata dictionaries
        :type metadata: list(dict)
        """
        metok, message = self.wellMetadataOk(metadata)
        if not metok:
            raise Plate.BadMetadata(str(message))

        # propagate parameters to the wells before deleting replicate groups
        for par in self.wells[0]._inheritableParameters.keys():
            self._explicitlySetParsInChildWells(par)

        # clear everything that depends on metadata
        self._clearMetadata()

        # set metadata of the wells
        wellit=self.wells.__iter__()
        for metdat in metadata:
            metdat['sample']=Plate.capitaliseId(metdat['sample'],self._capitaliseBackgroundIds)
            well=next(wellit)
            well._setMetadata(metdat)

        # create replicate groups based on sample ids and conditions
        self._createReplicateGroupsFromSampleIdsNConditions()
        # use guessed background sampleIds to set background of single well and replicate groups
        self._setBackgroundForAllReplicates(self._guessBackgroundSampleIds())

        # propagate parameters from the wells to the replicate groups (or plate) if possible
        for par in self.wells[0]._inheritableParameters.keys():
            self._reduceExplicitParameter(par)

    def wellMetadata(self):
        """
        Return the metadata of the wells.

        :return: list(dict) -- metadata
        """
        metadata=[]
        for well in self.wells:
            metadata.append(well._getMetadata())
        return metadata

    def _setupBackgroundIndices(self):
        """
        Set self._backgroundGroupIndices and self._backgroundWellIndices.

        Records the indices of tc.background (which are rpelicate
        groups) for all wells and replicateGroups and also the indices
        of the underlying background wells.

        For internal use only.
        """
        self._backgroundGroupIndices=set()
        self._backgroundWellIndices=set()

        if self.wells:
            for tc in self.wells:
                if tc.background:
                    self._backgroundGroupIndices.add(self._indexOfReplicateGroup(tc.background))
        if self.replicateGroups:
            for tc in self.replicateGroups:
                if tc.background:
                    self._backgroundGroupIndices.add(self._indexOfReplicateGroup(tc.background))

        for idx in self._backgroundGroupIndices:
            for chldidx in self.replicateGroups[idx].childWellIndices():
                self._backgroundWellIndices.add(chldidx)

    def _guessBackgroundSampleIds(self):
        """
        Guess sample ids of background wells ("BLANK" or "BACKGROUND")

        For internal use only.
        """
        backgroundKeys={}
        for tc in self.wells:
            if tc.sampleid == "BLANK" or tc.sampleid == "BACKGROUND":
                backgroundKeys[tc.sampleid]=1

        backgroundSampleIds=sorted(list(backgroundKeys.keys()))
        return backgroundSampleIds

    def _setBackgroundStatus(self):
        """
        Add conditions/samples for which no background was found to self._loadStatus

        This should be called when the background was set for some wells/replicate groups.

        For internal use only.
        """
        self._loadStatus.removeStatusesWithKey('No background samples:')
        self._loadStatus.removeStatusesWithKey('No background for some samples:')

        backgroundSampleIds=set()
        for idx in self.backgroundReplicateGroupIndices():
            backgroundSampleIds.add(self.replicateGroups[idx].sampleid)
        for idx in self.backgroundWellIndices():
            backgroundSampleIds.add(self.wells[idx].sampleid)

        if len(backgroundSampleIds) < 1:
            self._loadStatus.addStatus(
                StatusMessage(
                    key='No background samples:',shortmsg='plateinit:noBackground',
                    longmsg=('No background (blank) wells could be identified.'+
                             ' This means no growth parameters will be extracted'),
                    severity=Severity.warning)
                )
            return

        noBackground={}
        for tc in self.nonBackgroundWells():
            if tc.background is None:
                if tc.condition not in noBackground:
                    noBackground[tc.condition]={}
                if tc.sampleid not in noBackground[tc.condition]:
                    noBackground[tc.condition][tc.sampleid]=[]
                noBackground[tc.condition][tc.sampleid].append(tc)

        for tc in self.nonBackgroundReplicates():
            if tc.background is None:
                if tc.condition not in noBackground:
                    noBackground[tc.condition]={}
                if tc.sampleid not in noBackground[tc.condition]:
                    noBackground[tc.condition][tc.sampleid]=[]
                noBackground[tc.condition][tc.sampleid].append(tc)

        if len(noBackground.keys()):
            affected=''
            for condition in sorted(noBackground):
                if condition is None or condition == '':
                    affected+='no condition:'
                else:
                    affected+=condition+':'
                for sampleid in sorted(noBackground[condition]):
                    affected+=' '+sampleid
                affected+='\n'
            self._loadStatus.addStatus(
                StatusMessage(
                    key='No background for some samples:',shortmsg='plateinit:noBackgroundForSomeSamples',
                    longmsg=('For some conditions no background (blank) could be identified.'+
                             ' This means no growth parameters will be extracted. The affected samples are:\n'+
                             affected),
                    severity=Severity.warning)
                )

    def backgroundReplicateGroupIndices(self):
        """
        Return indices into self.replicateGroups for replicate groups being listed as background.

        :return: list(int) -- indices of background replicate groups
        """
        if self._backgroundGroupIndices is None:
            self._setupBackgroundIndices()
        return self._backgroundGroupIndices

    def backgroundReplicateGroups(self):
        """
        Return replicate groups being listed as background.

        :return: list(Replicate) -- replicate groups listed as background
        """
        tcs=[]
        for idx in self.backgroundReplicateGroupIndices():
            tcs.append(self.replicateGroups[idx])
        return tcs

    def backgroundWellIndices(self):
        """
        Return indices into self.wells for wells being listed as background.

        :return: list(int) -- indices of background wells
        """
        if self._backgroundWellIndices is None:
            self._setupBackgroundIndices()
        return self._backgroundWellIndices

    def backgroundWells(self):
        """
        Return wells being listed as background.

        :return: list(Replicate) -- wells listed as background
        """
        tcs=[]
        for idx in self.backgroundWellIndices():
            tcs.append(self.wells[idx])
        return tcs

    def _createSampleConditionToWellIndices(self):
        """
        Create a mapping to quickly find single-well objects based on sample id and condition.

        For internal use only.
        """
        # gather sampleids and conditions
        self._conditionToWellIdx={}
        self._sampleConditionToWellIdcs={}
        tcidx=0
        for tc in self.wells:
            # add well to the condition mapping
            if tc.condition not in self._conditionToWellIdx:
                self._conditionToWellIdx[tc.condition]=[]
            self._conditionToWellIdx[tc.condition].append(tcidx)
    
            # add well to the replicate mapping (sampleid and condition)
            if tc.sampleid not in self._sampleConditionToWellIdcs:
                self._sampleConditionToWellIdcs[tc.sampleid]={}
            if tc.condition not in self._sampleConditionToWellIdcs[tc.sampleid]:
                self._sampleConditionToWellIdcs[tc.sampleid][tc.condition]=[]
            self._sampleConditionToWellIdcs[tc.sampleid][tc.condition].append(tcidx)
    
            tcidx+=1

    def _createReplicateGroupsFromSampleIdsNConditions(self):
        """
        Create replicate groups by grouping wells of the same sample id and condition.

        For internal use only.
        """
        if self._sampleConditionToWellIdcs is None:
            self._createSampleConditionToWellIndices()
        sampleids=list(self._sampleConditionToWellIdcs.keys())
        sampleids.sort()
        self.replicateGroups=[]
        for sampleid in sampleids:
            conditions=list(self._sampleConditionToWellIdcs[sampleid].keys())
            conditions.sort()
            for condition in conditions:
                comptc=Replicate(self,self._sampleConditionToWellIdcs[sampleid][condition],
                                                      None,condition,isReplicateGroup=True)
                self.replicateGroups.append(comptc)
                # set parental replicate group of the children
                for childtc in comptc.childWells():
                    childtc._setReplicateGroupParent(comptc)

    def _createSampleConditionToReplicateGroupIndices(self):
        """
        Create a mapping to quickly find replicate groups based on sample id and condition.

        For internal use only.
        """
        self._sampleConditionToReplicateGroupIdcs={}
        coidx=0
        for tc in self.replicateGroups:
            if tc.sampleid not in self._sampleConditionToReplicateGroupIdcs:
                self._sampleConditionToReplicateGroupIdcs[tc.sampleid]={}
            if tc.condition not in self._sampleConditionToReplicateGroupIdcs[tc.sampleid]:
                self._sampleConditionToReplicateGroupIdcs[tc.sampleid][tc.condition]=[]
            self._sampleConditionToReplicateGroupIdcs[tc.sampleid][tc.condition].append(coidx)
            coidx+=1

    def _createConditionToReplicateGroupIndices(self):
        """
        Create a mapping to quickly find all replicate groups for a specific condition.

        For internal use only.
        """
        self._conditionToReplicateGroupIdx={}
        coidx=0
        for tc in self.replicateGroups:
            # add replicate group to the condition mapping
            if tc.condition not in self._conditionToReplicateGroupIdx:
                self._conditionToReplicateGroupIdx[tc.condition]=[]
            self._conditionToReplicateGroupIdx[tc.condition].append(coidx)
            coidx+=1

    def _setBackgroundForAllReplicates(self,backgroundSampleIds):
        """
        Set background replicate group for single-wells and replicate groups.

        Currently, if there are multiple background ids, an exception is raised.

        For internal use only.
        """
        self._backgroundWellIndices=None
        self._backgroundGroupIndices=None

        if backgroundSampleIds is None or not len(backgroundSampleIds):
            if self.wells is not None:
                for tc in self.wells:
                    tc._setBackgroundIndex(None)
            if self.replicateGroups is not None:
                for tc in self.replicateGroups:
                    tc._setBackgroundIndex(None)
            self._setBackgroundStatus()
            return
        if len(backgroundSampleIds) > 1:
            raise Plate.MultipleBackgroundIdsError(backgroundSampleIds)
        backgroundSampleId=backgroundSampleIds[0]

        if self._sampleConditionToReplicateGroupIdcs is None:
            self._createSampleConditionToReplicateGroupIndices()

        # set background index for the single (non-averaged) wells
        for tc in self.wells:
            if tc.sampleid not in backgroundSampleIds:
                if tc.condition in self._sampleConditionToReplicateGroupIdcs[backgroundSampleId]:
                    # NOTE there should be only one element in self._sampleConditionToReplicateGroupIdcs[backgroundSampleId][tc.condition]
                    tc._setBackgroundIndex(self._sampleConditionToReplicateGroupIdcs[backgroundSampleId][tc.condition][0])

        # set background for replicate groups
        for tc in self.replicateGroups:
            if tc.sampleid not in backgroundSampleIds:
                if tc.condition in self._sampleConditionToReplicateGroupIdcs[backgroundSampleId]:
                    # NOTE there should be only one element in self._sampleConditionToReplicateGroupIdcs[backgroundSampleId][tc.condition]
                    tc._setBackgroundIndex(self._sampleConditionToReplicateGroupIdcs[backgroundSampleId][tc.condition][0])

        # append warnings to self._loadStatus if for some replicates no background was set
        self._setBackgroundStatus()

    def replicateGroupIdxForSampleCondition(self,sampleid,condition):
        """
        Return index of replicate group with the given sample Id and condition.

        :param sampleid: Id of the sample.
        :type sampleid: string
        :param condition: Condition under which the sample was grown.
        :type condition: string

        :return: int -- Index (into self.replicateGroups) of Replicate with given id and condition.
        """
        if self._sampleConditionToReplicateGroupIdcs is None:
            self._createSampleConditionToReplicateGroupIndices()
        if sampleid not in self._sampleConditionToReplicateGroupIdcs:
            return None
        if condition not in self._sampleConditionToReplicateGroupIdcs[sampleid]:
            return None
        if len(self._sampleConditionToReplicateGroupIdcs[sampleid][condition]) != 1:
            raise RuntimeError('more than one replicate group for '+sampleid+' '+condition)
        return self._sampleConditionToReplicateGroupIdcs[sampleid][condition][0]

    def replicateGroupForSampleCondition(self,sampleid,condition):
        """
        Return index of replicate group with the given sample Id and condition.

        :param sampleid: Id of the sample.
        :type sampleid: string
        :param condition: Condition under which the sample was grown.
        :type condition: string

        :return: Replicate -- replicate group with given id and condition.
        """
        idx=self.replicateGroupIdxForSampleCondition(sampleid,condition)
        if idx is None:
            return None
        return self.replicateGroups[idx]

    def replicateGroupIdcsForCondition(self,condition):
        """
        Return a list of indices of replicate groups with the given condition.

        :param condition: Condition under which the samples were grown.
        :type condition: string

        :return: list(int) -- Indices (into self.replicateGroups) of replicate groups with the given condition.
        """
        if self._conditionToReplicateGroupIdx is None:
            self._createConditionToReplicateGroupIndices()
        if condition not in self._conditionToReplicateGroupIdx:
            return None
        return self._conditionToReplicateGroupIdx[condition]

    def replicateGroupsForCondition(self,condition):
        """
        Return a list of replicate groups with the given condition.

        :param condition: Condition under which the samples were grown.
        :type condition: string

        :return: list(Replicate) -- Replicate groups with given condition.
        """
        idcs=self.replicateGroupIdcsForCondition(condition)
        if idcs is None:
            return None
        tcs=[]
        for idx in idcs:
            tcs.append(self.replicateGroups[idx])
        return tcs

    def conditions(self):
        """
        Return a list of conditions.

        :return: list(str) -- Conditions.
        """
        if self._conditionToReplicateGroupIdx is None:
            self._createConditionToReplicateGroupIndices()
        conditions=list(self._conditionToReplicateGroupIdx.keys())
        conditions.sort()
        return conditions

    def nonBackgroundReplicates(self):
        """
        :return: list(Replicate) -- replicate groups that are not background samples.
        """
        backgroundIndices=self.backgroundReplicateGroupIndices()
        nbckg=[]
        idx=0
        for tc in self.replicateGroups:
            if idx not in backgroundIndices:
                nbckg.append(tc)
            idx+=1
        return nbckg

    def nonBackgroundReplicateIndices(self):
        """
        :return: list(Replicate) -- Indices of replicate groups that are not background samples.
        """
        backgroundIndices=self.backgroundReplicateGroupIndices()
        nbckgidcs=[]
        idx=0
        for tc in self.replicateGroups:
            if idx not in backgroundIndices:
                nbckgidcs.append(idx)
            idx+=1
        return nbckgidcs

    def nonBackgroundWells(self):
        """
        :return: list(Replicate) -- wells that are not background samples.
        """
        backgroundIndices=self.backgroundWellIndices()
        nbckg=[]
        idx=0
        for tc in self.wells:
            if idx not in backgroundIndices:
                nbckg.append(tc)
            idx+=1
        return nbckg

    def _indexOfReplicateGroup(self,ctc):
        """
        Determine the index of the given replicate group.

        For internal use only.

        :return: int -- Index of replicate group.
        """

        if self.replicateGroups is None:
            return None

        idx=0
        idxOfTc=None
        for ttc in self.replicateGroups:
            if ttc._wellIndices == ctc._wellIndices:
                if idxOfTc is not None:
                    raise RuntimeError("multiple similar replicate groups?")
                else:
                    idxOfTc=idx
            idx+=1

        return idxOfTc

    def _parametersUpdated(self,par=None):
        """
        Notify replicate(s) that a parameter changed and memoised results should be deleted.

        For internal use only.

        :param par: The name of the parameter that was changed.
        :type par: str

        The Replicate objects memoise some results that are expensive
        to calculate. When a parameter is updated, the results may not
        be valid anymore and should get removed from the "cache".
        If par is given, this method can decide which results should
        be removed.
        """

        # only needed for non-background replicate groups (as background does not depend on parameters)
        for tc in self.nonBackgroundWells():
            tc._parametersUpdated(par,dontRecurse=True)
        for tc in self.nonBackgroundReplicates():
            tc._parametersUpdated(par,dontRecurse=True)
        self.modified=True

    def _replicateChanged(self,tc,par=None):
        """
        Update replicates that depend on the given replicate.

        For internal use only.
        """

        if self.replicateGroups is None:
            # for startup code: there are no replicate groups yet
            return

        idxOfTc=self._indexOfReplicateGroup(tc)

        if idxOfTc is None:
            raise RuntimeError("no matching tc for "+tc.fullId())

        for ptc in self.wells:
            if ptc._backgroundIndex == idxOfTc:
                ptc._parametersUpdated(par='backgroundRawOd')
        for ctc in self.replicateGroups:
            if ctc._backgroundIndex == idxOfTc:
                ctc._parametersUpdated(par='backgroundRawOd')

    def _getDefaultParameter(self,par):
        """
        Get default value of parameter.

        For internal use only.

        :param par: The name of the parameter.
        :type par: str

        The Plate stores values of plate-wide parameters
        and default parameters.
        """
        if par not in self._inheritableParameters:
            raise RuntimeError('_getDefaultParameter: unknown parameter '+par)
        return self._inheritableParameters[par]

    def _getExplicitParameter(self,par):
        """
        Get explicit value of parameter (alias for _getDefaultParameter).

        For internal use only.

        :param par: The name of the parameter.
        :type par: str
        """
        return self._getDefaultParameter(par)

    def getParameter(self,par):
        """
        Return the requested parameter.

        :param par: The name of the parameter.
        :type par: str

        If the parameter is explicitly set for the plate, this value
        returned. Otherwise return None.

        See chapter :ref:`parameters <gat parameters>` for details of
        parameter handling and available parameters.
        """
        return self._getDefaultParameter(par)

    def parameterIsEditible(self,par):
        """
        Return True if this is a parameter can have a plate-wide default.

        :return: bool -- True if parameter can be edited.

        Some parameters can only be changed per Replicate, some only
        per Plate. This method is used to distinguish between them.

        See chapter :ref:`parameters <gat parameters>` for details of
        parameter handling and available parameters.
        """
        if par in Plate._isNotPlateParameter and Plate._isNotPlateParameter[par]:
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
        is the case. Since this object is a plate it tells whether a
        default value has been set.

        See chapter :ref:`parameters <gat parameters>` for details of
        parameter handling and available parameters.
        """
        return self._getExplicitParameter(par) is not None

    def activeChildReplicatesHaveExplicitParameter(self,par):
        """
        Return True if for at least one of the replicate groups the given parameter is explicitly set.

        :param par: The name of the parameter.
        :type par: str

        :return: bool -- True if parameter is explicitly set in one of the replicate groups.

        See chapter :ref:`parameters <gat parameters>` for details of
        parameter handling and available parameters.
        """
        for childtc in self.nonBackgroundReplicates():
            if childtc._getExplicitParameter(par) is not None:
                return True
            if childtc.activeChildReplicatesHaveExplicitParameter(par):
                return True
        return False

    def _setDefaultParameter(self,par,val):
        """
        Change the (default) value of the given parameter.

        For internal use only.

        :param par: The name of the parameter that will be changed.
        :type par: str
        :param val: The new value.
        """
        if par not in self._inheritableParameters:
            raise RuntimeError('_setDefaultParameter: unknown parameter '+par)
        self._inheritableParameters[par]=val
        self._parametersUpdated(par)

    def _setExplicitParameter(self,par,val):
        """
        Change the value of the given parameter (alias for _setDefaultParameter).

        For internal use only.

        :param par: The name of the parameter that will be changed.
        :type par: str
        :param val: The new value.
        """
        self._setDefaultParameter(par,val)

    def setMaxGrowthLowerTimeCutoff(self,t):
        """Set lower limit of interval in which the maximal growth should be searched."""
        self._setDefaultParameter('maxGrowthLowerTimeCutoff',t)

    def setMaxGrowthUpperTimeCutoff(self,t):
        """Set upper limit of interval in which the maximal growth should be searched."""
        self._setDefaultParameter('maxGrowthUpperTimeCutoff',t)

    def setLogOdCutoff(self,lod):
        """Set cutoff value of log(OD)."""
        self._setDefaultParameter('logOdCutoff',lod)

    def setLagAtLogOdEquals(self,lagat):
        """Set value of log(OD) used to define the lag time."""
        self._setDefaultParameter('lagAtLogOdEquals',lagat)

    def setHighDensityCorrectionLinear(self,hdCorrectionLinear=None):
        """Set coefficient of linear term of high density correction."""
        self._setDefaultParameter('hdCorrectionLinear',hdCorrectionLinear)

    def setHighDensityCorrectionQuadratic(self,hdCorrectionQuadratic=None):
        """Set coefficient of quadratic term of high density correction."""
        self._setDefaultParameter('hdCorrectionQuadratic',hdCorrectionQuadratic)

    def setHighDensityCorrectionCubic(self,hdCorrectionCubic=None):
        """Set coefficient of cubic term of high density correction."""
        self._setDefaultParameter('hdCorrectionCubic',hdCorrectionCubic)

    def setSmoothingK(self,k):
        """Set degree of the smoothing spline."""
        self._setDefaultParameter('smoothingK',k)

    def setSmoothingS(self,s):
        """Set smoothing factor used to choose the number of knots."""
        self._setDefaultParameter('smoothingS',s)

    def setSlidingWindowSize(self,win):
        """
        Set number of datapoints of sliding windows.
        
        The value that is used for local exponential fit (growth rate) and linear regression (growth yield).
        """
        self._setDefaultParameter('slidingWindowSize',win)

    @staticmethod
    def guessWellIds(numberOfWells):
        """
        Return well ids by guessing the plate layout based on number of wells.

        This function will return A1-P24 or A1-H12.

        :param numberOfWells: number of wells of the plate
        :type numberOfWells: int

        :return: list(str) -- the guessed well ids (None if layout could not be guessed)
        """
        # some "heuristics" about well ids: A1-P24 or A1-H12
        if numberOfWells == 384:
            labeldivisor=24
        elif numberOfWells == 96:
            labeldivisor=12
        else:
            return None
        rowlabels=[chr(x) for x in range(ord('A'), ord('P') + 1)]
        wellids=[]
        for i in range(numberOfWells):
            (lblchar,lblnum)=divmod(i, labeldivisor)
            wellids.append(str(rowlabels[lblchar])+str(lblnum+1))
        return wellids

    @staticmethod
    def availableColumnsForCsvExport(logOdDerivativeProperties=True):
        """
        List the available properties that can be chosen for csv export.

        :param logOdDerivativeProperties: include properties determined from log(OD) derivative
        :type logOdDerivativeProperties: bool

        :return: list(str), list(str) -- fixed columns (ids), properties

        The 'fixed columns' list contains the sample/condition tuples
        which should always be exported in order to identify the
        replicates. For the other properties (except 'wellids') the
        variance can be chosen by adding '_var' to the property name.
        """
        fixedcolumns=['sample','condition']
        columns=[]
        columns.extend(['slope_linear',
                        'intercept_linear',
                        'timeOfMax_linear',
                        'lag_linear'])
        columns.extend(['doublingtime_expfit',
                        'growthrate_expfit',
                        'od0_expfit',
                        'timeOfMax_expfit',
                        'lag_expfit'])
        if logOdDerivativeProperties:
            columns.extend(['doublingtime_local',
                            'growthrate_local',
                            'od0_local',
                            'timeOfMax_local',
                            'lag_local'])
        columns.extend(['yield',
                        'timeOfYield'])
        columns.extend(['wellids'])
        return fixedcolumns, columns

    def growthParametersToCsv(self,filename,addVarianceColumns=True,singleWells=False, columns=None, progressCall=None,
                              **csvkwargs):
        """
        Write a "comma seperated values" (csv) file of properties for all replicate groups.

        :param filename: Filename.
        :type filename: string
        :param columns: List of properties that shall get exported (in that order).
        :type columns: list(str)
        :param addVarianceColumns: For each entry in columns add the corresponding variance
        :type addVarianceColumns: bool
        :param singleWells: Export properties of single well replicates instead of replicate groups
        :type singleWells: bool
        :param progressCall: Function that will be called on each iteration.
        :type progressCall: @fun(int)
        :param csvkwargs: Parameters which are passed on to the csv module; defaults to { 'dialect': 'excel' }
        :type csvkwargs: dict()
        """
        if 'dialect' not in csvkwargs:
            csvkwargs['dialect']='excel'

        col2collabel={
            'lag_expfit': 'lag_expfit (ln(OD) == lagAtCutoff)',
            'lag_expfit_var': 'lag_expfit_var (ln(OD) == lagAtCutoff)',
            'lag_local': 'lag_local (ln(OD) == lagAtCutoff)',
            'lag_local_var': 'lag_local_var (ln(OD) == lagAtCutoff)',
            }

        if columns is None:
            columns, morecolumns=Plate.availableColumnsForCsvExport()
            columns.extend(morecolumns)
        if addVarianceColumns and not singleWells:
            newcolums=[]
            for col in columns:
                newcolums.append(col)
                if col in ['sample','condition','wellids']:
                    continue
                if not col.endswith('_var') and col+'_var' not in columns:
                    newcolums.append(col+'_var')
            columns=newcolums
        if singleWells:
            replicates=self.nonBackgroundWells()
        else:
            replicates=self.nonBackgroundReplicates()

        with CsvFileUnicodeWriter(filename,**csvkwargs) as sliwriter:
            descrow=[]
            for col in columns:
                if col in col2collabel:
                    descrow.append(col2collabel[col])
                else:
                    descrow.append(col)
            sliwriter.writerow(descrow)
    
            allcnt=-1
            for tc in replicates:
                allcnt+=1
                if progressCall is not None:
                    progressCall(allcnt)

                if tc.od() is not None:
                    doublingtime_ef=None
                    doublingtimevar_ef=None
                    doublingtime_nls=None
                    doublingtimevar_nls=None
                    lag_linear=None
                    lagVar_linear=None
                    mu_ef, mu_ef_var, od0_ef, od0_ef_var, maxt_ef, maxt_ef_var, lag_ef, lag_ef_var, method_ef, status = tc.maxGrowthrate()
                    mu_nls, mu_nls_var, od0_nls, od0_nls_var, maxt_nls, maxt_nls_var, lag_nls, lag_nls_var, method_nls, status = tc.maxGrowthrateFromLogOdDerivative()
                    growthyield, growthyield_var, tgrowthyield, tgrowthyield_var, status=tc.growthyield()
                    slope_linear, slopeVar_linear, intercept_linear, interceptVar_linear, timeOfMax_linear, timeOfMaxVar_linear, timeOfMaxIndices_linear, plainSlopeStatus=tc.odSlopemaxIntercept()
                    doublingtime_ef, doublingtimevar_ef=Replicate.growthrateToDoublingTime(mu_ef,mu_ef_var)
                    doublingtime_nls, doublingtimevar_nls=Replicate.growthrateToDoublingTime(mu_nls,mu_nls_var)
                    if slope_linear is not None and slope_linear != 0:
                        lag_linear=-intercept_linear/(slope_linear)
                        if slopeVar_linear is not None and interceptVar_linear is not None:
                            lagVar_linear=((intercept_linear/(slope_linear**2))**2 * slopeVar_linear +
                                            1/slope_linear**2 * interceptVar_linear)
                else:
                    (doublingtime_ef, doublingtimevar_ef, doublingtime_nls, doublingtimevar_nls)=(None,None,None,None)
                    (mu_ef, mu_ef_var, od0_ef, od0_ef_var, maxt_ef, maxt_ef_var, lag_ef, lag_ef_var)=([None,None,None,None,None,None,None,None])
                    (mu_nls, mu_nls_var, od0_nls, od0_nls_var, maxt_nls, maxt_nls_var, lag_nls, lag_nls_var)=([None,None,None,None,None,None,None,None])
                    (growthyield,growthyield_var,tgrowthyield,tgrowthyield_var)=([None,None,None,None])
                    (slope_linear, slopeVar_linear, intercept_linear, interceptVar_linear,
                     timeOfMax_linear, timeOfMaxVar_linear, lag_linear, lagVar_linear)=([None,None,None,None,None,None,None,None])
    
                thisrow=[]
                for col in columns:
                    if col == 'sample':
                        thisrow.append(tc.sampleid)
                    elif col == 'condition':
                        thisrow.append(tc.condition)
                    elif col == 'slope_linear':
                        thisrow.append(slope_linear)
                    elif col == 'slope_linear_var':
                        thisrow.append(slopeVar_linear)
                    elif col == 'intercept_linear':
                        thisrow.append(intercept_linear)
                    elif col == 'intercept_linear_var':
                        thisrow.append(interceptVar_linear)
                    elif col == 'timeOfMax_linear':
                        thisrow.append(timeOfMax_linear)
                    elif col == 'timeOfMax_linear_var':
                        thisrow.append(timeOfMaxVar_linear)
                    elif col == 'lag_linear':
                        thisrow.append(lag_linear)
                    elif col == 'lag_linear_var':
                        thisrow.append(lagVar_linear)
                    elif col == 'doublingtime_expfit':
                        thisrow.append(doublingtime_ef)
                    elif col == 'doublingtime_expfit_var':
                        thisrow.append(doublingtimevar_ef)
                    elif col == 'growthrate_expfit':
                        thisrow.append(mu_ef)
                    elif col == 'growthrate_expfit_var':
                        thisrow.append(mu_ef_var)
                    elif col == 'od0_expfit':
                        thisrow.append(od0_ef)
                    elif col == 'od0_expfit_var':
                        thisrow.append(od0_ef_var)
                    elif col == 'timeOfMax_expfit':
                        thisrow.append(maxt_ef)
                    elif col == 'timeOfMax_expfit_var':
                        thisrow.append(maxt_ef_var)
                    elif col == 'lag_expfit':
                        thisrow.append(lag_ef)
                    elif col == 'lag_expfit_var':
                        thisrow.append(lag_ef_var)
                    elif col == 'doublingtime_local':
                        thisrow.append(doublingtime_nls)
                    elif col == 'doublingtime_local_var':
                        thisrow.append(doublingtimevar_nls)
                    elif col == 'growthrate_local':
                        thisrow.append(mu_nls)
                    elif col == 'growthrate_local_var':
                        thisrow.append(mu_nls_var)
                    elif col == 'od0_local':
                        thisrow.append(od0_nls)
                    elif col == 'od0_local_var':
                        thisrow.append(od0_nls_var)
                    elif col == 'timeOfMax_local':
                        thisrow.append(maxt_nls)
                    elif col == 'timeOfMax_local_var':
                        thisrow.append(maxt_nls_var)
                    elif col == 'lag_local':
                        thisrow.append(lag_nls)
                    elif col == 'lag_local_var':
                        thisrow.append(lag_nls_var)
                    elif col == 'yield':
                        thisrow.append(growthyield)
                    elif col == 'yield_var':
                        thisrow.append(growthyield_var)
                    elif col == 'timeOfYield':
                        thisrow.append(tgrowthyield)
                    elif col == 'timeOfYield_var':
                        thisrow.append(tgrowthyield_var)
                    elif col == 'wellids':
                        thisrow.append(tc.activeChildWellIdStr())
                    else:
                        raise RuntimeError('unknown property '+col)
    
                sliwriter.writerow(thisrow)

    def timeseriesToCsv(self,filename,
                        addVarianceColumns=True,
                        singleWells=False,
                        columns=None,
                        fullId=False,
                        progressCall=None,
                        **csvkwargs):
        """
        Write a "comma seperated values" (csv) file of time series for all replicate groups.

        :param filename: Filename.
        :type filename: string
        :param columns: List of time series that shall get exported for each replicate.
        :type columns: list(str)
        :param addVarianceColumns: For each entry in columns add the corresponding variance
        :type addVarianceColumns: bool
        :param singleWells: Export time series of single well replicates instead of replicate groups
        :type singleWells: bool
        :param fullId: Label the columns with the full id (including well ids) instead of "sample condition"
        :type fullId: bool
        :param progressCall: Function that will be called on each iteration.
        :type progressCall: @fun(int)
        :param csvkwargs: Parameters which are passed on to the csv module; defaults to { 'dialect': 'excel' }
        :type csvkwargs: dict()
        """
        if 'dialect' not in csvkwargs:
            csvkwargs['dialect']='excel'
        col2collabel={
            'od': 'OD',
            'od_var': 'var(OD)',
            'lnod': 'ln(OD)',
            }
        if columns is None:
            columns=['od']
        if addVarianceColumns and not singleWells:
            newcolums=[]
            for col in columns:
                newcolums.append(col)
                if col in ['lnod']:
                    continue
                if not col.endswith('_var') and col+'_var' not in columns:
                    newcolums.append(col+'_var')
            columns=newcolums
        if singleWells:
            replicates=self.nonBackgroundWells()
        else:
            replicates=self.nonBackgroundReplicates()

        with CsvFileUnicodeWriter(filename,**csvkwargs) as sliwriter:
            # header
            descrow=['t']
            for tc in replicates:
                for col in columns:
                    if col in col2collabel:
                        lbl=col2collabel[col]
                    else:
                        lbl=col
                    if fullId:
                        lbl+=' '+tc.fullId()
                    else:
                        lbl+=' '+tc.sampleid+' '+tc.condition
                        if singleWells:
                            lbl+=' '+tc.activeChildWellIdStr()
                    descrow.append(lbl)
            sliwriter.writerow(descrow)
            # data
            allcnt=-1
            for ti in range(len(self.time)):
                allcnt+=1
                if progressCall is not None:
                    progressCall(allcnt)
                thisrow=[]
                thisrow.append(self.time[ti])
                for tc in replicates:
                    for col in columns:
                        if col == 'od':
                            if tc.od() is not None:
                                thisrow.append(tc.od()[ti])
                            else:
                                thisrow.append(None)
                        elif col == 'od_var':
                            if tc.odVar() is not None:
                                thisrow.append(tc.odVar()[ti])
                            else:
                                thisrow.append(None)
                        elif col == 'lnod':
                            if tc.logOd() is not None:
                                thisrow.append(tc.logOd()[ti])
                            else:
                                thisrow.append(None)
                        else:
                            raise RuntimeError('unknown property '+col)
                sliwriter.writerow(thisrow)

    @staticmethod
    def _numWellsToFormatString(numWells):
        """
        Return a string uniquely identifying a plate format.

        NOTE this function is subject to change.
        """
        if numWells == 100:
            return '100honeycomb'
        elif numWells == 200:
            return '200honeycomb'
        return str(numWells)

    @staticmethod
    def writeMetadata(filename,metadata,metadataKeys,plateformat='96',**csvkwargs):
        """
        :param metadata: the metadata
        :type metadata: list(dict)
        """
        columnMajorOrder=False
        if plateformat == '96':
            if len(metadata) != 96:
                raise RuntimeError('metadata is not of length 96')
            numcols=12
            numrows=8
        elif plateformat == '384':
            if len(metadata) != 384:
                raise RuntimeError('metadata is not of length 384')
            numcols=24
            numrows=16
        elif plateformat == '200honeycomb':
            if len(metadata) != 200:
                raise RuntimeError('metadata is not of length 200')
            columnMajorOrder=True
            numcols=20 # number of columns in the layout of the exported metadata
            numrows=10 # number of rows

        if plateformat == '96' or plateformat == '384':
            rowlabels=[chr(x) for x in range(ord('A'), ord('A') + numrows)]
            collabels=[str(i+1) for i in range(numcols)]
        elif plateformat == '200honeycomb':
            rowlabels=[str(i) for i in range(1,numrows+1)]
            collabels=[str(i+1) for i in range(0,len(metadata),numrows)]
        else:
            raise RuntimeError('not implemented for format other than 96, 384 or 200 honeycomb')

        if columnMajorOrder:
            reordered=[]
            # transpose
            for rowidx in range(numrows):
                for colidx in range(numcols):
                    metentryidx=rowidx + colidx * numrows
                    reordered.append(metadata[metentryidx])
        else:
            reordered=metadata # keep order

        if 'dialect' not in csvkwargs:
            csvkwargs['dialect']='excel'
        with CsvFileUnicodeWriter(filename,**csvkwargs) as writer:
            # header: just the name of the metadata
            for key in metadataKeys:
                row=[key]
                writer.writerow(row)
                # the column ids
                row=['<>']
                row.extend(collabels)
                writer.writerow(row)
                # now the data, divided into rows of numcols columns
                colit=reordered.__iter__()
                for rowlab in rowlabels:
                    row=[rowlab]
                    for j in range(numcols):
                        thismeta=next(colit)
                        val=thismeta[key] if key in thismeta else None
                        row.append(val)
                    writer.writerow(row)
                # an empty row
                row=[]
                writer.writerow(row)

    @staticmethod
    def readMetadata(filename,plateformat='96',**csvkwargs):
        """
        Read metadata from a csv file.

        For each metadata key a table is read. The table should be laid
        out as according to the plate layout. To get a template, call
        writeMetadata(outfile,[{} for i in range(numOfColumns)],Plate.metadataKeys)
        """
        columnMajorOrder=False
        if plateformat == '96':
            numcols=12
            numrows=8
        elif plateformat == '384':
            numcols=24
            numrows=16
        elif plateformat == '200honeycomb':
            columnMajorOrder=True
            numcols=20 # number of columns in the layout of the exported metadata
            numrows=10 # number of rows

        if plateformat == '96' or plateformat == '384':
            rowlabels=[chr(x) for x in range(ord('A'), ord('A') + numrows)]
            collabels=[str(i+1) for i in range(numcols)]
        elif plateformat == '200honeycomb':
            rowlabels=[str(i) for i in range(1,numrows+1)]
            collabels=[str(i+1) for i in range(0,numcols*numrows,numrows)]
        else:
            raise RuntimeError('not implemented for format other than 96, 384 or 200 honeycomb')

        # initialise the metadata list
        metadata=[{} for i in range(numcols*numrows)]
        if 'dialect' not in csvkwargs:
            csvkwargs['dialect']='excel'
        with CsvFileUnicodeReader(filename,**csvkwargs) as odreader:
            nextlinemode='nada'
            metkey=None
            lineno=0
            for row in odreader:
                lineno+=1
                if nextlinemode == 'nada':
                    if len(row) == 0 or row[0] == '':
                        # skip empty row
                        continue
                    else:
                        metkey=row[0]
                        nextlinemode='starttable'
                elif nextlinemode == 'starttable':
                    if len(row) == 0:
                        raise Plate.BadMetadata('Row at start of table is empty',lineno,filename=filename)
                    if row[0] != '<>' or row[1] != '1' or len(row) != numcols+1:
                        raise Plate.BadMetadata('This does not look like the beginning of a table'+
                                                ', expected row[0] == "<>", row[1] == "1" and len(row) == '+str(numcols+1)+
                                                ', but got row[0]=="'+str(row[0])+'", row[1] == "'+str(row[1])+
                                                '" and len(row) == '+str(len(row)),
                                                lineno,filename=filename)
                    nextlinemode='intable'
                    rowcnt=0
                    metit=metadata.__iter__()
                elif nextlinemode == 'intable':
                    rowcnt+=1
                    if len(row) == 0:
                        raise Plate.BadMetadata('Row '+str(rowcnt)+' is empty',lineno,filename=filename)
                    if row[0].upper() != rowlabels[rowcnt-1]:
                        raise Plate.BadMetadata('Row '+str(rowcnt)+' does not start with '+rowlabels[rowcnt-1]+
                                                ' (found "'+row[0].upper()+'")',lineno,filename=filename)
                    row.pop(0)
                    numOfValsThisRow=len(row)
                    if numOfValsThisRow > numcols:
                        numOfValsThisRow=numcols
                    # read the columns of this row
                    colit=row.__iter__()
                    for i in range(numOfValsThisRow):
                        val=next(colit)
                        if val == '':
                            # map empty sting to None
                            val=None
                        metentry=next(metit)
                        metentry[metkey]=val
                    # if the last columns are empty, fill them up with None
                    for i in range(numcols-numOfValsThisRow):
                        metentry=next(metit)
                        metentry[metkey]=None
                    if rowcnt == numrows:
                        nextlinemode='nada'

        if columnMajorOrder:
            reordered=[]
            for colidx in range(numcols):
                for rowidx in range(numrows):
                    metentryidx=rowidx * numcols + colidx
                    reordered.append(metadata[metentryidx])
            metadata=reordered

        return metadata


    class Error(Exception):
        """Base class for exceptions in this module."""
        pass

    class MultipleBackgroundIdsError(Error):
        """Exception raised if there are different IDs for background wells."""

        def __init__(self, backgroundSampleIds):
            self._backgroundSampleIds = backgroundSampleIds

        def __str__(self):
            return str('multiple keys were found that could be background (blank) samples, make sure there is only one.' + 
                       '\nThe keys are:\n'+str(self._backgroundSampleIds))

    class UnknownFileFormat(Error):
        """Exception raised when an unsupported serialisation format is opened."""

        def __init__(self,filename,serFormat=None,serFormatVersion=None,detailedError=None):
            self.filename = filename
            self.serFormat = serFormat
            self.serFormatVersion = serFormatVersion
            self.detailedError = detailedError

        def __str__(self):
            if self.serFormat is not None:
                if self.serFormat.startswith('clsplates'):
                    message= 'You tried to open a Chronological Life Span (CLS) file ("'+self.filename+'"), please use the CLS analyser for this'
                else:
                    message = 'Unsupported file format "'+self.serFormat+'"'
                    if self.serFormatVersion is not None:
                        message += ' version "'+self.serFormatVersion+'"'
                    message += ' in file "'+self.filename+'"'
            else:
                message = 'Unsupported file format in file "'+self.filename+'"'
            if self.detailedError is not None:
                message+=': '+self.detailedError+'.'
            else:
                message+='.'
            return message

    class BadMetadata(Error):
        """Exception raised when an unsupported serialisation format is opened."""

        def __init__(self,detailedError=None,lineno=None,filename=None):
            self.detailedError = detailedError
            self.filename = filename
            self.lineno = lineno

        def __str__(self):
            message = self.detailedError
            if self.lineno is not None:
                message += ' around line '+str(self.lineno)
            if self.filename is not None:
                message += ' in file '+str(self.filename)
            return message
