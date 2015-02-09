"""
This module implements the :py:class:`Cls` class for CATHODE.

Chronological life span Analysis Tool for High-throughput Optical
Density Experiments (CATHODE) Cls class.
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
import re

import numpy
import json
import bz2
import csv

from platereader.plate import Plate
from platereader.clsreplicate import ClsReplicate
from platereader.statusmessage import StatusMessage, Severity
from platereader.csvunicode import CsvFileUnicodeWriter

class Cls(object):
    """
    """

    def __init__(self,files=None,days=None,serialisedFilename=None):
        if serialisedFilename is not None:
            self._loadLightweight(serialisedFilename)
        else:
            if files is None:
                raise Cls.Error('need plate-files (GAT) that should be loaded')
            if days is None:
                raise Cls.Error('need days, one corresponding to each plate')
            if len(files) != len(days):
                raise Cls.Error('need the same number of files and days (one corresponding to each plate)')
            self._loadGatFiles(files,numpy.array(days,dtype=float))
        self.modified=False

    def _loadLightweight(self,filename):
        with bz2.BZ2File(filename, 'r') as rfile:
            pickled=rfile.read().decode("utf-8")

        try:
            unpickled = json.loads(pickled)
        except ValueError as err:
            raise Cls.UnknownFileFormat(filename,detailedError=str(err))

        # restore activation of child wells
        self._deserialiseLightweight(unpickled,directory=os.path.dirname(filename),filename=filename)

    def reload(self):
        """
        Reload the underlying :py:class:`Plates <.Plate>`.

        The activation of the ClsReplicate groups is preserved.
        """
        # store activation of child wells
        pp=self._serialiseLightweight()

        # reload plates, restore activation of child wells
        self._deserialiseLightweight(pp)

    def _serialiseLightweight(self,directory=None):
        sr=dict()
        sr['format']='clsplates'
        sr['formatversion']='1' # this is an unsigned integer
        # make sure paths are not stored absolute but rather relative to the directory of the saved file
        if directory is None or not os.path.isdir(directory):
            directory='.' # may happen when file was given on commandline without a directory
        relativeFiles=[]
        for f in self.files:
            frel=os.path.relpath(os.path.abspath(f),directory)
            frel=re.sub(r'\\', "/", frel) # windows file separator is changed to unix separator
            relativeFiles.append(frel)
        sr['files']=relativeFiles
        sr['days']=self.days.tolist()
        sr['clsWells']=[]
        for clstc in self.clsWells:
            sr['clsWells'].append(clstc._serialiseLightweight())
        sr['clsReplicateGroups']=[]
        for clstc in self.clsReplicateGroups:
            sr['clsReplicateGroups'].append(clstc._serialiseLightweight())

        return sr

    def saveLightweight(self,filename):
        """
        Partially serialise CLS data and save to file.

        :param filename: Name of the file.
        :type filename: str

        :return: StatusMessage/None -- non-fatal notifications.

        Allows to preserve the state of activated wells. The plates
        are not saved, therefore their relative file location must not
        change for the loading to work.
        """
        sr=self._serialiseLightweight(os.path.dirname(filename))
        status = None
        pickled = json.dumps(sr)

        with bz2.BZ2File(filename, 'w') as wfile:
            wfile.write(pickled.encode('utf-8'))

        self.modified=False
        return status

    def _deserialiseLightweight(self,unpickled,directory=None,filename=None):
        """
        Unpickle data and read the GAT plate-files.

        When re-reading files we would like to preserve the state of
        activated wells. This method does exactly this: an initialised
        Cls object (basically an array of re-read
        Plate objects) is amended with data that was
        previously saved.

        For internal use only.
        """
        if 'format' not in unpickled:
            raise Cls.UnknownFileFormat(filename,detailedError='no "format" keyword found in serialisation')
        serFormatVersion=unpickled['formatversion'] if 'formatversion' in unpickled else 'undefined'
        if unpickled['format'] != 'clsplates' or serFormatVersion != '1':
            raise Cls.UnknownFileFormat(filename,serFormat=unpickled['format'],serFormatVersion=serFormatVersion)

        # rewrite GAT paths to absolute paths
        if directory is None or not os.path.isdir(directory):
            directory='.' # may happen when file was given on commandline without a directory
        absoluteFiles=[]
        for f in unpickled['files']:
            absoluteFiles.append(os.path.abspath(directory+"/"+f))

        # load plates
        self._loadGatFiles(absoluteFiles,numpy.array(unpickled['days'],dtype=float))

        # check that GAT files have the same amount of wells/replicateGroups as this lightweight serialisation
        if len(unpickled['clsWells']) != len(self.clsWells):
            raise Cls.Error('number of cls wells different: '+len(self.clsWells)+' != '+len(unpickled['clsWells']))
        if len(unpickled['clsReplicateGroups']) != len(self.clsReplicateGroups):
            raise Cls.Error('number of cls replicate groups different: '+
                                  len(self.clsReplicateGroups)+' != '+len(unpickled['clsReplicateGroups']))
        clsidx=0
        for clstc in self.clsWells:
            clstc._deserialiseLightweight(unpickled['clsWells'][clsidx])
            clsidx+=1
        clsidx=0
        for clstc in self.clsReplicateGroups:
            clstc._deserialiseLightweight(unpickled['clsReplicateGroups'][clsidx])
            clsidx+=1

    def _loadGatFiles(self,files,days):
        """
        Reads data from GAT files and sets up 

        For internal use only.
        """
        self.days = days
        self.files = files

        self.plates=[]
        missingfiles=[]
        for f in self.files:
            if not os.path.exists(f):
                missingfiles.append(f)
        if len(missingfiles) > 0:
            raise Cls.PlateFileDoesNotExist(missingfiles)
        for f in self.files:
            self.plates.append(Plate(filename=f,fileformat='gat'))
            self.plates[-1].verbose=False

        self.clsWells=[]
        # we use the sampleid/conditions tuples from the first plate to initialise this cls analyser
        wellidx=0
        for tc0 in self.plates[0].wells:
            tc0fullId=tc0.fullId()
            listOfWells=[tc0]
            statuses=StatusMessage()
            for tcidx in range(1,len(self.plates)):
                listOfWells.append(self.plates[tcidx].wells[wellidx])
                if listOfWells[-1].fullId() != tc0fullId:
                    wellidstr=listOfWells[-1].activeChildWellIdStr()
                    raise Cls.Error('Plates do not match: could not find aquivalent to "'+tc0.fullId()+'" in "'
                                          +self.files[tcidx]+'" in well '+wellidstr+'.')
            self.clsWells.append(ClsReplicate(parent=self,listOfOdReplicates=listOfWells,days=self.days,
                                                                wellIndices=[wellidx],
                                                                wellids=tc0.wellids,sampleid=tc0.sampleid,condition=tc0.condition,
                                                                status=statuses))
            wellidx+=1

        self.clsReplicateGroups=[]
        # we use the sampleid/conditions tuples from the first plate to initialise this cls analyser
        for tc0 in self.plates[0].replicateGroups:
            listOfWells=[tc0]
            tc0fullId=tc0.fullId()
            statuses=StatusMessage()
            for tcidx in range(1,len(self.plates)):
                listOfWells.append(self.plates[tcidx].replicateGroupForSampleCondition(tc0.sampleid,tc0.condition))
                if listOfWells[-1] is None:
                    raise Cls.Error('Plates do not match: could not find aquivalent to "'+tc0.fullId()+'" in "'
                                          +self.files[tcidx]+'".')
                else:
                    diffstatus=Cls.differencesInChildActivations(tc0,listOfWells[-1],0,tcidx)
                    if diffstatus is not None:
                        statuses.addStatus(diffstatus)
            self.clsReplicateGroups.append(ClsReplicate(parent=self,listOfOdReplicates=listOfWells,days=self.days,
                                                                         status=statuses,
                                                                         wellIndices=tc0.childWellIndices(),isReplicateGroup=True,
                                                                         wellids=tc0.wellids,sampleid=tc0.sampleid,condition=tc0.condition))

        # NOTE we create a reference here in order to be able to use this class with ODplateModel
        self.replicateGroups=self.clsReplicateGroups

    def clsReplicateGroupIdxForSampleCondition(self,sampleid,condition):
        idx=self.plates[0].replicateGroupIdxForSampleCondition(sampleid,condition)
        if idx is None:
            return None
        return self.clsReplicateGroups[idx]

    def clsReplicateGroupIdcsForCondition(self,condition):
        idcs=self.plates[0].replicateGroupIdcsForCondition(condition)
        if idcs is None:
            return None
        tcs=[]
        for idx in idcs:
            tcs.append(self.clsReplicateGroups[idx])
        return tcs

    def conditions(self):
        return self.plates[0].conditions()

    def nonBackgroundClsIndices(self):
        """
        :return: list(ClsReplicate) -- ClsReplicate that are not background samples.
        """
        return self.plates[0].nonBackgroundReplicateIndices()

    def nonBackgroundCls(self):
        """
        :return: list(ClsReplicate) -- Indices of ClsReplicate that are not background samples.
        """
        nbckg=[]
        for idx in self.nonBackgroundClsIndices():
            nbckg.append(self.clsReplicateGroups[idx])
        return nbckg

    def numberOfNonBackgroundCls(self):
        return len(self.nonBackgroundClsIndices())

    def viabilities(self,sampleConditionTuples):
        if sampleConditionTuples is None:
            raise RuntimeError('sampleConditionTuples is None')

        viab=[]
        viabvar=[]
        for sct in sampleConditionTuples:
            sampleid, condition = sct
            clstc=self.clsReplicateGroupIdxForSampleCondition(sampleid,condition)
            dd, v, vvar, status = clstc.viability()
            viab.append(v)
            viabvar.append(vvar)

        return self.days,viab,viabvar

    def survivalToCsv(self,filename,showViabilities=True,columns=None,progressCall=None,
                      **csvkwargs):
        """
        :param csvkwargs: Parameters which are passed on to the csv module; defaults to { 'dialect': 'excel' }
        :type csvkwargs: dict()
        """
        if 'dialect' not in csvkwargs:
            csvkwargs['dialect']='excel'

        col2collabel={}

        with CsvFileUnicodeWriter(filename,**csvkwargs) as sliwriter:
            if columns is None:
                columns=['sample','condition','survivalIntegral','survivalIntegral_var']
                if showViabilities:
                    for d in self.days:
                        columns.extend(['viabilityDay%02d' % (d),'viabilityDay%02d_var' % (d)])
                columns.extend(['wellids'])
    
            descrow=[]
            for col in columns:
                if col in col2collabel:
                    descrow.append(col2collabel[col])
                else:
                    descrow.append(col)
            sliwriter.writerow(descrow)
    
            clsidx=-1
            for clstc in self.nonBackgroundCls():
                clsidx+=1
                if progressCall is not None:
                    progressCall(clsidx)

                si, sivar, status = clstc.survivalIntegral()
                valdict={}
                if showViabilities is not None:
                    days, viability, viabilityvar, status = clstc.viability()
                    idx=0
                    for d in self.days:
                        valdict['viabilityDay%02d' % (d)]=viability[idx] if viability is not None else None
                        valdict['viabilityDay%02d_var' % (d)]=viabilityvar[idx] if viabilityvar is not None else None
                        idx+=1
    
                thisrow=[]
                for col in columns:
                    if col == 'sample':
                        thisrow.append(clstc.sampleid)
                    elif col == 'condition':
                        thisrow.append(clstc.condition)
                    elif col == 'survivalIntegral':
                        thisrow.append(si)
                    elif col == 'survivalIntegral_var':
                        thisrow.append(sivar)
                    elif col == 'wellids':
                        thisrow.append(clstc.activeChildWellIdStr())
                    elif col in valdict:
                        thisrow.append(valdict[col])
                    else:
                        raise RuntimeError('unknown column '+col)
    
                sliwriter.writerow(thisrow)

    @staticmethod
    def differencesInChildActivations(tc1,tc2,tc1idx,tc2idx):
        """
        Return differences in sample/condition of wells of two plates.

        :return: StatusMessage -- None if no differences, otherwise a status message with the differences.
        """
        if tc1.sampleid != tc2.sampleid:
            return StatusMessage('diffSampleId',longmsg='different ids '+tc1.sampleid+' '+tc2.sampleid,severity=Severity.failed)
        if tc1.condition != tc2.condition:
            return StatusMessage('diffCond',longmsg='different conditions '+tc1.condition+' '+tc2.condition,severity=Severity.failed)
        if tc1.childWellIndices() != tc2.childWellIndices():
            msg=str(tc1idx)+'|'+str(tc2idx)+':\tdifferent child indices'
            only1=set(tc1.childWellIndices())-set(tc2.childWellIndices())
            only2=set(tc2.childWellIndices())-set(tc1.childWellIndices())
            if len(only1):
                msg+='  < '+str(only1)
            if len(only2):
                msg+='  > '+str(only2)
            return StatusMessage('diffChildIdcs'+str(tc1idx)+'|'+str(tc2idx),longmsg=msg,severity=Severity.message)
        if tc1.activeChildWellIndices() != tc2.activeChildWellIndices():
            msg=str(tc1idx)+'|'+str(tc2idx)+': different active children in underlying plates:'
            only1=set(tc1.activeChildWellIndices())-set(tc2.activeChildWellIndices())
            only2=set(tc2.activeChildWellIndices())-set(tc1.activeChildWellIndices())
            if len(only1):
                msg+='  < '+str(only1)
            if len(only2):
                msg+='  > '+str(only2)
            return StatusMessage('diffActChild'+str(tc1idx)+'|'+str(tc2idx),longmsg=msg,severity=Severity.message)

        return None

    class Error(Exception):
        """Base class for exceptions in this module."""
        pass

    class PlateFileDoesNotExist(Error):
        """
        Exception raised if a plate does not exist (anymore).

        Cls does not store the data in itself, but rather uses
        existing Plate files. If these got (re)moved we
        cannot initialise Cls.
        """

        def __init__(self, missingfiles):
            self.missingfiles = missingfiles

        def __str__(self):
            return str('Plate files that are needed for Chronological Life Span analysis are missing (maybe they were relocated?).' + 
                       ' Missing are:\n'+'\n'.join(self.missingfiles))

    class UnknownFileFormat(Error):
        """Exception raised when an unsupported serialisation format is opened."""

        def __init__(self,filename,serFormat=None,serFormatVersion=None,detailedError=None):
            self.filename = filename
            self.serFormat = serFormat
            self.serFormatVersion = serFormatVersion
            self.detailedError = detailedError

        def __str__(self):
            if self.serFormat is not None:
                if self.serFormat.startswith('opticaldensityplate'):
                    message = ('You tried to open an OD plate file ("'+self.filename
                                     +'"), please open a CLS file or create a new project with the plate')
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
