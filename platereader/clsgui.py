#! /usr/bin/env python
"""
This module implements a GUI for CATHODE.

Chronological life span Analysis Tool for High-throughput Optical
Density Experiments (CATHODE) GUI classes and entry point for scripts
using the GUI.
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
import sys
import re
import math
import traceback
import argparse
import numpy
import csv

import sip
sip.setapi('QVariant', 2)
sip.setapi('QString', 1)
from PyQt4 import Qt
from PyQt4 import QtCore
from PyQt4 import QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QTAgg as NavigationToolbar

from platereader.odgui import EmptyModel, ODplateItem, ODplateModel, MyMplCanvas, DeselectableTreeView
from platereader.odgui import getSaveFileNameDialogWithDefaultSuffix, MyProgressWidget, agpllicense

from platereader.plate import Plate
from platereader.statusmessage import StatusMessage, Severity
from platereader.cls import Cls
from platereader.clsreplicate import ClsReplicate
import platereader.clsplot
from platereader._version import __version__

########################################################################################################
def storeViewExpandState(view,model):
    """
    a very simple way to store the state of expanded ("opened") items in QTreeView
    """
    if model is None:
        return
    opened=[]
    for i in range(0,model.rowCount(view.rootIndex())):
        if view.isExpanded(model.index(i, 0, view.rootIndex())):
            opened.append(i)

    return opened

def restoreViewExpandState(view,model,opened):
    """
    restore expanded ("opened") items in QTreeView
    """
    if opened is None:
        return
    for i in opened:
        view.setExpanded(model.index(i, 0, view.rootIndex()), True)

########################################################################################################
class ClsWidgets(object):
    """
    Plot widgets, dialogs an the like needed for Chronological Life Span analysis
    """

    def __init__(self):
        self.tableModel=None
        self.clsView=DeselectableTreeView()
        self.mplcnvs=MyMplCanvas()
        self.mainstack=QtGui.QStackedWidget()
        self.mainstack.addWidget(self.mplcnvs)
        self.navigationWidget=QtGui.QWidget()
        self.ntb=NavigationToolbar(self.mplcnvs,self.navigationWidget)
        self.statusbar=None
        self.statusbarWidgets=[]
        self.clear()

    def clear(self):
        self.mplcnvs.fig.clear()
        self.mplcnvs.draw()

        self.statusToStatusBar(None)

        self.clsView.setModel(EmptyModel())
        self.tableModel=None

        self.platefiles=[]
        self.days=[]
        self.cls=None

    def disconnect(self):
        """
        NOTE this cannot be inside clear as disconnect fails if nothing is bound
        """
        if self.tableModel is not None:
            self.tableModel.dataChanged.disconnect()

    def setCls(self,cls,select={'row0': 0,'row1': None}):
        self.cls=cls
        self.tableModel=ODplateModel(self.cls)
        self.tableModel.dataChanged.connect(self.updatePlateFig)

        self.clsView.setModel(self.tableModel)
        self.clsView.header().setResizeMode(0, Qt.QHeaderView.ResizeToContents)
        self.clsView.header().setResizeMode(1, Qt.QHeaderView.ResizeToContents)
        self.clsView.selectionModel().selectionChanged.connect(self.selectionChanged)
        self.clsView.setFocus()
        self._setSample(self.cls)
        if select is not None:
            if select['row1'] is None:
                selectIndex=self.tableModel.index(select['row0'],0,self.clsView.currentIndex())
            else:
                parentIndex=self.tableModel.index(select['row0'],0,self.clsView.currentIndex())
                selectIndex=self.tableModel.index(select['row1'],0,parentIndex)
            self.clsView.setCurrentIndex(selectIndex)

    def loadSerialisedCls(self,catfilename):
        try:
            self.cls=Cls(None,None,serialisedFilename=catfilename)
            self.platefiles=self.cls.files
            self.days=self.cls.days
            self.setCls(self.cls)
        except Exception as err:
            self.cls=None
            self.platefiles=None
            self.days=None
            raise

    def loadCls(self,files,days,select={'row0': 0,'row1': None}):
        if files is None or days is None:
            return
        if len(files) != len(days):
            QtGui.QMessageBox.warning(Qt.qApp.activeWindow(),'Could not load plates',
                                      'Could not load plates because number of plates does not equal number of days')
        else:
            self.platefiles=files
            self.days=days
            self.cls=Cls(files,days)
            self.setCls(self.cls,select)

    def currentViewSelectionAsDict(self):
        sel=self.clsView.selectedIndexes()
        select=None
        if sel is not None and len(sel) > 0:
            if sel[0].parent().isValid():
                select={ 'row0': sel[0].parent().row(), 'row1': sel[0].row()}
            else:
                select={ 'row0': sel[0].row(), 'row1': None}
        return select

    def reload(self):
        # let mouse cursor appear busy
        QtGui.QApplication.setOverrideCursor(Qt.QCursor(QtCore.Qt.WaitCursor))

        # store selection of QTreeView
        select=self.currentViewSelectionAsDict()
        self.disconnect()
        # store state of QTreeView
        expState=storeViewExpandState(self.clsView,self.tableModel)

        # reload and restore selection
        self.cls.reload()
        self.setCls(self.cls,select)

        # restore state of QTreeView
        restoreViewExpandState(self.clsView,self.tableModel,expState)

        # mouse cursor should not be busy anymore
        QtGui.QApplication.restoreOverrideCursor()

    def selectionChanged(self,selected,deselected):
        """
        the model view selection changed
        """
        if len(selected.indexes()):
            selindex=selected.indexes()[0]
            tc=selindex.internalPointer()._itemData
            self._setSample(tc)
        else:
            self._setSample(self.cls)

    def _setSample(self,tc):
        self.selectedClsObject=tc
        self.updatePlateFig()

    def statusToStatusBar(self,status):
        if self.statusbar is None:
            return

        # clear previous statuses
        for st in self.statusbarWidgets:
            self.statusbar.removeWidget(st)
        self.statusbarWidgets=[]

        if status is None:
            return

        type2status=status.type2status()

        thekeys=list(type2status.keys())
        thekeys.sort()
        for sttype in thekeys:
            if type2status[sttype].severity() is Severity.warning:
                continue # handled below
            ql=QtGui.QLabel()
            ql.setText(type2status[sttype].shortmessage())
            ql.setToolTip(type2status[sttype].longmessage())
            self.statusbar.addPermanentWidget(ql)
            self.statusbarWidgets.append(ql)

        for sttype in thekeys:
            if type2status[sttype].severity() is not Severity.warning:
                continue # handled above
            ql=QtGui.QLabel()
            ql.setText(type2status[sttype].shortmessage())
            ql.setToolTip(type2status[sttype].longmessage())
            ql.setStyleSheet("color: red")
            self.statusbar.addPermanentWidget(ql)
            self.statusbarWidgets.append(ql)

    def updatePlateFig(self):
        self.mplcnvs.fig.clear()
        self.statusToStatusBar(None)
        if self.selectedClsObject is None:
            return
        QtGui.QApplication.setOverrideCursor(Qt.QCursor(QtCore.Qt.WaitCursor))

        try:
            status=StatusMessage()
            if isinstance(self.selectedClsObject, ClsReplicate):
                status=platereader.clsplot.viabilityToMatplotlib(self.selectedClsObject,self.mplcnvs.fig)
                self.mplcnvs.fig.tight_layout()
    
            self.mplcnvs.draw()
            self.statusToStatusBar(status)
            QtGui.QApplication.restoreOverrideCursor()
        except Exception as err:
            msg=str(err)+"\n"
            if msg is None:
                msg=''
            exc_type, exc_value, exc_tb = sys.exc_info()
            for line in traceback.format_exception(exc_type, exc_value, exc_tb):
                msg+=line
            self.mplcnvs.fig.clear()
            self.mplcnvs.draw()
            self.statusToStatusBar(None)
            QtGui.QApplication.restoreOverrideCursor()
            QtGui.QMessageBox.warning(self.mainstack,'Unknown error',msg)
            return

    class ClsFromPlateFilesDialog(QtGui.QDialog):
        def __init__(self,clswidgets,files=[],days=[]):
            super(ClsWidgets.ClsFromPlateFilesDialog, self).__init__()
            self.setWindowTitle('CLS files and days')

            self.clswidgets=clswidgets
            self.platefiles=files
            self.days=days

            self.lastSelectedDirectory=RememberLastVisistedDirectory('')
            #
            self.tableOfFileNDayLayout=QtGui.QVBoxLayout()
            self.setupTableOfFiles()

            self.addButtonLayout=QtGui.QHBoxLayout()
            self.addButton=QtGui.QPushButton("+")
            self.addButton.setMinimumWidth(40)
            self.addButton.setMinimumWidth(40)
            self.addButton.clicked.connect(self.addFileDay)
            self.addButtonLayout.addWidget(self.addButton)
            self.addButtonLayout.addStretch()

            self.cancelButton=QtGui.QPushButton("Cancel")
            self.okButton=QtGui.QPushButton("Ok")
            self.cancelButton.setMinimumWidth(80)
            self.cancelButton.setMaximumWidth(80)
            self.okButton.setMinimumWidth(80)
            self.okButton.setMaximumWidth(80)
            self.cancelButton.clicked.connect(self.close)
            self.okButton.clicked.connect(self.saveSettingsAndClose)
            self.okCancelLayout=QtGui.QHBoxLayout()
            self.okCancelLayout.addStretch()
            self.okCancelLayout.addWidget(self.cancelButton)
            self.okCancelLayout.addWidget(self.okButton)

            self.mainlayout=QtGui.QVBoxLayout()
            self.mainlayout.addLayout(self.tableOfFileNDayLayout)
            self.mainlayout.addLayout(self.addButtonLayout)
            self.mainlayout.addLayout(self.okCancelLayout)

            self.setLayout(self.mainlayout)

        def setupTableOfFiles(self):
            if self.platefiles is not None and self.days is not None and len(self.platefiles) == len(self.days):
                for idx in range(0,len(self.platefiles)):
                    self.tableOfFileNDayLayout.addLayout(FileLineEditChooserNDayEdit(self.lastSelectedDirectory,
                                                                                     self.platefiles[idx],
                                                                                     str(self.days[idx])))
                    if os.path.exists(self.platefiles[idx]):
                        self.lastSelectedDirectory.setDir(os.path.dirname(self.platefiles[idx]))
            self.tableOfFileNDayLayout.addLayout(FileLineEditChooserNDayEdit(self.lastSelectedDirectory))

        def addFileDay(self):
            self.tableOfFileNDayLayout.addLayout(FileLineEditChooserNDayEdit(self.lastSelectedDirectory))

        def saveSettingsAndClose(self):
            self.platefiles=[]
            self.days=[]
            for i in range(0,self.tableOfFileNDayLayout.count()):
                filename=str(self.tableOfFileNDayLayout.itemAt(i).fileLineEdit.existingFilename())
                day, ok=self.tableOfFileNDayLayout.itemAt(i).dayEdit.text().toDouble()
                if filename is not None and ok:
                    self.platefiles.append(filename)
                    self.days.append(day)
            self.accept()

    class PdfDialog(QtGui.QDialog):
        def __init__(self):
            super(ClsWidgets.PdfDialog, self).__init__()
            self.setWindowTitle('Save pdf options')

            self.currentOrAllRadios=QtGui.QButtonGroup()
            r1=QtGui.QRadioButton('Current selection')
            r2=QtGui.QRadioButton('All')
            r1.setChecked(True)
            self.currentOrAllRadios.addButton(r1)
            self.currentOrAllRadios.addButton(r2)
            self.currentOrAllRadios.setId(r1,0)
            self.currentOrAllRadios.setId(r2,1)

            currentOrAllGroupbox=QtGui.QGroupBox('CLS to export')
            grly=QtGui.QVBoxLayout()
            grly.addWidget(r1)
            grly.addWidget(r2)
            currentOrAllGroupbox.setLayout(grly)
            self.showTitle=QtGui.QCheckBox("Add title")
            self.showTitle.setToolTip("Add a title consisting of 'sample-condition' to figures.")
            showLayout=QtGui.QVBoxLayout()
            showLayout.addWidget(self.showTitle)
            showGroup=QtGui.QGroupBox('Figure options')
            showGroup.setLayout(showLayout)

            allLayout=QtGui.QHBoxLayout()
            allLayout.addWidget(currentOrAllGroupbox)
            allLayout.addWidget(showGroup)
            allGroup=QtGui.QGroupBox('')
            allGroup.setLayout(allLayout)

            buttonlayout=QtGui.QHBoxLayout()
            okButton=QtGui.QPushButton("Ok")
            okButton.setDefault(True)
            cancelButton=QtGui.QPushButton("Cancel")
            buttonlayout.addWidget(cancelButton)
            buttonlayout.addWidget(okButton)

            cancelButton.clicked.connect(self.reject)
            okButton.clicked.connect(self.accept)

            self.mainlayout=QtGui.QVBoxLayout()
            self.mainlayout.addWidget(allGroup)
            self.mainlayout.addLayout(buttonlayout)
            self.setLayout(self.mainlayout)

class FileLineEditChooserNDayEdit(QtGui.QHBoxLayout):
    def __init__(self,dirref=None,f=None,d=None):
        super(QtGui.QHBoxLayout, self).__init__()
        self.dirref=dirref
        self.fileLineEdit=FileLineEditChooser(self.dirref)
        self.fileLineEdit.setMinimumWidth(400)
        if f is not None:
            self.fileLineEdit.setFilename(f)
        self.dayEdit=Qt.QLineEdit()
        dayval=QtGui.QDoubleValidator()
        self.dayEdit.setValidator(dayval)
        self.dayEdit.setFixedWidth(40)
        if d is not None:
            self.dayEdit.setText(d)
        self.daysEdited()
        self.dayEdit.textEdited.connect(self.daysEdited)

        self.addWidget(self.fileLineEdit)
        self.addWidget(Qt.QLabel('Day:'))
        self.addWidget(self.dayEdit)

    def daysEdited(self):
        if self.dayEdit.hasAcceptableInput():
            self.dayEdit.setStyleSheet("")
        else:
            self.dayEdit.setStyleSheet("border: 2px solid red")

class FileLineEditChooser(QtGui.QWidget):
    def __init__(self,dirref=None):
        super(QtGui.QWidget, self).__init__()
        if dirref is not None:
            self.dirref=dirref
        else:
            self.dirref=RememberLastVisistedDirectory('')
        self.fedit=Qt.QLineEdit()
        self.filenameEdited()
        self.fedit.textEdited.connect(self.filenameEdited)
        self.browsebutton=QtGui.QPushButton("Browse")
        self.browsebutton.clicked.connect(self.chooseFile)

        self.mainlayout=QtGui.QHBoxLayout()
        self.mainlayout.addWidget(self.fedit)
        self.mainlayout.addWidget(self.browsebutton)

        self.setLayout(self.mainlayout)

    def filenameEdited(self):
        if os.path.exists(self.fedit.text()):
            self.fedit.setStyleSheet("")
        else:
            self.fedit.setStyleSheet("border: 2px solid red")

    def chooseFile(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Add file', self.dirref.dir(), 'Optical density plate (*.gat)')
        if os.path.exists(fname):
            self.fedit.setText(fname)
            self.dirref.setDir(os.path.dirname(str(fname)))
        self.filenameEdited()

    def setMinimumWidth(self,val):
        self.fedit.setMinimumWidth(val)

    def setFilename(self,filename):
        self.fedit.setText(filename)
        self.filenameEdited()

    def filename(self):
        return self.fedit.text()

    def existingFilename(self):
        if os.path.exists(self.fedit.text()):
            return self.fedit.text()
        return None

class RememberLastVisistedDirectory(object):
    def __init__(self,d=None):
        self.d=''
    def dir(self):
        return self.d
    def setDir(self,d):
        self.d=d

class ODMainWindow(QtGui.QMainWindow):
    def __init__(self):
        super(ODMainWindow, self).__init__()
        self.initUI()
        
    def initUI(self):
        self.clswidgets=ClsWidgets()
        self.clswidgets.statusbar=self.statusBar()

        self.newAction = QtGui.QAction('&New', self)
        self.newAction.setShortcut(QtGui.QKeySequence.New)
        self.newAction.setStatusTip('Create new project from plate files.')
        self.newAction.triggered.connect(self.newClsFromPlatefilesChoose)

        self.openAction = QtGui.QAction('&Open', self)
        self.openAction.setShortcut(QtGui.QKeySequence.Open)
        self.openAction.setStatusTip('Open file.')
        self.openAction.triggered.connect(self.openSerialisedClsChoose)

        self.reloadAction = QtGui.QAction('&Reload', self)
        self.reloadAction.setShortcut(QtGui.QKeySequence.Refresh)
        self.reloadAction.setStatusTip('Reload plate files.')
        self.reloadAction.triggered.connect(self.clswidgets.reload)

        self.saveAsAction = QtGui.QAction('Save as', self)
        self.saveAsAction.setShortcut(QtGui.QKeySequence.SaveAs)
        self.saveAsAction.setStatusTip('Save data.')
        self.saveAsAction.triggered.connect(self.saveAs)

        self.saveAction = QtGui.QAction('Save', self)
        self.saveAction.setShortcut(QtGui.QKeySequence.Save)
        self.saveAction.setStatusTip('Save data.')
        self.saveAction.triggered.connect(self.save)

        self.savePdfAction = QtGui.QAction('Export &figure', self)        
        self.savePdfAction.setStatusTip('Save a figure of viabilities/survival integrals.')
        self.savePdfAction.triggered.connect(self.savePdf)

        self.saveCsvAction = QtGui.QAction('Export &properties', self)        
        self.saveCsvAction.setStatusTip('Save a table containing viabilities/survival integrals.')
        self.saveCsvAction.triggered.connect(self.saveCsv)

        self.closeAction = QtGui.QAction('&Close', self)
        self.closeAction.setShortcut(QtGui.QKeySequence.Close)
        self.closeAction.setStatusTip('Close file.')
        self.closeAction.triggered.connect(self.close)

        self.exitAction = QtGui.QAction('&Exit', self)
        self.exitAction.setShortcut(QtGui.QKeySequence.Quit)
        self.exitAction.setStatusTip('Exit.')
        self.exitAction.triggered.connect(self.quit)

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(self.newAction)
        fileMenu.addAction(self.openAction)
        fileMenu.addAction(self.reloadAction)
        fileMenu.addAction(self.saveAsAction)
        fileMenu.addAction(self.saveAction)
        fileMenu.addAction(self.saveCsvAction)
        fileMenu.addAction(self.savePdfAction)
        fileMenu.addAction(self.closeAction)
        fileMenu.addAction(self.exitAction)

        self.aboutAction = QtGui.QAction('&About', self)        
        self.aboutAction.triggered.connect(self.showAbout)
        helpMenu = menubar.addMenu('&Help')
        helpMenu.addAction(self.aboutAction)


        self.mainwidget=QtGui.QWidget(self)
        mainwidgetLayout=QtGui.QVBoxLayout()
        self.mainwidget.setLayout(mainwidgetLayout)

        self.setCentralWidget(self.mainwidget)

        self.graphsnchooser = QtGui.QSplitter(QtCore.Qt.Horizontal,self.mainwidget)
        mainwidgetLayout.addWidget(self.graphsnchooser)
        self.graphsnchooser.setSizePolicy(Qt.QSizePolicy.Expanding,
                                          Qt.QSizePolicy.Expanding)

        self.graphsnchooser.addWidget(self.clswidgets.clsView)
        self.graphsnchooser.addWidget(self.clswidgets.mainstack)
        self.graphsnchooser.setStretchFactor(0,1)
        self.graphsnchooser.setStretchFactor(1,4)

        self.setGeometry(100, 100, 1024, 600)

        self.clear()
        self.show()

    def clear(self):
        self.newAction.setEnabled(True)
        self.openAction.setEnabled(True)
        self.reloadAction.setEnabled(False)
        self.saveAsAction.setEnabled(False)
        self.saveAction.setEnabled(False)
        self.saveCsvAction.setEnabled(False)
        self.savePdfAction.setEnabled(False)
        self.closeAction.setEnabled(False)
        self.exitAction.setEnabled(True)
        self.saveFilename=None
        self.clswidgets.clear()
        self._setWindowTitle()

    def enDisableActionsAfterLoading(self):
        self.newAction.setEnabled(False)
        self.openAction.setEnabled(False)
        self.reloadAction.setEnabled(True)
        self.saveAsAction.setEnabled(True)
        self.saveAction.setEnabled(True)
        self.saveCsvAction.setEnabled(True)
        self.savePdfAction.setEnabled(True)
        self.closeAction.setEnabled(True)
        self.exitAction.setEnabled(True)

    def _setWindowTitle(self):
        wtitle=''
        if self.clswidgets.cls is not None:
            if self.saveFilename is not None:
                wtitle=os.path.basename(self.saveFilename)+' - '
            elif (self.clswidgets.platefiles is not None and len(self.clswidgets.platefiles)
                  and self.clswidgets.days is not None and len(self.clswidgets.days)):
                wtitle=os.path.basename(self.clswidgets.platefiles[0])+' day '+'{:g}'.format(self.clswidgets.days[0])+', ... - '
        self.setWindowTitle(wtitle+'Chronological life span Analysis Tool (CATHODE)')

    def newClsFromPlatefilesChoose(self):
        dialog=ClsWidgets.ClsFromPlateFilesDialog(self.clswidgets,self.clswidgets.platefiles,self.clswidgets.days)
        if dialog.exec_() == QtGui.QDialog.Accepted:
            self.newClsFromPlatefiles(dialog.platefiles,dialog.days)

    def newClsFromPlatefiles(self,platefiles,days):
        try:
            self.clswidgets.loadCls(platefiles,days)
            self.enDisableActionsAfterLoading()
        except Cls.Error as err:
            QtGui.QMessageBox.warning(self,'Error loading file',str(err))
            return
        except Exception as err:
            msg=str(err)+"\n"
            if msg is None:
                msg=''
            exc_type, exc_value, exc_tb = sys.exc_info()
            for line in traceback.format_exception(exc_type, exc_value, exc_tb):
                msg+=line
            QtGui.QMessageBox.warning(self,'Unknown error',msg)
        self._setWindowTitle()

    def openSerialisedClsChoose(self):
        self.close()
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '', 'Cls data (*.cat)')
        if len(fname) == 0:
            return
        condition=""
        self.openSerialisedCls(str(fname))

    def openSerialisedCls(self,catfilename):
        self.close()
        if not os.path.exists(catfilename):
            QtGui.QMessageBox.warning(self,'Error loading file',catfilename+' does not exist.')
            self.clear()
            return
        try:
            self.clswidgets.loadSerialisedCls(catfilename)
            self.saveFilename=catfilename
            self.enDisableActionsAfterLoading()
        except (Cls.Error, Plate.Error) as err:
            QtGui.QMessageBox.warning(self,'Error loading file',str(err))
            return
        except Exception as err:
            msg=str(err)+"\n"
            if msg is None:
                msg=''
            exc_type, exc_value, exc_tb = sys.exc_info()
            for line in traceback.format_exception(exc_type, exc_value, exc_tb):
                msg+=line
            QtGui.QMessageBox.warning(self,'Unknown error',msg)
        self._setWindowTitle()

    def saveAs(self):
        if self.clswidgets.cls is None:
            return
        fname = str(getSaveFileNameDialogWithDefaultSuffix(self, 'Save file', '', 'CLS data (*.cat)', 'cat'))
        if len(fname) == 0:
            return
        self.saveFilename=fname
        self.saveAction.setEnabled(True)
        self._setWindowTitle()
        return self._save()

    def save(self):
        if self.saveFilename is None:
            return self.saveAs()
        else:
            return self._save()

    def _save(self):
        try:
            savestatus=self.clswidgets.cls.saveLightweight(self.saveFilename)
            if savestatus is not None:
                QtGui.QMessageBox.warning(self,savestatus.messageType(),savestatus.longmessage())
            return True
        except (Cls.Error, Plate.Error) as err:
            QtGui.QMessageBox.warning(self,'Error saving datafile',str(err))
            odplate_=None
        except Exception as err:
            msg=str(err)+"\n"
            if msg is None:
                msg=''
            exc_type, exc_value, exc_tb = sys.exc_info()
            for line in traceback.format_exception(exc_type, exc_value, exc_tb):
                msg+=line
            QtGui.QMessageBox.warning(self,'Unknown error',msg)
            odplate_=None
        return False

    def close(self):
        if self.clswidgets.cls is not None and self.clswidgets.cls.modified:
            ret=QtGui.QMessageBox.question(self,"CLS Analysis Tool (CATHODE)",
                                           "The document has been modified.\n"
                                           +"Do you want to save your changes?",
                                           QtGui.QMessageBox.Save | QtGui.QMessageBox.Discard | QtGui.QMessageBox.Cancel,
                                           QtGui.QMessageBox.Save)
            if ret == QtGui.QMessageBox.Save:
                if not self.save():
                    return False
            elif ret == QtGui.QMessageBox.Cancel:
                return False
            #elif ret == QtGui.QMessageBox.Discard:

        self.clear()
        return True

    def quit(self):
        if self.close():
            QtGui.qApp.quit()

    def saveCsv(self):
        if self.clswidgets.cls is None:
            return
        fname = str(getSaveFileNameDialogWithDefaultSuffix(self, 'Save file', '', 'csv (*.csv)', 'csv'))
        if len(fname) == 0:
            return

        # create progress bar, do not allow to press buttons, etc
        self._setupLongrunningCalculation('Saving table',len(self.clswidgets.cls.clsReplicateGroups))
        # run the real calculation
        try:
            self.clswidgets.cls.survivalToCsv(fname,progressCall=self._updateProgress)
        except Exception as err:
            exc_type, exc_value, exc_tb = sys.exc_info()
            msg=str(err)+"\n"
            if msg is None:
                msg=''
            for line in traceback.format_exception(exc_type, exc_value, exc_tb):
                msg+=line
            QtGui.QMessageBox.warning(self,'An error occurred',msg)
        # clean up: remove progress bar, enable buttons etc
        self._longrunningCalculationFinished()

    def savePdf(self):
        if self.clswidgets.cls is None:
            return
        # ask what to export: currently selected one or all CLS replicate groups
        dialog=ClsWidgets.PdfDialog()
        dialog.show()
        if not dialog.exec_():
            return
        # ask for filename
        fname = str(getSaveFileNameDialogWithDefaultSuffix(self, 'Save file', '', 'pdf (*.pdf)', 'pdf'))
        if len(fname) == 0:
            return

        if dialog.currentOrAllRadios.checkedId() == 0:
            # currently selected
            elementaryIndices=[]
            replicateGroupIndices=[]
            select=self.clswidgets.currentViewSelectionAsDict()
            addWellIdsToTitle=False
            if select['row1'] is None:
                replicateGroupIndices=[select['row0']]
            else:
                elementaryIndices=[self.clswidgets.cls.clsReplicateGroups[select['row0']].childWellIndices()[select['row1']]]
                addWellIdsToTitle=True
            platereader.clsplot.viabilitiesToPdf(self.clswidgets.cls,
                                                 fname,
                                                 showTitle=dialog.showTitle.checkState(),
                                                 addWellIdsToTitle=addWellIdsToTitle,
                                                 creator=os.path.basename(sys.argv[0]),
                                                 replicateGroupIndices=replicateGroupIndices,
                                                 elementaryIndices=elementaryIndices)
        elif dialog.currentOrAllRadios.checkedId() == 1:
            # all replicates
            # create progress bar, do not allow to press buttons, etc
            self._setupLongrunningCalculation('Saving figures',self.clswidgets.cls.numberOfNonBackgroundCls())
            try:
                platereader.clsplot.viabilitiesToPdf(self.clswidgets.cls,
                                                     fname,
                                                     showTitle=dialog.showTitle.checkState(),
                                                     creator=os.path.basename(sys.argv[0]),
                                                     replicateGroupIndices=[],
                                                     addWellIdsToTitle=False,progressCall=self._updateProgress)
            except Exception as err:
                exc_type, exc_value, exc_tb = sys.exc_info()
                msg=str(err)+"\n"
                if msg is None:
                    msg=''
                for line in traceback.format_exception(exc_type, exc_value, exc_tb):
                    msg+=line
                QtGui.QMessageBox.warning(self,'An error occurred',msg)
            # clean up: remove progress bar, enable buttons etc
            self._longrunningCalculationFinished()
        else:
            raise RuntimeError('PdfDialog: unknown selection')

    def _setupLongrunningCalculation(self,msg,numForProgress):
        self.setEnabled(False)
        QtGui.QApplication.setOverrideCursor(Qt.QCursor(QtCore.Qt.WaitCursor))
        # create a progress that is shown in the main stack (where usually the plot is shown)
        self.progresswidget=MyProgressWidget(msg,numForProgress)
        self.clswidgets.mainstack.addWidget(self.progresswidget)
        self.clswidgets.mainstack.setCurrentIndex(1)
        self.clswidgets.statusbar.setVisible(False)

    def _longrunningCalculationFinished(self):
        self.clswidgets.statusbar.setVisible(True)
        self.clswidgets.mainstack.setCurrentIndex(0)
        self.clswidgets.mainstack.removeWidget(self.progresswidget)
        self.progresswidget=None
        self.setEnabled(True)
        self.enDisableActionsAfterLoading()
        QtGui.QApplication.restoreOverrideCursor()

    def _updateProgress(self,intval):
        self.progresswidget.updateProgress(intval)
        QtGui.QApplication.processEvents()

    def showAbout(self):
        agpllicense=("This program is free software: you can redistribute it and/or modify "
                     +"it under the terms of the GNU Affero General Public License as "
                     +"published by the Free Software Foundation, either version 3 of the "
                     +"License, or (at your option) any later version.")
        msgBox=QtGui.QMessageBox.about(self,'Chronological life span Analysis Tool (CATHODE)',
                                       'Chronological life span Analysis Tool for High-throughput Optical Density Experiments (CATHODE)\nversion '
                                       +__version__
                                       +'\nCopyright (C) 2014 Nils Christian\n'+
                                       agpllicense)


########################################################################################################
########################################################################################################

def clsgui_main():
    commandline=' '.join(sys.argv)
    executablename = os.path.basename(sys.argv[0])
    parser = argparse.ArgumentParser(description=executablename+': GUI for Chronological Life Span analysis')
    parser.add_argument('--version', action='store_true', default=None, help='show version and exit')
    parser.add_argument('infiles', metavar='', nargs='*',
                        help='CATHODE file (.cat) or multiple GATHODE files (.gat)')
    parser.add_argument('--day', metavar='float', action='append', type=float, dest='days',
                        help='day of lag/growthrate measurement of the corresponding plate file')

    args = parser.parse_args()

    if args.version:
        print(executablename+' '+__version__)
        sys.exit(0)

    app = QtGui.QApplication(sys.argv)
    qod=ODMainWindow()

    if args.infiles:
        if args.infiles[0].endswith(".cat"):
            qod.openSerialisedCls(args.infiles[0])
        else:
            qod.newClsFromPlatefiles(args.infiles,args.days)

    sys.exit(app.exec_())

if __name__ == "__main__":
    clsgui_main()
