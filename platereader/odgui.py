#!/usr/bin/env python
"""
This module implements a GUI for GATHODE.

Growth Analysis Tool for High-throughput Optical Density Experiments
(GATHODE) GUI classes and entry point for scripts using the GUI.
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

import sys
import os.path
import traceback
import argparse

import sip
sip.setapi('QVariant', 2)
sip.setapi('QString', 1)
from PyQt4 import Qt
from PyQt4 import QtCore
from PyQt4 import QtGui

import matplotlib.figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

from platereader.plate import Plate
from platereader.odplot import replicateToFig, plotFullOdPlate, odPlateOverviewToAxes, replicateNderivativePlotAsPdf
from platereader.statusmessage import StatusMessage, Severity
from platereader._version import __version__

agpllicense=("This program is free software: you can redistribute it and/or modify "
             +"it under the terms of the GNU Affero General Public License as "
             +"published by the Free Software Foundation, either version 3 of the "
             +"License, or (at your option) any later version.")

class EmptyModel(Qt.QAbstractTableModel):
    def rowCount(self,parent=None):
        return 0

    def columnCount(self,parent=None):
        return 0

    def data(self,index,role=QtCore.Qt.DisplayRole):
        return None

class ODplateItem(object):
    def __init__(self,tc,parentItem):
        self._itemData=tc
        self._parentItem=parentItem
        self._childItems=[]

    def appendChild(self,childItem):
        self._childItems.append(childItem)

    def child(self,row):
        if row >= len(self._childItems):
            return None
        return self._childItems[row]

    def childCount(self):
        return len(self._childItems)

    def row(self):
        """ this item's row in the parent """
        if self._parentItem is not None:
            return self._parentItem._childItems.index(self)

        return 0

    def columnCount(self):
        return 3

    def flags(self):
        if self._parentItem is None or self._parentItem._itemData is None or self._parentItem._parentItem is None:
            #print self._itemData.sampleid, self._itemData.condition, "not checkable"
            return QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEnabled

        #print self._itemData.sampleid, self._itemData.condition, "checkable"
        return QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsUserCheckable

    def _relativeParentIndex(self):
        """
        returns the index of this replicate in its parental replicate group

        FIXME this belongs to Plate.py
        """
        if self._parentItem is None:
            return None
        # find out this child items relative index in parent
        ownIndex=self._itemData.childWellIndices()
        if len(ownIndex) != 1:
            return None
        parentsChildIdcs=self._parentItem._itemData.childWellIndices()
        return parentsChildIdcs.index(ownIndex[0])

    def isChecked(self):
        if self._parentItem is None or self._parentItem._itemData is None:
            return None

        relparidx=self._relativeParentIndex()
        activeIdcs=self._parentItem._itemData.activeChildWellIndices()

        if relparidx is not None and activeIdcs is not None:
            return activeIdcs.count(relparidx) == 1

        return None

    def setChecked(self,value):
        if self._parentItem is None or self._parentItem._itemData is None:
            return None

        relparidx=self._relativeParentIndex()
        self._parentItem._itemData.activateChildWellIndex(relparidx,value)

    def data(self,column):
        if column == 2:
            return self._itemData.activeChildWellIdStr()
        if column == 1:
            return self._itemData.condition
        elif column == 0:
            return self._itemData.sampleid

        return None

    def parent(self):
        return self._parentItem

class ODplateModel(Qt.QAbstractItemModel):
    def __init__(self,odplate):
        super(Qt.QAbstractItemModel, self).__init__()
        self.odplate=odplate
        self._rootItem=ODplateItem(None,None)
        self._setup()

    def _setup(self):
        self._wellToIdx={}
        replicateGroupIdx=0
        for tc in self.odplate.replicateGroups:
            tcitem=ODplateItem(tc,self._rootItem)
            self._rootItem.appendChild(tcitem)
            stat=tc.sampleid+" "+tc.condition+":"
            childIdx=0
            for stc in tc.childWells():
                tcitem.appendChild(ODplateItem(stc,tcitem))
                self._wellToIdx[stc]={'replicateGroupIdx': replicateGroupIdx, 'childIdx': childIdx}
                childIdx+=1
            replicateGroupIdx+=1

    def wellIndex(self,well,column=0):
        """
        Return the QModelIndex of the given well.
        """
        if well not in self._wellToIdx:
            print('well not in dict: '+well.fullId())
            return QtCore.QModelIndex()
        idxGroup = self.index(self._wellToIdx[well]['replicateGroupIdx'],0,QtCore.QModelIndex()); 
        return self.index(self._wellToIdx[well]['childIdx'],column,idxGroup); 

    def rowCount(self,parent):
         if parent.column() > 0:
             return 0;

         if not parent.isValid():
             parentItem=self._rootItem;
         else:
             parentItem=parent.internalPointer()

         return parentItem.childCount();

    def columnCount(self,parent=QtCore.QModelIndex()):
        if parent.isValid():
            return parent.internalPointer().columnCount()
        else:
            return self._rootItem.columnCount();

    def index(self,row,column,parent):
        if not self.hasIndex(row,column,parent):
            print("plate.index("+str(row)+","+str(column)+") hasIndex False")
            return QtCore.QModelIndex();

        if not parent.isValid():
            parentItem = self._rootItem;
        else:
            parentItem = parent.internalPointer();

        childItem=parentItem.child(row)
        if childItem is not None:
            return self.createIndex(row, column, childItem);
        else:
            return QtCore.QModelIndex();

    def parent(self,index):
        """returns index to parent"""
        if not index.isValid():
            return QtCore.QModelIndex();

        childItem=index.internalPointer()
        parentItem=childItem.parent()

        if parentItem == self._rootItem:
            return QtCore.QModelIndex();

        return self.createIndex(parentItem.row(), 0, parentItem);

    def flags(self,index):
        if not index.isValid():
            return None

        item=index.internalPointer()

        return item.flags()

    def data(self,index,role=QtCore.Qt.DisplayRole):
        if not index.isValid():
            return None

        item=index.internalPointer()

        if role == QtCore.Qt.CheckStateRole and index.column() == 0:
            ischk=item.isChecked()
            if ischk is None:
                return None
            return int( QtCore.Qt.Checked if item.isChecked() else QtCore.Qt.Unchecked )

        if role != QtCore.Qt.DisplayRole:
            return None

        return item.data(index.column())

    def setData(self,index,value,role=QtCore.Qt.EditRole):
        if index.column() == 0:
            if role != QtCore.Qt.CheckStateRole:
                return False
            item = index.internalPointer()
            item.setChecked(value)
            self.dataChanged.emit(index, index)
            return True

    def headerData(self,section,orientation,role=QtCore.Qt.DisplayRole):
        if orientation == QtCore.Qt.Horizontal:
            if section > 2:
                return None;
            if role != QtCore.Qt.DisplayRole:
                return None;
            if section == 0:
                return "sample"
            elif section == 1:
                return "condition"
            else:
                return "Well Ids"

        return None

class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = matplotlib.figure.Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self,
                                   Qt.QSizePolicy.Expanding,
                                   Qt.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def sizeHint(self):
        w, h = self.get_width_height()
        return Qt.QSize(w, h)

    def minimumSizeHint(self):
        return Qt.QSize(10, 10)

class DeselectableTreeView(QtGui.QTreeView):
    """
    allows to deselect a row (QTreeView does not allow this in
    SingleSelection mode)

    NOTE hardcoded deselecting full row, which is what we want here
    but it is not a generic solution.
    """
    def mousePressEvent(self, event):
        itemIdx=self.indexAt(event.pos())
        deselect=False
        # only if the user didn't click on the 'expand' arrow we may want to deselect this
        if itemIdx.isValid() and event.pos().x() >= self.visualRect(itemIdx).x() - self.visualRect(self.rootIndex()).x():
            # deselect if this was previsously selected
            deselect=self.selectionModel().isSelected(itemIdx)
        QtGui.QTreeView.mousePressEvent(self, event)
        if deselect:
            self.selectionModel().select(itemIdx, QtGui.QItemSelectionModel.Deselect | QtGui.QItemSelectionModel.Rows)

class ElementChooserWidget(QtGui.QWidget):
    """
    Widget presenting two QListWidgets and arrow buttons to move
    elements from one list to the other.
    """
    def __init__(self,deselected,selected=[],labelNotSelected='',labelSelected=''):
        super(ElementChooserWidget, self).__init__()

        if labelNotSelected != '':
            notselectedLabel=QtGui.QLabel(labelNotSelected)
        self.notselected=QtGui.QListWidget()
        addButton=QtGui.QPushButton('>>')
        remButton=QtGui.QPushButton('<<')
        if labelSelected != '':
            selectedLabel=QtGui.QLabel('Selected for export:')
        self.selected=QtGui.QListWidget()
        self.selected.setDragDropMode(QtGui.QAbstractItemView.InternalMove);

        for item in deselected:
            self.notselected.addItem(item)
        for item in selected:
            self.selected.addItem(item)

        addButton.clicked.connect(self.onAdd)
        remButton.clicked.connect(self.onRemove)

        notselectedLayout=QtGui.QVBoxLayout()
        if labelNotSelected != '':
            notselectedLayout.addWidget(notselectedLabel)
        notselectedLayout.addWidget(self.notselected)

        addRemoveButtonLayout=QtGui.QVBoxLayout()
        addRemoveButtonLayout.addWidget(addButton)
        addRemoveButtonLayout.addWidget(remButton)

        selectedLayout=QtGui.QVBoxLayout()
        if labelSelected != '':
            selectedLayout.addWidget(selectedLabel)
        selectedLayout.addWidget(self.selected)

        addRemoveLayout=QtGui.QHBoxLayout()
        addRemoveLayout.addLayout(notselectedLayout)
        addRemoveLayout.addLayout(addRemoveButtonLayout)
        addRemoveLayout.addLayout(selectedLayout)

        self.setLayout(addRemoveLayout)

    def onAdd(self):
        for item in self.notselected.selectedItems():
            self.selected.addItem(self.notselected.takeItem(self.notselected.row(item)))

    def onRemove(self):
        for item in self.selected.selectedItems():
            self.notselected.addItem(self.selected.takeItem(self.selected.row(item)))


class QCheckBoxRichText(QtGui.QWidget):
    def __init__(self,label): 
        super(QCheckBoxRichText, self).__init__()
        self.checkbox = QtGui.QCheckBox('')
        self.checkbox.stateChanged.connect(self.stateChanged.emit)
        self.label = QtGui.QLabel(label)

        self._myLayout = QtGui.QHBoxLayout()
        self._myLayout.setContentsMargins(QtCore.QMargins(0,0,5,0))
        self._myLayout.addWidget(self.checkbox)
        self._myLayout.addWidget(self.label)

        self.setLayout(self._myLayout)

    stateChanged = QtCore.pyqtSignal(int)

    def setCheckState(self,checkstate):
        self.checkbox.setCheckState(checkstate)

    def checkState(self):
        return self.checkbox.checkState()


class OdReplicateWidgets(object):
    """
    Widgets showing different aspects of Plate/Replicate.
    """

    def __init__(self):
        self.tableModel=None
        self.odplateView=DeselectableTreeView()
        self.overviewmplcnvs=MyMplCanvas()
        self.mplcnvs=MyMplCanvas()
        self._overviewWellRectangles=None
        self._lastTooltipWellIdx=None
        self.overviewmplcnvs.mpl_connect('motion_notify_event', self._updateTooltip)
        self.overviewmplcnvs.mpl_connect('button_press_event',self._buttonPressEvent)
        self.mainstack=QtGui.QStackedWidget()
        self.mainstack.addWidget(self.overviewmplcnvs)
        self.mainstack.addWidget(self.mplcnvs)
        self.navigationWidget = QtGui.QWidget()
        self.ntb=NavigationToolbar(self.mplcnvs,self.navigationWidget)
        self.statusbar=None
        self.statusbarWidgets=[]

        self.hdCorrectionLinearEdit=Qt.QLineEdit()
        self.hdCorrectionLinearEdit.setValidator(QtGui.QDoubleValidator())
        self.hdCorrectionLinearEdit.setFixedWidth(50)
        self.hdCorrectionQuadraticEdit=Qt.QLineEdit()
        self.hdCorrectionQuadraticEdit.setValidator(QtGui.QDoubleValidator())
        self.hdCorrectionQuadraticEdit.setFixedWidth(50)
        self.hdCorrectionCubicEdit=Qt.QLineEdit()
        self.hdCorrectionCubicEdit.setValidator(QtGui.QDoubleValidator())
        self.hdCorrectionCubicEdit.setFixedWidth(50)
        self.hdCorrectionGroup=QtGui.QGroupBox('High density correction')
        self.hdCorrectionLayout=QtGui.QHBoxLayout()
        self.hdCorrectionLayout.addWidget(QtGui.QLabel("OD<sub>cor</sub> ="))
        self.hdCorrectionLayout.addWidget(self.hdCorrectionLinearEdit)
        self.hdCorrectionLayout.addWidget(QtGui.QLabel("OD<sub>obs</sub> +"))
        self.hdCorrectionLayout.addWidget(self.hdCorrectionQuadraticEdit)
        self.hdCorrectionLayout.addWidget(QtGui.QLabel("OD<sup>2</sup><sub>obs</sub> +"))
        self.hdCorrectionLayout.addWidget(self.hdCorrectionCubicEdit)
        self.hdCorrectionLayout.addWidget(QtGui.QLabel("OD<sup>3</sup><sub>obs</sub>"))
        self.hdCorrectionGroup.setLayout(self.hdCorrectionLayout)
        self.hdCorrectionGroup.setSizePolicy(Qt.QSizePolicy.Maximum,
                                             Qt.QSizePolicy.Maximum)

        self.slidingWindowSizeEdit=Qt.QLineEdit()
        intval=QtGui.QIntValidator()
        intval.setRange(3,1042)
        self.slidingWindowSizeEdit.setValidator(intval);
        self.slidingWindowSizeEdit.setFixedWidth(40)
        self.slidingWindowSizeGroup=QtGui.QGroupBox('Fit Window')
        self.slidingWindowSizeLayout=QtGui.QHBoxLayout()
        self.slidingWindowSizeLayout.addWidget(QtGui.QLabel("w ="))
        self.slidingWindowSizeLayout.addWidget(self.slidingWindowSizeEdit)
        self.slidingWindowSizeGroup.setLayout(self.slidingWindowSizeLayout)
        self.slidingWindowSizeGroup.setSizePolicy(Qt.QSizePolicy.Maximum,
                                                  Qt.QSizePolicy.Maximum)


        self.lagAtLogOdEqualsEdit=Qt.QLineEdit()
        lagAtLogOdEqualsval=QtGui.QDoubleValidator()
        self.lagAtLogOdEqualsEdit.setValidator(lagAtLogOdEqualsval);
        self.lagAtLogOdEqualsEdit.setFixedWidth(40)
        self.lagAtLogOdEqualsGroup=QtGui.QGroupBox('Lag at')
        self.lagAtLogOdEqualsLayout=QtGui.QHBoxLayout()
        self.lagAtLogOdEqualsLayout.addWidget(QtGui.QLabel("ln(OD) ="))
        self.lagAtLogOdEqualsLayout.addWidget(self.lagAtLogOdEqualsEdit)
        self.lagAtLogOdEqualsGroup.setLayout(self.lagAtLogOdEqualsLayout)
        self.lagAtLogOdEqualsGroup.setSizePolicy(Qt.QSizePolicy.Maximum,
                                                 Qt.QSizePolicy.Maximum)

        self.cutoffEdit=Qt.QLineEdit()
        cutoffval=QtGui.QDoubleValidator()
        self.cutoffEdit.setValidator(cutoffval);
        self.cutoffEdit.setFixedWidth(40)
        self.cutoffGroup=QtGui.QGroupBox('ln cutoff')
        self.cutoffLayout=QtGui.QHBoxLayout()
        self.cutoffLayout.addWidget(QtGui.QLabel("ln(OD) >"))
        self.cutoffLayout.addWidget(self.cutoffEdit)
        self.cutoffGroup.setLayout(self.cutoffLayout)
        self.cutoffGroup.setSizePolicy(Qt.QSizePolicy.Maximum,
                                       Qt.QSizePolicy.Maximum)

        self.maxGrowthLowerTimeCutoffEdit=Qt.QLineEdit()
        maxGrowthLowerCutoffval=QtGui.QDoubleValidator()
        self.maxGrowthLowerTimeCutoffEdit.setValidator(maxGrowthLowerCutoffval);
        self.maxGrowthLowerTimeCutoffEdit.setFixedWidth(40)

        self.maxGrowthUpperTimeCutoffEdit=Qt.QLineEdit()
        maxGrowthUpperCutoffval=QtGui.QDoubleValidator()
        self.maxGrowthUpperTimeCutoffEdit.setValidator(maxGrowthUpperCutoffval);
        self.maxGrowthUpperTimeCutoffEdit.setFixedWidth(40)

        self.allowMaxGrowthrateAtLowerCutoffCheckbox=QtGui.QCheckBox("allow at\nlow. cutoff")

        self.maxGrowthTimeCutoffGroup=QtGui.QGroupBox('Max growth cutoff')
        self.maxGrowthTimeCutoffLayout=QtGui.QHBoxLayout()
        self.maxGrowthTimeCutoffLayout.addWidget(self.maxGrowthLowerTimeCutoffEdit)
        self.maxGrowthTimeCutoffLayout.addWidget(QtGui.QLabel("< t <"))
        self.maxGrowthTimeCutoffLayout.addWidget(self.maxGrowthUpperTimeCutoffEdit)
        self.maxGrowthTimeCutoffLayout.addWidget(self.allowMaxGrowthrateAtLowerCutoffCheckbox)
        self.maxGrowthTimeCutoffGroup.setLayout(self.maxGrowthTimeCutoffLayout)
        self.maxGrowthTimeCutoffGroup.setSizePolicy(Qt.QSizePolicy.Maximum,
                                                    Qt.QSizePolicy.Maximum)

        self.allowGrowthyieldSlopeNStderrAwayFromZeroSpinbox=QtGui.QSpinBox()
        self.allowGrowthyieldSlopeNStderrAwayFromZeroSpinbox.setFixedWidth(50)

        self.allowGrowthyieldSlopeNStderrAwayFromZeroGroup=QtGui.QGroupBox('Growthyield')
        self.allowGrowthyieldSlopeNStderrAwayFromZeroLayout=QtGui.QHBoxLayout()
        itslabels=QtGui.QLabel("allow\nn stdErr")
        itslabels.setToolTip('Allow growth yield\'s slope to be zero within n standard errors.')
        self.allowGrowthyieldSlopeNStderrAwayFromZeroLayout.addWidget(itslabels)
        self.allowGrowthyieldSlopeNStderrAwayFromZeroLayout.addWidget(self.allowGrowthyieldSlopeNStderrAwayFromZeroSpinbox)
        self.allowGrowthyieldSlopeNStderrAwayFromZeroGroup.setLayout(self.allowGrowthyieldSlopeNStderrAwayFromZeroLayout)
        self.allowGrowthyieldSlopeNStderrAwayFromZeroGroup.setSizePolicy(Qt.QSizePolicy.Maximum,
                                                                         Qt.QSizePolicy.Maximum)

        sval=QtGui.QDoubleValidator()
        sval.setRange(0.0,314.0,9)
        self.smoothingSEdit=Qt.QLineEdit()
        self.smoothingSEdit.setValidator(sval)
        self.smoothingSEdit.setFixedWidth(40)
        self.smoothingSEdit.setToolTip('Positive smoothing factor used to choose the number of knots.')
        kval=QtGui.QIntValidator()
        kval.setRange(1,5)
        self.smoothingKEdit=Qt.QLineEdit()
        self.smoothingKEdit.setValidator(kval)
        self.smoothingKEdit.setFixedWidth(40)
        self.smoothingKEdit.setToolTip('Degree of the smoothing spline. Must be <= 5.')
        self.smoothingGroup=QtGui.QGroupBox('Smoothing')
        self.smoothingLayout=QtGui.QHBoxLayout()
        self.smoothingLayout.addWidget(QtGui.QLabel("k ="))
        self.smoothingLayout.addWidget(self.smoothingKEdit)
        self.smoothingLayout.addWidget(QtGui.QLabel(" s ="))
        self.smoothingLayout.addWidget(self.smoothingSEdit)
        self.smoothingGroup.setLayout(self.smoothingLayout)
        self.smoothingGroup.setSizePolicy(Qt.QSizePolicy.Maximum,
                                          Qt.QSizePolicy.Maximum)


        self.showRawCheck=QtGui.QCheckBox("raw")
        self.showRawCheck.setCheckState(QtCore.Qt.Unchecked)
        self.showRawCheck.setToolTip("Show raw optical density (high density correction not applied, background not subtracted).")
        self.showBackgroundCheck=QtGui.QCheckBox("background")
        self.showBackgroundCheck.setCheckState(QtCore.Qt.Unchecked)
        self.showBackgroundCheck.setToolTip("Show raw optical density of background wells.")
        self.showSingleCheck=QtGui.QCheckBox("individual")
        self.showSingleCheck.setCheckState(QtCore.Qt.Unchecked)
        self.showSingleCheck.setToolTip("Show optical density (OD) of individual wells that make up the averaged OD.")
        self.showSmoothedCheck=QtGui.QCheckBox("smoothed")
        self.showSmoothedCheck.setCheckState(QtCore.Qt.Checked)
        self.showSmoothedCheck.setToolTip("Show optical density processed by a smoothing spline.")

        self.showLogCheck=QtGui.QCheckBox("ln(OD)")
        self.showLogCheck.setCheckState(QtCore.Qt.Checked)
        self.showLogCheck.setToolTip("Show natural logarithm of optical density.")
        self.showLogSmoothedCheck=QtGui.QCheckBox("ln(OD) smoothed")
        self.showLogSmoothedCheck.setCheckState(QtCore.Qt.Unchecked)
        self.showLogSmoothedCheck.setToolTip("Show natural logarithm of optical density processed by a smoothing spline.")

        self.showExpFitsOd0MuCheck=QCheckBoxRichText('&mu; from &mu;,OD<sub>0</sub> fit')
        self.showExpFitsOd0MuCheck.setCheckState(QtCore.Qt.Checked)
        self.showExpFitsOd0MuCheck.setToolTip("Show growth rate derived from fits (in window around data point) of exponential functions to the data.")
        self.showMaxGrowthrateCheck=QCheckBoxRichText('&mu;<sub>max</sub> from &mu;,OD<sub>0</sub> fit')
        self.showMaxGrowthrateCheck.setCheckState(QtCore.Qt.Checked)
        self.showMaxGrowthrateCheck.setToolTip("Show maximal growth rate derived from fits (in window around data point) of exponential functions to the data.")

        self.showGrowthyieldCheck=QtGui.QCheckBox("yield")
        self.showGrowthyieldCheck.setCheckState(QtCore.Qt.Checked)
        self.showGrowthyieldCheck.setToolTip("Show growth yield.")

        self.showLogDerCheck=QtGui.QCheckBox("d[ln(OD)]/dt")
        self.showLogDerCheck.setCheckState(QtCore.Qt.Unchecked)
        self.showLogDerCheck.setToolTip("Show derivative of log(OD).")
        self.showLogDerNonLogCheck=QtGui.QCheckBox("1/OD d[OD]/dt")
        self.showLogDerNonLogCheck.setCheckState(QtCore.Qt.Unchecked)
        self.showLogDerNonLogCheck.setToolTip("Show growth rate (derivative) derived by chain rule of logarithm of optical density.")
        self.showLogDerNonLogSmoothedCheck=QtGui.QCheckBox("1/OD d[OD]/dt smoothed")
        self.showLogDerNonLogSmoothedCheck.setCheckState(QtCore.Qt.Unchecked)
        self.showLogDerNonLogSmoothedCheck.setToolTip("Show growth rate (derivative) derived by chain rule of smoothed logarithm of optical density.")
        self.showMaxGrowthrateFromLogOdCheck=QCheckBoxRichText(Qt.QString('&mu;<sub>max</sub> from 1/OD d[OD]/dt smoothed'))
        self.showMaxGrowthrateFromLogOdCheck.setToolTip("Show maximal growth (derivative) derived by chain rule of smoothed logarithm of optical density.")
        self.showMaxGrowthrateFromLogOdCheck.setCheckState(QtCore.Qt.Unchecked)

        self.showAlternativeGrowthrateExtraction=True

        self.showAllowMaxGrowthrateAtLowerCutoff=False
        self.hasAllowMaxGrowthrateAtLowerCutoff=False
        self.showAllowGrowthyieldSlopeNStderrAwayFromZero=False
        self.hasAllowGrowthyieldSlopeNStderrAwayFromZero=False

        self.selectedReplicate=None

        self.clear()

    def setAlternativeGrowthrateExtractionVisible(self,visible):
        self.showAlternativeGrowthrateExtraction=visible
        if not visible:
            # make sure check state of these is unchecked
            self.showLogDerNonLogCheck.setCheckState(QtCore.Qt.Unchecked)
            self.showLogDerNonLogSmoothedCheck.setCheckState(QtCore.Qt.Unchecked)
            self.showMaxGrowthrateFromLogOdCheck.setCheckState(QtCore.Qt.Unchecked)
        self.showLogDerNonLogCheck.setVisible(visible)
        self.showLogDerNonLogSmoothedCheck.setVisible(visible)
        self.showMaxGrowthrateFromLogOdCheck.setVisible(visible)

    def clear(self):
        self.slidingWindowSizeEdit.setEnabled(False)
        self.cutoffEdit.setEnabled(False)
        self.lagAtLogOdEqualsEdit.setEnabled(False)
        self.maxGrowthLowerTimeCutoffEdit.setEnabled(False)
        self.maxGrowthUpperTimeCutoffEdit.setEnabled(False)
        self.allowMaxGrowthrateAtLowerCutoffCheckbox.setEnabled(False)
        self.hdCorrectionLinearEdit.setEnabled(False)
        self.hdCorrectionQuadraticEdit.setEnabled(False)
        self.hdCorrectionCubicEdit.setEnabled(False)
        self.smoothingKEdit.setEnabled(False)
        self.smoothingSEdit.setEnabled(False)

        self.slidingWindowSizeEdit.clear()
        self.cutoffEdit.clear()
        self.lagAtLogOdEqualsEdit.clear()
        self.maxGrowthLowerTimeCutoffEdit.clear()
        self.maxGrowthUpperTimeCutoffEdit.clear()
        self.allowMaxGrowthrateAtLowerCutoffCheckbox.blockSignals(True)
        self.allowMaxGrowthrateAtLowerCutoffCheckbox.setCheckState(QtCore.Qt.Unchecked)
        self.allowMaxGrowthrateAtLowerCutoffCheckbox.setTristate(False)
        self.allowMaxGrowthrateAtLowerCutoffCheckbox.blockSignals(False)
        self.allowMaxGrowthrateAtLowerCutoffCheckbox.setVisible(False)

        self.allowGrowthyieldSlopeNStderrAwayFromZeroSpinbox.blockSignals(True)
        self.allowGrowthyieldSlopeNStderrAwayFromZeroSpinbox.setValue(1)
        self.allowGrowthyieldSlopeNStderrAwayFromZeroSpinbox.blockSignals(False)
        self.allowGrowthyieldSlopeNStderrAwayFromZeroGroup.setVisible(False)

        self.hdCorrectionLinearEdit.clear()
        self.hdCorrectionQuadraticEdit.clear()
        self.hdCorrectionCubicEdit.clear()
        self.smoothingKEdit.clear()
        self.smoothingSEdit.clear()

        self.showAllowMaxGrowthrateAtLowerCutoff=False
        self.hasAllowMaxGrowthrateAtLowerCutoff=False
        self.showAllowGrowthyieldSlopeNStderrAwayFromZero=False
        self.hasAllowGrowthyieldSlopeNStderrAwayFromZero=False

        self.mplcnvs.fig.clear()
        self.mplcnvs.draw()
        self.overviewmplcnvs.fig.clear()
        self.overviewmplcnvs.draw()

        self.statusToStatusBar(None)

        self.odplateView.setModel(EmptyModel())
        self.tableModel=None
        self.odplate=None

    def disconnect(self):
        """
        NOTE this cannot be inside clear as disconnect fails if nothing is bound
        """
        self.slidingWindowSizeEdit.textEdited.disconnect()
        self.cutoffEdit.textEdited.disconnect()
        self.lagAtLogOdEqualsEdit.textEdited.disconnect()
        self.maxGrowthLowerTimeCutoffEdit.textEdited.disconnect()
        self.maxGrowthUpperTimeCutoffEdit.textEdited.disconnect()
        self.allowMaxGrowthrateAtLowerCutoffCheckbox.stateChanged.disconnect()
        self.allowGrowthyieldSlopeNStderrAwayFromZeroSpinbox.valueChanged.disconnect()
        self.hdCorrectionLinearEdit.textEdited.disconnect()
        self.hdCorrectionQuadraticEdit.textEdited.disconnect()
        self.hdCorrectionCubicEdit.textEdited.disconnect()
        self.smoothingKEdit.textEdited.disconnect()
        self.smoothingSEdit.textEdited.disconnect()

        self.showRawCheck.stateChanged.disconnect()
        self.showSingleCheck.stateChanged.disconnect()
        self.showBackgroundCheck.stateChanged.disconnect()
        self.showSmoothedCheck.stateChanged.disconnect()
        self.showLogCheck.stateChanged.disconnect()
        self.showLogSmoothedCheck.stateChanged.disconnect()
        self.showMaxGrowthrateCheck.stateChanged.disconnect()
        self.showMaxGrowthrateFromLogOdCheck.stateChanged.disconnect()
        self.showLogDerCheck.stateChanged.disconnect()
        self.showLogDerNonLogCheck.stateChanged.disconnect()
        self.showLogDerNonLogSmoothedCheck.stateChanged.disconnect()
        self.showExpFitsOd0MuCheck.stateChanged.disconnect()
        self.showGrowthyieldCheck.stateChanged.disconnect()

        self.tableModel.dataChanged.disconnect()

    def setPlate(self,odplate):
        self.odplate=odplate
        self.tableModel=ODplateModel(self.odplate)
        self.slidingWindowSizeEdit.setEnabled(True)
        self.cutoffEdit.setEnabled(True)
        self.lagAtLogOdEqualsEdit.setEnabled(True)
        self.maxGrowthLowerTimeCutoffEdit.setEnabled(True)
        self.maxGrowthUpperTimeCutoffEdit.setEnabled(True)
        self.allowMaxGrowthrateAtLowerCutoffCheckbox.setEnabled(True)
        self.hdCorrectionLinearEdit.setEnabled(True)
        self.hdCorrectionQuadraticEdit.setEnabled(True)
        self.hdCorrectionCubicEdit.setEnabled(True)
        self.smoothingKEdit.setEnabled(True)
        self.smoothingSEdit.setEnabled(True)

        self.slidingWindowSizeEdit.textEdited.connect(self.slidingWindowSizeEdited)
        self.cutoffEdit.textEdited.connect(self.cutoffEdited)
        self.lagAtLogOdEqualsEdit.textEdited.connect(self.lagAtEdited)
        self.maxGrowthLowerTimeCutoffEdit.textEdited.connect(self.maxGrowthLowerTimeCutoffEdited)
        self.maxGrowthUpperTimeCutoffEdit.textEdited.connect(self.maxGrowthUpperTimeCutoffEdited)
        self.allowMaxGrowthrateAtLowerCutoffCheckbox.stateChanged.connect(self.allowMaxGrowthrateAtLowerCutoffEdited)
        self.allowGrowthyieldSlopeNStderrAwayFromZeroSpinbox.valueChanged.connect(self.allowGrowthyieldSlopeNStderrAwayFromZeroEdited)
        self.hdCorrectionLinearEdit.textEdited.connect(self.hdCorrectionLinearEdited)
        self.hdCorrectionQuadraticEdit.textEdited.connect(self.hdCorrectionQuadraticEdited)
        self.hdCorrectionCubicEdit.textEdited.connect(self.hdCorrectionCubicEdited)
        self.smoothingSEdit.textEdited.connect(self.smoothingSEdited)
        self.smoothingKEdit.textEdited.connect(self.smoothingKEdited)

        self.showRawCheck.stateChanged.connect(self.updatePlateFig)
        self.showSingleCheck.stateChanged.connect(self.updatePlateFig)
        self.showBackgroundCheck.stateChanged.connect(self.updatePlateFig)
        self.showSmoothedCheck.stateChanged.connect(self.updatePlateFig)
        self.showLogCheck.stateChanged.connect(self.updatePlateFig)
        self.showLogSmoothedCheck.stateChanged.connect(self.updatePlateFig)
        self.showMaxGrowthrateCheck.stateChanged.connect(self.updatePlateFig)
        self.showMaxGrowthrateFromLogOdCheck.stateChanged.connect(self.updatePlateFig)
        self.showLogDerCheck.stateChanged.connect(self.updatePlateFig)
        self.showLogDerNonLogCheck.stateChanged.connect(self.updatePlateFig)
        self.showLogDerNonLogSmoothedCheck.stateChanged.connect(self.updatePlateFig)
        self.showExpFitsOd0MuCheck.stateChanged.connect(self.updatePlateFig)
        self.showGrowthyieldCheck.stateChanged.connect(self.updatePlateFig)

        self.hasAllowMaxGrowthrateAtLowerCutoff=self.plateHasAllowMaxGrowthrateAtLowerCutoff()
        if self.showAllowMaxGrowthrateAtLowerCutoff or self.hasAllowMaxGrowthrateAtLowerCutoff:
            self.allowMaxGrowthrateAtLowerCutoffCheckbox.setVisible(True)
        # NOTE if the plate default for allowMaxGrowthrateAtLowerCutoff is True
        # (which should only happen for old or manipulated datasets),
        # make this tristate
        if self.odplate.getParameter('allowMaxGrowthrateAtLowerCutoff'):
            self.allowMaxGrowthrateAtLowerCutoffCheckbox.setTristate(True)

        self.hasAllowGrowthyieldSlopeNStderrAwayFromZero=self.plateHasAllowGrowthyieldSlopeNStderrAwayFromZero()
        if self.hasAllowGrowthyieldSlopeNStderrAwayFromZero:
            self.allowGrowthyieldSlopeNStderrAwayFromZeroGroup.setVisible(True)

        self.tableModel.dataChanged.connect(self.updatePlateFig)

        self.odplateView.setModel(self.tableModel)
        self.odplateView.header().setResizeMode(0, Qt.QHeaderView.ResizeToContents)
        self.odplateView.header().setResizeMode(1, Qt.QHeaderView.ResizeToContents)

        self.odplateView.selectionModel().selectionChanged.connect(self.selectionChanged)
        self.odplateView.setFocus()
        firstIndex=self.tableModel.index(0,0,self.odplateView.currentIndex())
        self._setSample(self.odplate)

    def plateHasAllowMaxGrowthrateAtLowerCutoff(self):
        """
        Sets self.hasAllowMaxGrowthrateAtLowerCutoff to True if at
        least one replicate has the parameter
        allowMaxGrowthrateAtLowerCutoff set to True.
        """
        self.hasAllowMaxGrowthrateAtLowerCutoff=False

        if self.odplate is None:
            return

        if self.odplate.getParameter('allowMaxGrowthrateAtLowerCutoff') is not None and self.odplate.getParameter('allowMaxGrowthrateAtLowerCutoff'):
            self.hasAllowMaxGrowthrateAtLowerCutoff=True
            return True
        for tc in self.odplate.wells:
            if tc.getParameter('allowMaxGrowthrateAtLowerCutoff') is not None and tc.getParameter('allowMaxGrowthrateAtLowerCutoff'):
                self.hasAllowMaxGrowthrateAtLowerCutoff=True
            if self.hasAllowMaxGrowthrateAtLowerCutoff:
                break
        for tc in self.odplate.replicateGroups:
            if tc.getParameter('allowMaxGrowthrateAtLowerCutoff') is not None and tc.getParameter('allowMaxGrowthrateAtLowerCutoff'):
                self.hasAllowMaxGrowthrateAtLowerCutoff=True
            if self.hasAllowMaxGrowthrateAtLowerCutoff:
                break

        return self.hasAllowMaxGrowthrateAtLowerCutoff

    def plateHasAllowGrowthyieldSlopeNStderrAwayFromZero(self):
        """
        Sets self.hasAllowGrowthyieldSlopeNStderrAwayFromZero to True
        if at least one replicate has the parameter
        allowGrowthyieldSlopeNStderrAwayFromZero set to a value
        different from one.
        """
        self.hasAllowGrowthyieldSlopeNStderrAwayFromZero=False

        if self.odplate is None:
            return

        if (self.odplate.getParameter('allowGrowthyieldSlopeNStderrAwayFromZero') is not None
            and self.odplate.getParameter('allowGrowthyieldSlopeNStderrAwayFromZero') != 1):
            self.hasAllowGrowthyieldSlopeNStderrAwayFromZero=True
            return True
        for tc in self.odplate.wells:
            if (tc.getParameter('allowGrowthyieldSlopeNStderrAwayFromZero') is not None
                and tc.getParameter('allowGrowthyieldSlopeNStderrAwayFromZero') != 1):
                self.hasAllowGrowthyieldSlopeNStderrAwayFromZero=True
            if self.hasAllowGrowthyieldSlopeNStderrAwayFromZero:
                break

        return self.hasAllowGrowthyieldSlopeNStderrAwayFromZero

    def selectionChanged(self,selected,deselected):
        """
        the model view selection changed
        """
        if len(selected.indexes()):
            selindex=selected.indexes()[0]
            tc=selindex.internalPointer()._itemData
            self._setSample(tc)
        else:
            self._setSample(self.odplate)

    def updateParameterLineEditColor(self,tc,lineedit,par):
        parval=tc.getParameter(par)
        if parval is not None:
            if not tc.parameterIsExplicitlySet(par):
                if tc.activeChildReplicatesHaveExplicitParameter(par):
                    lineedit.setStyleSheet("background-color: qlineargradient(x1: 0, y1: 0.25, x2: 0, y2: 0.75, stop: 0 lightgrey, stop: 1 white)")
                else:
                    ## NOTE setting the background of a QSpinBox to plain lightgrey results in a weird frame (at least on Debian wheezy)
                    #lineedit.setStyleSheet("background-color: lightgrey;")
                    lineedit.setStyleSheet("background-color: qlineargradient(x1: 0, y1: 0.25, x2: 0, y2: 0.75, stop: 0 lightgrey, stop: 1 lightgrey)")
            else:
                if tc.activeChildReplicatesHaveExplicitParameter(par):
                    lineedit.setStyleSheet("background-color: qlineargradient(x1: 0.25, y1: 0, x2: 0.75, y2: 0, stop: 0 lightgrey, stop: 1 white)")
                else:
                    ## NOTE setting the background of a QSpinBox to plain white results in a weird frame (at least on Debian wheezy)
                    #lineedit.setStyleSheet("background-color: white")
                    lineedit.setStyleSheet("background-color: qlineargradient(x1: 0.25, y1: 0, x2: 0.75, y2: 0, stop: 0 white, stop: 1 white)")
        else:
            if tc.activeChildReplicatesHaveExplicitParameter(par):
                lineedit.setStyleSheet("background-color: qlineargradient(x1: 0, y1: 0.25, x2: 0, y2: 0.75, stop: 0 lightgrey, stop: 1 white)")
            else:
                ## NOTE setting the background of a QSpinBox to plain lightgrey results in a weird frame (at least on Debian wheezy)
                #lineedit.setStyleSheet("background-color: lightgrey")
                lineedit.setStyleSheet("background-color: qlineargradient(x1: 0, y1: 0.25, x2: 0, y2: 0.75, stop: 0 lightgrey, stop: 1 lightgrey)")

    def updateParameterLineEdit(self,tc,lineedit,par):
        parval=tc.getParameter(par)

        #lineedit.setReadOnly(not tc.parameterIsEditible(par))
        lineedit.setEnabled(tc.parameterIsEditible(par))
        if parval is not None:
            lineedit.setText(str(parval))
            lineedit.setCursorPosition(0)
        else:
            lineedit.clear()
        self.updateParameterLineEditColor(tc,lineedit,par)

    def updateParameterCheckboxTristate(self,tc,checkbox,par):
        parval=tc.getParameter(par)

        checkbox.setEnabled(tc.parameterIsEditible(par))
        checkbox.blockSignals(True)
        if parval is None:
            checkbox.setCheckState(QtCore.Qt.PartiallyChecked)
        elif parval is True:
            checkbox.setCheckState(QtCore.Qt.Checked)
        elif parval is False:
            checkbox.setCheckState(QtCore.Qt.Unchecked)
        self.updateParameterLineEditColor(tc,checkbox,par)
        checkbox.blockSignals(False)

    def updateParameterCheckboxUncheckedIsNone(self,tc,checkbox,par):
        parval=tc.getParameter(par)

        checkbox.setEnabled(tc.parameterIsEditible(par))
        checkbox.blockSignals(True)
        if parval:
            checkbox.setCheckState(QtCore.Qt.Checked)
        else:
            checkbox.setCheckState(QtCore.Qt.Unchecked)
        self.updateParameterLineEditColor(tc,checkbox,par)
        checkbox.blockSignals(False)

    def updateParameterSpinbox(self,tc,spinbox,par):
        parval=tc.getParameter(par)

        spinbox.setEnabled(tc.parameterIsEditible(par))
        spinbox.blockSignals(True)
        if parval is not None:
            spinbox.setValue(parval)
        else:
            spinbox.clear()
        self.updateParameterLineEditColor(tc,spinbox,par)
        spinbox.blockSignals(False)

    def _setSample(self,tc):
        self.selectedReplicate=tc

        self.updateParameterLineEdit(tc,self.hdCorrectionLinearEdit,"hdCorrectionLinear")
        self.updateParameterLineEdit(tc,self.hdCorrectionQuadraticEdit,"hdCorrectionQuadratic")
        self.updateParameterLineEdit(tc,self.hdCorrectionCubicEdit,"hdCorrectionCubic")
        self.updateParameterLineEdit(tc,self.smoothingKEdit,"smoothingK")
        self.updateParameterLineEdit(tc,self.smoothingSEdit,"smoothingS")
        self.updateParameterLineEdit(tc,self.slidingWindowSizeEdit,"slidingWindowSize")
        self.updateParameterLineEdit(tc,self.cutoffEdit,"logOdCutoff")
        self.updateParameterLineEdit(tc,self.lagAtLogOdEqualsEdit,"lagAtLogOdEquals")
        self.updateParameterLineEdit(tc,self.maxGrowthLowerTimeCutoffEdit,"maxGrowthLowerTimeCutoff")
        self.updateParameterLineEdit(tc,self.maxGrowthUpperTimeCutoffEdit,"maxGrowthUpperTimeCutoff")
        if self.allowMaxGrowthrateAtLowerCutoffCheckbox.isTristate():
            self.updateParameterCheckboxTristate(tc,self.allowMaxGrowthrateAtLowerCutoffCheckbox,"allowMaxGrowthrateAtLowerCutoff")
        else:
            self.updateParameterCheckboxUncheckedIsNone(tc,self.allowMaxGrowthrateAtLowerCutoffCheckbox,"allowMaxGrowthrateAtLowerCutoff")
        self.updateParameterSpinbox(tc,self.allowGrowthyieldSlopeNStderrAwayFromZeroSpinbox,"allowGrowthyieldSlopeNStderrAwayFromZero")

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
        self._overviewWellRectangles=None
        if self._lastTooltipWellIdx is not None:
            QtGui.QToolTip.hideText()
            self._lastTooltipWellIdx=None
        self.statusToStatusBar(None)
        if self.selectedReplicate is None:
            return
        QtGui.QApplication.setOverrideCursor(Qt.QCursor(QtCore.Qt.WaitCursor))

        try:
            if isinstance(self.selectedReplicate, Plate):
                # this is the whole plate, so show an overview of all wells
                self.overviewmplcnvs.fig.clear()
                ax=self.overviewmplcnvs.fig.add_subplot(1,1,1,axisbg='#eeeeee')
                self._overviewWellRectangles=odPlateOverviewToAxes(ax,self.odplate)
                if self.odplate._loadStatus is not None:
                    status=self.odplate._loadStatus
                else:
                    status=StatusMessage()
                self.overviewmplcnvs.fig.tight_layout()
                self.overviewmplcnvs.draw()
                self.mainstack.setCurrentWidget(self.overviewmplcnvs)
                # coordinates are wrong, so don't show them;
                self.ntb.coordinates=False
                self.ntb.setEnabled(False)

            else:
                self.ntb.setEnabled(True)
                self.ntb.coordinates=True
                self.mplcnvs.fig.clear()
                tc=self.selectedReplicate
                addToTitle=tc.activeChildWellIdStr()
                status=replicateToFig(self.mplcnvs.fig,tc,
                                      showTitle=False,
                                      showRaw=(self.showRawCheck.checkState()  == QtCore.Qt.Checked),
                                      showSingle=(self.showSingleCheck.checkState()  == QtCore.Qt.Checked),
                                      showBackground=(self.showBackgroundCheck.checkState() == QtCore.Qt.Checked),
                                      showSmoothed=(self.showSmoothedCheck.checkState()  == QtCore.Qt.Checked),
                                      showSmoothedDerivativeLinear=(self.showSmoothedCheck.checkState() == QtCore.Qt.Checked),
                                      showLogod=(self.showLogCheck.checkState()  == QtCore.Qt.Checked),
                                      showLogodSmoothed=(self.showLogSmoothedCheck.checkState()  == QtCore.Qt.Checked),
                                      showMaxGrowthrate=(self.showMaxGrowthrateCheck.checkState()  == QtCore.Qt.Checked),
                                      showMaxGrowthrateFromLogOdDerivative=(self.showMaxGrowthrateFromLogOdCheck.checkState()  == QtCore.Qt.Checked),
                                      showLogOdDerivative=(self.showLogDerCheck.checkState()  == QtCore.Qt.Checked),
                                      showLogOdDerivativeFromNonLog=(self.showLogDerNonLogCheck.checkState()  == QtCore.Qt.Checked),
                                      showLogOdDerivativeFromNonLogSmoothed=(self.showLogDerNonLogSmoothedCheck.checkState()  == QtCore.Qt.Checked),
                                      showExpFitsOd0Mu=(self.showExpFitsOd0MuCheck.checkState()  == QtCore.Qt.Checked),
                                      showGrowthyield=(self.showGrowthyieldCheck.checkState()  == QtCore.Qt.Checked),
                                      addToTitle=addToTitle)
                self.mplcnvs.fig.tight_layout()
                self.mplcnvs.draw()
                self.mainstack.setCurrentWidget(self.mplcnvs)

            if status is not None:
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
            self.overviewmplcnvs.fig.clear()
            self.overviewmplcnvs.draw()
            self.statusToStatusBar(None)
            QtGui.QApplication.restoreOverrideCursor()
            QtGui.QMessageBox.warning(self.mainstack,'Unknown error',msg)
            return

    def _eventToOverviewRectangle(self,event):
        # FIXME all this could probably be done more efficiently
        if self._overviewWellRectangles is not None and event.inaxes:
            x, y = event.xdata, event.ydata
            # check whether mouse pointer position is within a reactangle
            for rec in self._overviewWellRectangles:
                if rec['x1'] < x and x < rec['x2'] and rec['y1'] < y and y < rec['y2']:
                    return rec
        return None

    def _updateTooltip(self, event):
        """
        Slot to react to mouse movements and show a tooltip if wanted.
        """
        rec = self._eventToOverviewRectangle(event)
        if rec is not None:
            if self._lastTooltipWellIdx is None or self._lastTooltipWellIdx != rec['wellidx']:
                QtGui.QToolTip.showText(QtGui.QCursor.pos(), rec['tooltip'])
                self._lastTooltipWellIdx = rec['wellidx']
            return
        if self._lastTooltipWellIdx is not None:
            QtGui.QToolTip.hideText()
            self._lastTooltipWellIdx = None

    def _buttonPressEvent(self, event):
        """
        Slot to react to mouse clicks, select the well under the mouse cursor if in overview mode.
        """
        rec = self._eventToOverviewRectangle(event)
        if rec is not None:
            wellsel=Qt.QItemSelection(self.tableModel.wellIndex(rec['well'],0),
                                      self.tableModel.wellIndex(rec['well'],2))
            self.odplateView.scrollTo(self.tableModel.wellIndex(rec['well'],0))
            self.odplateView.selectionModel().select(wellsel, Qt.QItemSelectionModel.Select)

    def hdCorrectionLinearEdited(self,newval):
        # we neglect newval as we get all the values from the QLineEdits
        lval, ok=self.hdCorrectionLinearEdit.text().toDouble()
        if ok:
            self.odplate.setHighDensityCorrectionLinear(lval)
        else:
            self.odplate.setHighDensityCorrectionLinear(None)
        self.updatePlateFig()

    def hdCorrectionQuadraticEdited(self,newval):
        # we neglect newval as we get all the values from the QLineEdits
        qval, ok=self.hdCorrectionQuadraticEdit.text().toDouble()
        if ok:
            self.odplate.setHighDensityCorrectionQuadratic(qval)
        else:
            self.odplate.setHighDensityCorrectionQuadratic(None)
        self.updatePlateFig()

    def hdCorrectionCubicEdited(self,newval):
        # we neglect newval as we get all the values from the QLineEdits
        cval, ok=self.hdCorrectionCubicEdit.text().toDouble()
        if ok:
            self.odplate.setHighDensityCorrectionCubic(cval)
        else:
            self.odplate.setHighDensityCorrectionCubic(None)
        self.updatePlateFig()

    def smoothingSEdited(self,newval):
        # we neglect newval as we get all the values from the QLineEdits
        sval, ok=self.smoothingSEdit.text().toDouble()
        if ok and sval >= 0.:
            self.odplate.setSmoothingS(sval)
        else:
            self.odplate.setSmoothingS(None)
        self.updatePlateFig()

    def smoothingKEdited(self,newval):
        # we neglect newval as we get all the values from the QLineEdits
        kval, ok=self.smoothingKEdit.text().toInt()
        if ok and kval >= 0:
            self.odplate.setSmoothingK(kval)
        else:
            self.odplate.setSmoothingK(None)
        self.updatePlateFig()

    def slidingWindowSizeEdited(self,newval):
        val, ok=newval.toInt()
        if ok and val > 2:
            self.odplate.setSlidingWindowSize(val)
        else:
            self.odplate.setSlidingWindowSize(None)
        self.updatePlateFig()

    def cutoffEdited(self,newval):
        if self.selectedReplicate is None:
            return
        tc=self.selectedReplicate

        val, ok=newval.toDouble()
        if ok:
            tc.setLogOdCutoff(val)
        else:
            tc.setLogOdCutoff(None)
        self.updateParameterLineEditColor(tc,self.cutoffEdit,"logOdCutoff")
        self.updatePlateFig()

    def lagAtEdited(self,newval):
        if self.selectedReplicate is None:
            return
        tc=self.selectedReplicate

        val, ok=newval.toDouble()
        if ok:
            tc.setLagAtLogOdEquals(val)
        else:
            tc.setLagAtLogOdEquals(None)
        self.updateParameterLineEditColor(tc,self.lagAtLogOdEqualsEdit,"lagAtLogOdEquals")
        self.updatePlateFig()

    def maxGrowthLowerTimeCutoffEdited(self,newval):
        if self.selectedReplicate is None:
            return

        tc=self.selectedReplicate

        val, ok=newval.toDouble()
        if ok:
            tc.setMaxGrowthLowerTimeCutoff(val)
        else:
            tc.setMaxGrowthLowerTimeCutoff(None)
        self.updateParameterLineEditColor(tc,self.maxGrowthLowerTimeCutoffEdit,"maxGrowthLowerTimeCutoff")
        self.updatePlateFig()

    def maxGrowthUpperTimeCutoffEdited(self,newval):
        if self.selectedReplicate is None:
            return

        tc=self.selectedReplicate

        val, ok=newval.toDouble()
        if ok:
            tc.setMaxGrowthUpperTimeCutoff(val)
        else:
            tc.setMaxGrowthUpperTimeCutoff(None)
        self.updateParameterLineEditColor(tc,self.maxGrowthUpperTimeCutoffEdit,"maxGrowthUpperTimeCutoff")
        self.updatePlateFig()

    def allowMaxGrowthrateAtLowerCutoffEdited(self):
        if self.selectedReplicate is None:
            return

        tc=self.selectedReplicate
        if self.allowMaxGrowthrateAtLowerCutoffCheckbox.checkState() == QtCore.Qt.Checked:
            tc.setAllowMaxGrowthrateAtLowerCutoff(True)
        elif self.allowMaxGrowthrateAtLowerCutoffCheckbox.checkState() == QtCore.Qt.Unchecked:
            if self.allowMaxGrowthrateAtLowerCutoffCheckbox.isTristate():
                tc.setAllowMaxGrowthrateAtLowerCutoff(False)
            else:
                tc.setAllowMaxGrowthrateAtLowerCutoff(None)
        elif self.allowMaxGrowthrateAtLowerCutoffCheckbox.checkState() == QtCore.Qt.PartiallyChecked:
            tc.setAllowMaxGrowthrateAtLowerCutoff(None)
        self.updateParameterLineEditColor(tc,self.allowMaxGrowthrateAtLowerCutoffCheckbox,'allowMaxGrowthrateAtLowerCutoff')
        self.updatePlateFig()

    def allowGrowthyieldSlopeNStderrAwayFromZeroEdited(self,newval):
        if self.selectedReplicate is None:
            return

        self.selectedReplicate.setAllowGrowthyieldSlopeNStderrAwayFromZero(newval)
        self.updateParameterLineEditColor(self.selectedReplicate,
                                          self.allowGrowthyieldSlopeNStderrAwayFromZeroSpinbox,"allowGrowthyieldSlopeNStderrAwayFromZero")
        self.updatePlateFig()

    class PropDialog(QtGui.QDialog):
        def __init__(self,odreplicatewidgets):
            self.odreplicatewidgets=odreplicatewidgets

            # update variables denoting whether options are set for one of the plate's replicates
            hasAllowMaxGrowthrateAtLowerCutoff=odreplicatewidgets.plateHasAllowMaxGrowthrateAtLowerCutoff()
            hasAllowGrowthyieldSlopeNStderrAwayFromZero=odreplicatewidgets.plateHasAllowGrowthyieldSlopeNStderrAwayFromZero()

            super(OdReplicateWidgets.PropDialog, self).__init__()
            self.setWindowTitle('Show configuration options')

            self.showAlternativeGrowthrateExtractionCheckbox=QtGui.QCheckBox('Show alternative methods to extract the growthrate')
            if odreplicatewidgets.showAlternativeGrowthrateExtraction:
                self.showAlternativeGrowthrateExtractionCheckbox.setCheckState(QtCore.Qt.Checked)
            else:
                self.showAlternativeGrowthrateExtractionCheckbox.setCheckState(QtCore.Qt.Unchecked)

            self.showAllowMaxGrowthrateAtLowerCutoffCheckbox=QtGui.QCheckBox('Show "allow growthrate at cutoff"')
            if hasAllowMaxGrowthrateAtLowerCutoff is True:
                self.showAllowMaxGrowthrateAtLowerCutoffCheckbox.setCheckState(QtCore.Qt.Checked)
                self.showAllowMaxGrowthrateAtLowerCutoffCheckbox.setEnabled(False)
            elif odreplicatewidgets.showAllowMaxGrowthrateAtLowerCutoff is True:
                self.showAllowMaxGrowthrateAtLowerCutoffCheckbox.setCheckState(QtCore.Qt.Checked)
            else:
                self.showAllowMaxGrowthrateAtLowerCutoffCheckbox.setCheckState(QtCore.Qt.Unchecked)

            self.showAllowGrowthyieldSlopeNStderrAwayFromZeroCheckbox=QtGui.QCheckBox('Show "allow growth yield\'s slope to be zero within n standard errors"')
            if hasAllowGrowthyieldSlopeNStderrAwayFromZero is True:
                self.showAllowGrowthyieldSlopeNStderrAwayFromZeroCheckbox.setCheckState(QtCore.Qt.Checked)
                self.showAllowGrowthyieldSlopeNStderrAwayFromZeroCheckbox.setEnabled(False)
            elif odreplicatewidgets.showAllowGrowthyieldSlopeNStderrAwayFromZero is True:
                self.showAllowGrowthyieldSlopeNStderrAwayFromZeroCheckbox.setCheckState(QtCore.Qt.Checked)
            else:
                self.showAllowGrowthyieldSlopeNStderrAwayFromZeroCheckbox.setCheckState(QtCore.Qt.Unchecked)

            self.cancelButton=QtGui.QPushButton("Cancel")
            self.okButton=QtGui.QPushButton("Ok")
            self.buttonlayout=QtGui.QHBoxLayout()
            self.buttonlayout.addWidget(self.cancelButton)
            self.buttonlayout.addWidget(self.okButton)

            self.cancelButton.clicked.connect(self.close)
            self.okButton.clicked.connect(self.saveSettingsAndClose)

            self.mainlayout=QtGui.QVBoxLayout()
            self.mainlayout.addWidget(self.showAlternativeGrowthrateExtractionCheckbox)
            self.mainlayout.addWidget(self.showAllowMaxGrowthrateAtLowerCutoffCheckbox)
            self.mainlayout.addWidget(self.showAllowGrowthyieldSlopeNStderrAwayFromZeroCheckbox)
            self.mainlayout.addLayout(self.buttonlayout)

            self.setLayout(self.mainlayout)

        def saveSettingsAndClose(self):
            if self.showAlternativeGrowthrateExtractionCheckbox.checkState() == QtCore.Qt.Checked:
                self.odreplicatewidgets.showAlternativeGrowthrateExtraction=True
                self.odreplicatewidgets.setAlternativeGrowthrateExtractionVisible(True)
            else:
                self.odreplicatewidgets.showAlternativeGrowthrateExtraction=False
                self.odreplicatewidgets.setAlternativeGrowthrateExtractionVisible(False)

            if self.showAllowMaxGrowthrateAtLowerCutoffCheckbox.checkState() == QtCore.Qt.Checked:
                self.odreplicatewidgets.showAllowMaxGrowthrateAtLowerCutoff=True
                self.odreplicatewidgets.allowMaxGrowthrateAtLowerCutoffCheckbox.setVisible(True)
            else:
                self.odreplicatewidgets.showAllowMaxGrowthrateAtLowerCutoff=False
                self.odreplicatewidgets.allowMaxGrowthrateAtLowerCutoffCheckbox.setVisible(False)

            if self.showAllowGrowthyieldSlopeNStderrAwayFromZeroCheckbox.checkState() == QtCore.Qt.Checked:
                self.odreplicatewidgets.showAllowGrowthyieldSlopeNStderrAwayFromZero=True
                self.odreplicatewidgets.allowGrowthyieldSlopeNStderrAwayFromZeroGroup.setVisible(True)
            else:
                self.odreplicatewidgets.showAllowGrowthyieldSlopeNStderrAwayFromZero=False
                self.odreplicatewidgets.allowGrowthyieldSlopeNStderrAwayFromZeroGroup.setVisible(False)

            self.accept()

    class SavePdfDialog(QtGui.QDialog):
        def __init__(self):
            super(OdReplicateWidgets.SavePdfDialog, self).__init__()
            self.setWindowTitle('Save Pdf options')

            self.selectPlots=QtGui.QButtonGroup()
            self.selectPlots.setExclusive(True)
            btnSelected=QtGui.QRadioButton('Currently selected')
            btnAllGroups=QtGui.QRadioButton('All replicate groups')
            btnAllWells=QtGui.QRadioButton('All single wells')
            self.selectPlots.addButton(btnSelected)
            self.selectPlots.addButton(btnAllGroups)
            self.selectPlots.addButton(btnAllWells)
            self.selectPlots.setId(btnSelected,0)
            self.selectPlots.setId(btnAllGroups,1)
            self.selectPlots.setId(btnAllWells,2)
            selectLayout=QtGui.QVBoxLayout()
            selectLayout.addWidget(btnSelected)
            selectLayout.addWidget(btnAllGroups)
            selectLayout.addWidget(btnAllWells)
            selectGroup=QtGui.QGroupBox('Replicates to export')
            selectGroup.setLayout(selectLayout)

            self.showDerivatives=QtGui.QCheckBox("Show derivatives in second panel")
            self.includeBackground=QtGui.QCheckBox("Also show background replicates")
            self.includeBackground.setToolTip('When saving all replicate groups/all wells do not skip the background replicates.')
            self.showTitle=QtGui.QCheckBox("Add title")
            self.showTitle.setToolTip("Add a title consisting of 'sample-condition' to figures.")
            showLayout=QtGui.QVBoxLayout()
            showLayout.addWidget(self.showDerivatives)
            showLayout.addWidget(self.includeBackground)
            showLayout.addWidget(self.showTitle)
            showGroup=QtGui.QGroupBox('Figure options')
            showGroup.setLayout(showLayout)

            btnSelected.setChecked(True)
            self.includeBackground.setEnabled(False)
            self.selectPlots.buttonClicked.connect(self.onSelectPlotsChanged)

            allLayout=QtGui.QHBoxLayout()
            allLayout.addWidget(selectGroup)
            allLayout.addWidget(showGroup)
            allGroup=QtGui.QGroupBox('')
            allGroup.setLayout(allLayout)

            self.cancelButton=QtGui.QPushButton("Cancel")
            self.okButton=QtGui.QPushButton("Ok")
            self.buttonlayout=QtGui.QHBoxLayout()
            self.buttonlayout.addWidget(self.cancelButton)
            self.buttonlayout.addWidget(self.okButton)

            self.cancelButton.clicked.connect(self.reject)
            self.okButton.clicked.connect(self.accept)

            self.mainlayout=QtGui.QVBoxLayout()
            self.mainlayout.addWidget(allGroup)
            self.mainlayout.addLayout(self.buttonlayout)

            self.setLayout(self.mainlayout)

        def onSelectPlotsChanged(self,idx):
            if self.selectPlots.checkedId() == 0:
                self.includeBackground.setEnabled(False)
            else:
                self.includeBackground.setEnabled(True)

    class SaveCsvDialog(QtGui.QDialog):
        def __init__(self,deselected,selected,singleWellsNotReplicateGroups,addVarianceColumns):
            super(OdReplicateWidgets.SaveCsvDialog, self).__init__()
            self.setWindowTitle('Save Csv options')

            ###
            self.addRemoveWidget=ElementChooserWidget(deselected,selected,'Available properties:','Selected for export:')

            ###
            cancelButton=QtGui.QPushButton('Cancel')
            okButton=QtGui.QPushButton('Ok')
            okcancelbuttonlayout=QtGui.QHBoxLayout()
            okcancelbuttonlayout.addWidget(cancelButton)
            okcancelbuttonlayout.addWidget(okButton)

            cancelButton.clicked.connect(self.reject)
            okButton.clicked.connect(self.accept)

            ###
            self.singleWellsNotReplicateGroupsCombo=QtGui.QComboBox()
            self.singleWellsNotReplicateGroupsCombo.addItem('Replicate Groups')
            self.singleWellsNotReplicateGroupsCombo.addItem('Single Wells')
            if singleWellsNotReplicateGroups:
                self.singleWellsNotReplicateGroupsCombo.setCurrentIndex(1)
            else:
                self.singleWellsNotReplicateGroupsCombo.setCurrentIndex(0)
            self.singleWellsNotReplicateGroupsCombo.currentIndexChanged.connect(self.onSingleWellsNotReplicateGroupsChanged)

            ###
            self.varianceCheckbox=QtGui.QCheckBox('Export with variance for each property')
            if addVarianceColumns:
                self.varianceCheckbox.setCheckState(QtCore.Qt.Checked)
            else:
                self.varianceCheckbox.setCheckState(QtCore.Qt.Unchecked)
            if singleWellsNotReplicateGroups:
                self.varianceCheckbox.setEnabled(False)

            ###
            maingroup=QtGui.QGroupBox('')
            mainlayout=QtGui.QVBoxLayout()
            mainlayout.addWidget(self.addRemoveWidget)
            mainlayout.addWidget(self.singleWellsNotReplicateGroupsCombo)
            mainlayout.addWidget(self.varianceCheckbox)
            mainlayout.addLayout(okcancelbuttonlayout)

            self.setLayout(mainlayout)

        def onSingleWellsNotReplicateGroupsChanged(self,idx):
            if idx == 1:
                self.varianceCheckbox.setEnabled(False)
            else:
                self.varianceCheckbox.setEnabled(True)

    class SaveTimeseriesDialog(QtGui.QDialog):
        def __init__(self,singleWellsNotReplicateGroups,addVarianceColumns):
            super(OdReplicateWidgets.SaveTimeseriesDialog, self).__init__()
            self.setWindowTitle('Save time series options')

            ###
            cancelButton=QtGui.QPushButton('Cancel')
            okButton=QtGui.QPushButton('Ok')
            okcancelbuttonlayout=QtGui.QHBoxLayout()
            okcancelbuttonlayout.addWidget(cancelButton)
            okcancelbuttonlayout.addWidget(okButton)

            cancelButton.clicked.connect(self.close)
            okButton.clicked.connect(self.accept)

            ###
            self.singleWellsNotReplicateGroupsCombo=QtGui.QComboBox()
            self.singleWellsNotReplicateGroupsCombo.addItem('Replicate Groups')
            self.singleWellsNotReplicateGroupsCombo.addItem('Single Wells')
            if singleWellsNotReplicateGroups:
                self.singleWellsNotReplicateGroupsCombo.setCurrentIndex(1)
            else:
                self.singleWellsNotReplicateGroupsCombo.setCurrentIndex(0)
            self.singleWellsNotReplicateGroupsCombo.currentIndexChanged.connect(self.onSingleWellsNotReplicateGroupsChanged)

            ###
            self.varianceCheckbox=QtGui.QCheckBox('Export with variance for each property')
            if addVarianceColumns:
                self.varianceCheckbox.setCheckState(QtCore.Qt.Checked)
            else:
                self.varianceCheckbox.setCheckState(QtCore.Qt.Unchecked)
            if singleWellsNotReplicateGroups:
                self.varianceCheckbox.setEnabled(False)

            ###
            maingroup=QtGui.QGroupBox('')
            mainlayout=QtGui.QVBoxLayout()
            mainlayout.addWidget(self.singleWellsNotReplicateGroupsCombo)
            mainlayout.addWidget(self.varianceCheckbox)
            mainlayout.addLayout(okcancelbuttonlayout)

            self.setLayout(mainlayout)

        def onSingleWellsNotReplicateGroupsChanged(self,idx):
            if idx == 1:
                self.varianceCheckbox.setEnabled(False)
            else:
                self.varianceCheckbox.setEnabled(True)

    class SaveMetadataTemplateDialog(QtGui.QDialog):
        def __init__(self):
            super(OdReplicateWidgets.SaveMetadataTemplateDialog, self).__init__()
            self.setWindowTitle('Save metadata template')

            self.selectPlateformat=QtGui.QButtonGroup()
            self.selectPlateformat.setExclusive(True)
            btn096=QtGui.QRadioButton('96 wells')
            btn384=QtGui.QRadioButton('384 wells')
            btn100honeycomb=QtGui.QRadioButton('100 wells honeycomb')
            btn200honeycomb=QtGui.QRadioButton('200 wells honeycomb')
            btn096.setChecked(True)
            self.selectPlateformat.addButton(btn096)
            self.selectPlateformat.addButton(btn384)
            self.selectPlateformat.addButton(btn100honeycomb)
            self.selectPlateformat.addButton(btn200honeycomb)
            self.selectPlateformat.setId(btn096,0)
            self.selectPlateformat.setId(btn384,1)
            self.selectPlateformat.setId(btn100honeycomb,2)
            self.selectPlateformat.setId(btn200honeycomb,3)
            selectLayout=QtGui.QVBoxLayout()
            selectLayout.addWidget(btn096)
            selectLayout.addWidget(btn384)
            selectLayout.addWidget(btn100honeycomb)
            selectLayout.addWidget(btn200honeycomb)

            okButton=QtGui.QPushButton('Ok')
            cancelButton=QtGui.QPushButton('Cancel')
            okButton.setDefault(True)
            cancelButton.clicked.connect(self.close)
            okButton.clicked.connect(self.accept)
            okcancelbuttonlayout=QtGui.QHBoxLayout()
            okcancelbuttonlayout.addWidget(cancelButton)
            okcancelbuttonlayout.addWidget(okButton)

            ###
            maingroup=QtGui.QGroupBox('')
            mainlayout=QtGui.QVBoxLayout()
            mainlayout.addLayout(selectLayout)
            mainlayout.addLayout(okcancelbuttonlayout)

            self.setLayout(mainlayout)

def getSaveFileNameDialogWithDefaultSuffix(parent, caption, directory, filters, defaultSuffix=None):
    fileDialog=QtGui.QFileDialog(parent, caption, directory, filters)
    fileDialog.setAcceptMode(QtGui.QFileDialog.AcceptSave)
    fileDialog.setFileMode(QtGui.QFileDialog.AnyFile)
    if defaultSuffix is not None:
        fileDialog.setDefaultSuffix(defaultSuffix)
    if fileDialog.exec_():
        return fileDialog.selectedFiles()[0]
    return Qt.QString()

class ODMainWindow(QtGui.QMainWindow):

    def __init__(self):
        super(ODMainWindow, self).__init__()

        self.defaultLogOdCutoff=None
        self.defaultLagAtLogOdEquals=None
        self.defaultMaxGrowthLowerTimeCutoff=None
        self.defaultMaxGrowthUpperTimeCutoff=None
        self.defaultHdCorrectionLinear=1.
        self.defaultHdCorrectionQuadratic=0.
        self.defaultHdCorrectionCubic=0.
        self.defaultSmoothingK=5
        self.defaultSmoothingS=0.01
        # which columns to show when exporting csv table
        columns, morecolumns=Plate.availableColumnsForCsvExport()
        self.presetColumns=['doublingtime_expfit',
                            'growthrate_expfit',
                            'lag_expfit',
                            'yield',
                            'wellids']
        self.addVarianceColumns=True
        self.singleWellsNotReplicateGroups=False
        self.initUI()
        
    def initUI(self):
        self.openAction = QtGui.QAction('&Open', self)
        self.openAction.setShortcut(QtGui.QKeySequence.Open)
        self.openAction.setStatusTip('Open file.')
        self.openAction.triggered.connect(self.openODfile)

        self.loadMetadataAction = QtGui.QAction('Load metadata', self)
        self.loadMetadataAction.setStatusTip('Load metadata (sample, condition) from csv.')
        self.loadMetadataAction.triggered.connect(self.loadMetadata)

        self.saveAsAction = QtGui.QAction('Save as', self)
        self.saveAsAction.setShortcut(QtGui.QKeySequence.SaveAs)
        self.saveAsAction.setStatusTip('Save data.')
        self.saveAsAction.triggered.connect(self.saveAs)

        self.saveAction = QtGui.QAction('Save', self)
        self.saveAction.setShortcut(QtGui.QKeySequence.Save)
        self.saveAction.setStatusTip('Save data.')
        self.saveAction.triggered.connect(self.save)

        self.saveMetadataAction = QtGui.QAction('Save metadata', self)
        self.saveMetadataAction.setStatusTip('Save metadata (sample, condition) as csv.')
        self.saveMetadataAction.triggered.connect(self.saveMetadata)

        self.savePdfAction = QtGui.QAction('Export &figure', self)
        self.savePdfAction.setStatusTip('Save a figure containing growth curves and their derivatives.')
        self.savePdfAction.triggered.connect(self.savePdf)

        self.saveCsvAction = QtGui.QAction('Export &properties', self)
        self.saveCsvAction.setStatusTip('Save a table containing growth properties such as growth rate, yield, etc.')
        self.saveCsvAction.triggered.connect(self.saveCsv)

        self.saveTimeseriesAction = QtGui.QAction('Export &time series', self)
        self.saveTimeseriesAction.setStatusTip('Save a table containing time series.')
        self.saveTimeseriesAction.triggered.connect(self.saveTimeseries)

        # alt-return or shift-Ctrl-D
        self.propertiesAction = QtGui.QAction('Settings', self)
        self.propertiesAction.setStatusTip('Enable or disable advanced features.')
        self.propertiesAction.triggered.connect(self.showProperties)

        self.closeAction = QtGui.QAction('&Close', self)
        self.closeAction.setShortcut(QtGui.QKeySequence.Close)
        self.closeAction.setStatusTip('Close file.')
        self.closeAction.triggered.connect(self.closeODfile)

        self.exitAction = QtGui.QAction('&Exit', self)
        self.exitAction.setShortcut(QtGui.QKeySequence.Quit)
        self.exitAction.setStatusTip('Exit.')
        self.exitAction.triggered.connect(self.quit)

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(self.openAction)
        fileMenu.addAction(self.loadMetadataAction)
        fileMenu.addAction(self.saveAsAction)
        fileMenu.addAction(self.saveAction)
        fileMenu.addAction(self.saveMetadataAction)
        fileMenu.addAction(self.savePdfAction)
        fileMenu.addAction(self.saveCsvAction)
        fileMenu.addAction(self.saveTimeseriesAction)
        fileMenu.addAction(self.propertiesAction)
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

        self.platewidgets=OdReplicateWidgets()
        self.platewidgets.statusbar=self.statusBar()

        parametersToolbar=self.addToolBar('Parameters')
        parametersToolbar.addWidget(self.platewidgets.maxGrowthTimeCutoffGroup)
        # NOTE this is a quick hack: visibility of a widget in a
        # toolbar has to be changed from the corresponding QAction
        # (see qtoolbar.html#addWidget)
        self.platewidgets.allowGrowthyieldSlopeNStderrAwayFromZeroGroup = parametersToolbar.addWidget(self.platewidgets.allowGrowthyieldSlopeNStderrAwayFromZeroGroup)
        parametersToolbar.addWidget(self.platewidgets.hdCorrectionGroup)
        parametersToolbar.addWidget(self.platewidgets.cutoffGroup)
        parametersToolbar.addWidget(self.platewidgets.lagAtLogOdEqualsGroup)
        parametersToolbar.addWidget(self.platewidgets.slidingWindowSizeGroup)
        parametersToolbar.addWidget(self.platewidgets.smoothingGroup)
        parametersToolbar.addWidget(self.platewidgets.ntb)

        self.addToolBarBreak()
        plotoptsToolbar=self.addToolBar('Plot options')
        plotoptsToolbar.addWidget(self.platewidgets.showRawCheck)
        plotoptsToolbar.addWidget(self.platewidgets.showSingleCheck)
        plotoptsToolbar.addWidget(self.platewidgets.showBackgroundCheck)
        plotoptsToolbar.addWidget(self.platewidgets.showSmoothedCheck)
        plotoptsToolbar.addWidget(self.platewidgets.showLogCheck)
        plotoptsToolbar.addWidget(self.platewidgets.showLogSmoothedCheck)
        plotoptsToolbar.addWidget(self.platewidgets.showExpFitsOd0MuCheck)
        plotoptsToolbar.addWidget(self.platewidgets.showMaxGrowthrateCheck)
        plotoptsToolbar.addWidget(self.platewidgets.showGrowthyieldCheck)
        # saving actions here, see below at setAlternativeGrowthrateExtractionVisible why
        self.showLogDerNonLogAction=plotoptsToolbar.addWidget(self.platewidgets.showLogDerNonLogCheck)
        self.showLogDerNonLogSmoothedAction=plotoptsToolbar.addWidget(self.platewidgets.showLogDerNonLogSmoothedCheck)
        self.showMaxGrowthrateFromLogOdAction=plotoptsToolbar.addWidget(self.platewidgets.showMaxGrowthrateFromLogOdCheck)

        self.graphsnchooser = QtGui.QSplitter(QtCore.Qt.Horizontal,self.mainwidget)
        mainwidgetLayout.addWidget(self.graphsnchooser)
        self.graphsnchooser.setSizePolicy(Qt.QSizePolicy.Expanding,
                                          Qt.QSizePolicy.Expanding)

        self.graphsnchooser.addWidget(self.platewidgets.odplateView)
        statusngraphs = QtGui.QSplitter(QtCore.Qt.Vertical,self.graphsnchooser)
        statusngraphs.addWidget(self.platewidgets.mainstack)
        self.graphsnchooser.setStretchFactor(0,1)
        self.graphsnchooser.setStretchFactor(1,4)


        self.setGeometry(100, 100, 1024, 600)

        self.setAlternativeGrowthrateExtractionVisible(False)
        self.clear()
        self.show()

    def setAlternativeGrowthrateExtractionVisible(self,visible):
        """
        NOTE addWidget documentation of QToolBar says:
        You should use QAction::setVisible() to change the visibility
        of the widget. Using QWidget::setVisible(), QWidget::show()
        and QWidget::hide() does not work.
        """
        self.platewidgets.setAlternativeGrowthrateExtractionVisible(visible)
        self.showLogDerNonLogAction.setVisible(visible)
        self.showLogDerNonLogSmoothedAction.setVisible(visible)
        self.showMaxGrowthrateFromLogOdAction.setVisible(visible)

    def clear(self):
        self.openAction.setEnabled(True)
        self.loadMetadataAction.setEnabled(False)
        self.saveAsAction.setEnabled(False)
        self.saveAction.setEnabled(False)
        self.saveMetadataAction.setEnabled(True)
        self.saveMetadataAction.setText('Save metadata template')
        self.saveMetadataAction.setStatusTip('Save metadata template as csv.')
        self.closeAction.setEnabled(False)
        self.savePdfAction.setEnabled(False)
        self.saveCsvAction.setEnabled(False)
        self.saveTimeseriesAction.setEnabled(False)
        self.propertiesAction.setEnabled(False)
        self.platewidgets.clear()
        self.saveFilename=None
        self.importFilename=None
        self._setWindowTitle()

    def enDisableActionsAfterLoading(self):
        self.openAction.setEnabled(False)
        self.loadMetadataAction.setEnabled(True)
        self.saveAsAction.setEnabled(True)
        self.saveAction.setEnabled(True)
        self.saveMetadataAction.setEnabled(True)
        self.saveMetadataAction.setText('Save metadata')
        self.saveMetadataAction.setStatusTip('Save metadata as csv.')
        self.savePdfAction.setEnabled(True)
        self.saveCsvAction.setEnabled(True)
        self.saveTimeseriesAction.setEnabled(True)
        self.propertiesAction.setEnabled(True)
        self.closeAction.setEnabled(True)
        self.exitAction.setEnabled(True)

    def _setWindowTitle(self):
        wtitle=''
        if self.platewidgets.odplate is not None:
            if self.saveFilename is not None:
                wtitle=os.path.basename(self.saveFilename)+' - '
            elif self.importFilename is not None:
                wtitle=os.path.basename(self.importFilename)+' - '
            else:
                wtitle=self.platewidgets.odplate.plateId
        self.setWindowTitle(wtitle+'Growth Analysis Tool (GATHODE)')

    def openODfile(self):
        self.closeODfile()
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '', 'Growth analysis data (*.asc *.txt *.gat)')
        if len(fname) == 0:
            return
        condition=""
        self._openODfile(str(fname))

    def _openODfile(self,fname):
        try:
            odplate_=Plate(filename=fname)
            if odplate_.readfileformat == 'gat':
                self.saveFilename=fname
                self.saveAction.setEnabled(True)
            else:
                self.importFilename=fname
                odplate_.setHighDensityCorrectionLinear(self.defaultHdCorrectionLinear)
                odplate_.setHighDensityCorrectionQuadratic(self.defaultHdCorrectionQuadratic)
                odplate_.setHighDensityCorrectionCubic(self.defaultHdCorrectionCubic)
                odplate_.setSmoothingS(self.defaultSmoothingS)
                odplate_.setSmoothingK(self.defaultSmoothingK)
                if self.defaultLogOdCutoff is not None:
                    odplate_.setLogOdCutoff(self.defaultLogOdCutoff)
                if self.defaultMaxGrowthLowerTimeCutoff is not None:
                    odplate_.setMaxGrowthLowerTimeCutoff(self.defaultMaxGrowthLowerTimeCutoff)
                if self.defaultMaxGrowthUpperTimeCutoff is not None:
                    odplate_.setMaxGrowthUpperTimeCutoff(self.defaultMaxGrowthUpperTimeCutoff)
                if self.defaultLagAtLogOdEquals is not None:
                    odplate_.setLagAtLogOdEquals(self.defaultLagAtLogOdEquals)
                odplate_.modified=False # needed as one of above set* methods might have set this to True
        except (Plate.Error, IOError) as err:
            QtGui.QMessageBox.warning(self,'Error reading datafile',str(err))
            odplate_=None
            return
        except Exception as err:
            msg=str(err)+"\n"
            if msg is None:
                msg=''
            exc_type, exc_value, exc_tb = sys.exc_info()
            for line in traceback.format_exception(exc_type, exc_value, exc_tb):
                msg+=line
            QtGui.QMessageBox.warning(self,'Unknown error',msg)
            odplate_=None
            return

        self.enDisableActionsAfterLoading()
        msg=odplate_._loadStatus.message()
        self.platewidgets.setPlate(odplate_)
        self._setWindowTitle()
        if msg != '':
            QtGui.QMessageBox.warning(self,'Warning',msg+'.')

    def saveAs(self):
        if self.platewidgets.odplate is None:
            return
        fname = str(getSaveFileNameDialogWithDefaultSuffix(self, 'Save file', '', 'Growth analysis data (*.gat)', 'gat'))
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
            savestatus=self.platewidgets.odplate.save(self.saveFilename)
            if savestatus is not None:
                QtGui.QMessageBox.warning(self,savestatus.messageType(),savestatus.longmessage())
            return True
        except Plate.Error as err:
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

    def saveMetadata(self):
        fname = str(getSaveFileNameDialogWithDefaultSuffix(self, 'Save metadata', '', 'GATHODE metadata (*.csv)', 'csv'))
        if len(fname) == 0:
            return
        try:
            if self.platewidgets.odplate is None:
                dialog=OdReplicateWidgets.SaveMetadataTemplateDialog()
                if not dialog.exec_():
                    return
                if dialog.selectPlateformat.checkedId() == 0:
                    numWells=96
                elif dialog.selectPlateformat.checkedId() == 1:
                    numWells=384
                elif dialog.selectPlateformat.checkedId() == 2:
                    numWells=100
                elif dialog.selectPlateformat.checkedId() == 3:
                    numWells=200
                Plate.writeMetadata(fname,[{} for i in range(numWells)],['sample','condition'],
                                    plateformat=Plate._numWellsToFormatString(numWells))
            else:
                metadata=self.platewidgets.odplate.wellMetadata()
                Plate.writeMetadata(fname,metadata,['sample','condition'],
                                    plateformat=Plate._numWellsToFormatString(len(metadata)))
            return True
        except Plate.Error as err:
            QtGui.QMessageBox.warning(self,'Error saving metadata',str(err))
        except Exception as err:
            msg=str(err)+"\n"
            if msg is None:
                msg=''
            exc_type, exc_value, exc_tb = sys.exc_info()
            for line in traceback.format_exception(exc_type, exc_value, exc_tb):
                msg+=line
            QtGui.QMessageBox.warning(self,'Unknown error',msg)
        return False

    def loadMetadata(self):
        if self.platewidgets.odplate is None:
            return
        fname = str(QtGui.QFileDialog.getOpenFileName(self, 'Open file', '', 'GATHODE metadata (*.csv)'))
        if len(fname) == 0:
            return
        try:
            metadata=Plate.readMetadata(fname,Plate._numWellsToFormatString(len(self.platewidgets.odplate.wellMetadata())))
            metok, metstatus = self.platewidgets.odplate.wellMetadataOk(metadata)
            if not metok:
                raise Plate.BadMetadata(metstatus.message(),filename=fname)

            ret=QtGui.QMessageBox.question(self,"Overwriting metadata",
                                           "You are overwriting the metadata (sample ids and conditions).\n"
                                           +"Do you want to proceed?",
                                           QtGui.QMessageBox.Ok | QtGui.QMessageBox.Cancel,
                                           QtGui.QMessageBox.Cancel)
            if ret == QtGui.QMessageBox.Cancel:
                return False

            self.platewidgets.odplate.setWellMetadata(metadata)
            self.platewidgets.setPlate(self.platewidgets.odplate)
            self._setWindowTitle()
            msg=self.platewidgets.odplate._loadStatus.message()
            if msg != '':
                QtGui.QMessageBox.warning(self,'Warning',msg+'.')
            return True
        except Plate.Error as err:
            QtGui.QMessageBox.warning(self,'Error loading metadata',str(err))
        except Exception as err:
            msg=str(err)+"\n"
            if msg is None:
                msg=''
            exc_type, exc_value, exc_tb = sys.exc_info()
            for line in traceback.format_exception(exc_type, exc_value, exc_tb):
                msg+=line
            QtGui.QMessageBox.warning(self,'Unknown error loading metadata',msg)
        return False

    def closeODfile(self):
        if self.platewidgets.odplate is not None and self.platewidgets.odplate.modified:
            ret=QtGui.QMessageBox.question(self,"Growth Analysis Tool (GATHODE)",
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

        if self.platewidgets.odplate is not None:
            self.platewidgets.disconnect()
        self.clear()
        return True

    def quit(self):
        if self.closeODfile():
            QtGui.qApp.quit()

    def savePdf(self):
        if self.platewidgets.odplate is None:
            return
        dialog=OdReplicateWidgets.SavePdfDialog()
        if not dialog.exec_():
            return
        if dialog.selectPlots.checkedId() == 0:
            if isinstance(self.platewidgets.selectedReplicate, Plate):
                QtGui.QMessageBox.warning(self,'No replicate selected','When choosing "Currently selected" you should have selected a replicate.')
                return
        fname = str(getSaveFileNameDialogWithDefaultSuffix(self, 'Save file', '', 'pdf (*.pdf)', 'pdf'))
        if len(fname) == 0:
            return
        plotopts={
            'showTitle': (dialog.showTitle.checkState() == QtCore.Qt.Checked),
            'creator': os.path.basename(sys.argv[0]),
            'showDerivatives': (dialog.showDerivatives.checkState() == QtCore.Qt.Checked),
            'includeBackground': (dialog.includeBackground.checkState() == QtCore.Qt.Checked),
            'showRaw': (self.platewidgets.showRawCheck.checkState() == QtCore.Qt.Checked),
            'showSingle': (self.platewidgets.showSingleCheck.checkState() == QtCore.Qt.Checked),
            'showBackground': (self.platewidgets.showBackgroundCheck.checkState() == QtCore.Qt.Checked),
            'showSmoothed': (self.platewidgets.showSmoothedCheck.checkState() == QtCore.Qt.Checked),
            'showSmoothedDerivativeLinear': (self.platewidgets.showSmoothedCheck.checkState() == QtCore.Qt.Checked),
            'showLogod': (self.platewidgets.showLogCheck.checkState() == QtCore.Qt.Checked),
            'showLogodSmoothed': (self.platewidgets.showLogSmoothedCheck.checkState() == QtCore.Qt.Checked),
            'showMaxGrowthrate': (self.platewidgets.showMaxGrowthrateCheck.checkState() == QtCore.Qt.Checked),
            'showMaxGrowthrateFromLogOdDerivative': (self.platewidgets.showMaxGrowthrateFromLogOdCheck.checkState() == QtCore.Qt.Checked),
            'showLogOdDerivative': (self.platewidgets.showLogDerCheck.checkState() == QtCore.Qt.Checked),
            'showLogOdDerivativeFromNonLog': (self.platewidgets.showLogDerNonLogCheck.checkState() == QtCore.Qt.Checked),
            'showLogOdDerivativeFromNonLogSmoothed': (self.platewidgets.showLogDerNonLogSmoothedCheck.checkState()
                                                   == QtCore.Qt.Checked),
            'showExpFitsOd0Mu': (self.platewidgets.showExpFitsOd0MuCheck.checkState() == QtCore.Qt.Checked),
            'showGrowthyield': (self.platewidgets.showGrowthyieldCheck.checkState() == QtCore.Qt.Checked),
            }
        if dialog.selectPlots.checkedId() == 0:
            # selected figure
            del plotopts['includeBackground']
            replicateNderivativePlotAsPdf(fname,self.platewidgets.selectedReplicate,**plotopts)
        elif dialog.selectPlots.checkedId() == 1:
            # all replicate groups
            plotopts['showReplicateGroups']=True
            if dialog.includeBackground.checkState() == QtCore.Qt.Checked:
                numPlots=len(self.platewidgets.odplate.replicateGroups)
            else:
                numPlots=len(self.platewidgets.odplate.nonBackgroundReplicates())
            self._setupLongrunningCalculation('Saving figures',numPlots)
            try:
                plotFullOdPlate(self.platewidgets.odplate,pdfout=fname,progressCall=self._updateProgress,**plotopts)
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
        elif dialog.selectPlots.checkedId() == 2:
            # all single wells
            plotopts['showReplicateGroups']=False
            plotopts['showWellString']=True
            if dialog.includeBackground.checkState() == QtCore.Qt.Checked:
                numPlots=len(self.platewidgets.odplate.wells)
            else:
                numPlots=len(self.platewidgets.odplate.nonBackgroundWells())
            self._setupLongrunningCalculation('Saving figures',numPlots)
            try:
                plotFullOdPlate(self.platewidgets.odplate,pdfout=fname,progressCall=self._updateProgress,**plotopts)
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

    def saveCsv(self):
        if self.platewidgets.odplate is None:
            return

        # in presetColumns is saved which columns to show when exporting csv table
        # here we generate the list of the rest of columns
        columns, morecolumns=Plate.availableColumnsForCsvExport(self.platewidgets.showAlternativeGrowthrateExtraction)
        available=[]
        for col in morecolumns:
            if col not in self.presetColumns:
                available.append(col)
        # we recalculate presetColumns here because showAlternativeGrowthrateExtraction may have changed
        newPresetColumns=[]
        for col in self.presetColumns:
            if col in morecolumns:
                newPresetColumns.append(col)
            self.presetColumns=newPresetColumns

        dialog=OdReplicateWidgets.SaveCsvDialog(available,self.presetColumns,self.singleWellsNotReplicateGroups,self.addVarianceColumns)
        if not dialog.exec_():
            return
        self.presetColumns=[]
        for i in range(dialog.addRemoveWidget.selected.count()):
            el = str(dialog.addRemoveWidget.selected.item(i).data(QtCore.Qt.DisplayRole))
            self.presetColumns.append(el)
            columns.append(el)
        self.addVarianceColumns=dialog.varianceCheckbox.checkState() == QtCore.Qt.Checked
        self.singleWellsNotReplicateGroups=dialog.singleWellsNotReplicateGroupsCombo.currentIndex() == 1

        fname = str(getSaveFileNameDialogWithDefaultSuffix(self, 'Save file', '', 'csv (*.csv)', 'csv'))
        if len(fname) == 0:
            return

        if self.singleWellsNotReplicateGroups:
            numForProgress=len(self.platewidgets.odplate.wells)
        else:
            numForProgress=len(self.platewidgets.odplate.replicateGroups)
        self._setupLongrunningCalculation('Saving table',numForProgress)
        try:
            self.platewidgets.odplate.growthParametersToCsv(fname,progressCall=self._updateProgress,
                                                            columns=columns,addVarianceColumns=self.addVarianceColumns,
                                                            singleWells=self.singleWellsNotReplicateGroups)
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

    def saveTimeseries(self):
        if self.platewidgets.odplate is None:
            return
        dialog=OdReplicateWidgets.SaveTimeseriesDialog(self.singleWellsNotReplicateGroups,self.addVarianceColumns)
        if not dialog.exec_():
            return
        self.addVarianceColumns=dialog.varianceCheckbox.checkState() == QtCore.Qt.Checked
        self.singleWellsNotReplicateGroups=dialog.singleWellsNotReplicateGroupsCombo.currentIndex() == 1

        fname = str(getSaveFileNameDialogWithDefaultSuffix(self, 'Save file', '', 'csv (*.csv)', 'csv'))
        if len(fname) == 0:
            return
        if self.singleWellsNotReplicateGroups:
            numForProgress=len(self.platewidgets.odplate.wells)
        else:
            numForProgress=len(self.platewidgets.odplate.replicateGroups)
        self._setupLongrunningCalculation('Saving time series',numForProgress)
        try:
            self.platewidgets.odplate.timeseriesToCsv(fname,progressCall=self._updateProgress,
                                                      addVarianceColumns=self.addVarianceColumns,
                                                      singleWells=self.singleWellsNotReplicateGroups)
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

    def _setupLongrunningCalculation(self,msg,numForProgress):
        self.setEnabled(False)
        QtGui.QApplication.setOverrideCursor(Qt.QCursor(QtCore.Qt.WaitCursor))
        # create a progress that is shown in the main stack (where usually the plot is shown)
        self.progresswidget=MyProgressWidget(msg,numForProgress)
        self._mainstackIndexSaved=self.platewidgets.mainstack.currentIndex()
        self.platewidgets.mainstack.addWidget(self.progresswidget)
        self.platewidgets.mainstack.setCurrentIndex(self.platewidgets.mainstack.count()-1)
        self.platewidgets.statusbar.setVisible(False)

    def _longrunningCalculationFinished(self):
        self.platewidgets.statusbar.setVisible(True)
        self.platewidgets.mainstack.setCurrentIndex(self._mainstackIndexSaved)
        self.platewidgets.mainstack.removeWidget(self.progresswidget)
        self.progresswidget=None
        self.setEnabled(True)
        self.enDisableActionsAfterLoading()
        QtGui.QApplication.restoreOverrideCursor()

    def _updateProgress(self,intval):
        self.progresswidget.updateProgress(intval)
        QtGui.QApplication.processEvents()

    def showProperties(self):
        dialog=OdReplicateWidgets.PropDialog(self.platewidgets)
        if dialog.exec_():
            self.setAlternativeGrowthrateExtractionVisible(self.platewidgets.showAlternativeGrowthrateExtraction)

    def showAbout(self):
        msgBox=QtGui.QMessageBox.about(self,'Growth Analysis Tool (GATHODE)',
                                       'Growth Analysis Tool for High-throughput Optical Density Experiments (GATHODE)\nversion '
                                       +__version__
                                       +'\nCopyright (C) 2014 Nils Christian\n'+
                                       agpllicense)

class MyProgressWidget(QtGui.QWidget):
    def __init__(self,title='',maximumValue=0):
        super(MyProgressWidget, self).__init__()

        mainlayout=QtGui.QVBoxLayout()

        #self.setWindowTitle(title)
        self.progress=QtGui.QProgressBar(self)
        mainlayout.addStretch()
        mainlayout.addWidget(Qt.QLabel(title))
        mainlayout.addWidget(self.progress)
        mainlayout.addStretch()
        mainlayout.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
        self.progress.setMinimum(0)
        self.progress.setMaximum(maximumValue)
        self.progress.setMinimumWidth(200)

        self.setLayout(mainlayout)

    def updateProgress(self, value):
        if self.progress is not None:
            self.progress.setValue(value)

def gui_main():
    commandline=' '.join(sys.argv)
    executablename = os.path.basename(sys.argv[0])
    parser = argparse.ArgumentParser(description=executablename+': GATHODE-GUI (Growth Analysis Tool).')
    parser.add_argument('--version', action='store_true', default=None, help='show version and exit')
    parser.add_argument('infile', metavar='file.gat|export.asc', action='store', default=None, nargs='?',
                        help='file to be loaded (gat-file or TECAN ASCII format)')
    parser.add_argument('--hdlin', action='store', type=float, default=1.,
                        help='high density correction linear term')
    parser.add_argument('--hdquad', action='store', type=float, default=0.,
                        help='high density correction quadratic term')
    parser.add_argument('--hdcub', action='store', type=float, default=0.,
                        help='high density correction cubic term')
    parser.add_argument('--logodcutoff', action='store', type=float, default=-5,
                        help='cutoff applied to ln(OD) when calculating various observables')
    parser.add_argument('--maxgrowthlowertimecutoff', action='store', type=float, default=1.,
                        help='lower cutoff applied to t when calculating maximal growth rate')
    parser.add_argument('--maxgrowthuppertimecutoff', action='store', type=float, default=None,
                        help='upper cutoff applied to t when calculating maximal growth rate')
    parser.add_argument('--lagatlogodequals', action='store', type=float, default=-5,
                        help='the lag will be defined as the intersection point of the linear equation '
                        +'of the maximal growth with the given ln(OD) value')
    args = parser.parse_args()

    if args.version:
        print(executablename+' '+__version__)
        sys.exit(0)

    app = QtGui.QApplication(sys.argv)
    qod=ODMainWindow()
    qod.defaultLogOdCutoff=args.logodcutoff
    qod.defaultMaxGrowthLowerTimeCutoff=args.maxgrowthlowertimecutoff
    qod.defaultMaxGrowthUpperTimeCutoff=args.maxgrowthuppertimecutoff
    qod.defaultLagAtLogOdEquals=args.lagatlogodequals
    qod.defaultHdCorrectionLinear=args.hdlin
    qod.defaultHdCorrectionQuadratic=args.hdquad
    qod.defaultHdCorrectionCubic=args.hdcub

    if args.infile is not None:
        qod._openODfile(args.infile)

    sys.exit(app.exec_())

if __name__ == '__main__':
    gui_main()
