"""
This module implements the plot functions for CATHODE.

Chronological life span Analysis Tool for High-throughput Optical
Density Experiments (CATHODE) plot functions.
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

import re
import math

import numpy
import matplotlib
import textwrap

import contextlib
import matplotlib.backends.backend_pdf

from platereader.plate import Plate
from platereader.odplot import twentysixColours
from platereader.statusmessage import StatusMessage, Severity
from platereader.numpytools import nonNanSqrt

def viabilityToMatplotlib(repl,fig,showTitle=False,addWellIdsToTitle=False,verbose=False,color=None,viabilityInPercent=True):
    """
    Show viabilities of a ClsReplicate on a matplotlib figure.

    FIXME enhance documentation of parameters.

    :return: StatusMessage -- Statuses of this dataset.
    """
    legendprop={'size': 0}
    ax=fig.add_subplot(1,1,1)
    #ax.tick_params(bottom=False,top=False,left=False,right=False,labelbottom=False,labeltop=False,labelleft=False,labelright=False)
    ax.set_xlabel("time/days")
    if viabilityInPercent:
        ax.set_ylabel("viability (%)")
    else:
        ax.set_ylabel("viability")

    ## title
    if showTitle:
        title=repl.sampleid+" "+repl.condition
        if addWellIdsToTitle:
            title+=' '+repl.activeChildWellIdStr()
        title="\n".join(textwrap.wrap(title, 100))
        ax.set_title(title)

    days, viability, viabilityvar, status = repl.viability()
    if verbose:
        print('viability   '+str(viability))
    if viability is not None:
        if viabilityInPercent:
            viability = viability * 100.
        if viabilityvar is not None:
            viabilityerr = nonNanSqrt(viabilityvar)
            if viabilityInPercent:
                viabilityerr = viabilityerr * 100.
            valpluserr=numpy.nan_to_num(viability) + numpy.nan_to_num(viabilityerr)
            ymax=numpy.max(numpy.nan_to_num(viability) + numpy.nan_to_num(viabilityerr))
        else:
            viabilityerr=None
            ymax=numpy.nanmax(viability)
        if verbose:
            print('viabilityerr '+str(viabilityerr))
            print('ymax '+str(ymax))
        if viabilityInPercent:
            ax.axis((0,days[-1]+1,-10,ymax+10))
        else:
            ax.axis((0,days[-1]+1,-0.1,ymax+.1))
        thecolor=color if color is not None else 'red'
        ax.errorbar(days,viability,
                    yerr=viabilityerr if viabilityerr is not None else None,
                    marker='+',markeredgewidth=4,color=thecolor)
        ax.fill_between(days, viability, numpy.zeros([len(viability)]),
                        alpha=0.2,edgecolor=thecolor,facecolor=thecolor)

    return status

def viabilitiesToPdf(clsplate,pdfout,
                     replicateGroupIndices=[],elementaryIndices=[],
                     creator=None,
                     showTitle=True, addWellIdsToTitle=True,
                     progressCall=None):
    """
    Create a multi-page pdf with viability plots.
    """
    if replicateGroupIndices == [] and elementaryIndices == []:
        # nothing given, we will create a pdf containing all non-background replicate groups
        replicateGroupIndices=clsplate.nonBackgroundClsIndices()

    num=len(elementaryIndices)+len(replicateGroupIndices)

    with contextlib.closing(matplotlib.backends.backend_pdf.PdfPages(pdfout)) as pdfp:
        figall = matplotlib.figure.Figure()
        canvas = matplotlib.backends.backend_pdf.FigureCanvasPdf(figall)

        allcnt=-1

        # elementary wells
        pdftitle=''
        for clstcidx in elementaryIndices:
            allcnt+=1
            if progressCall is not None:
                progressCall(allcnt)
    
            figclstc = matplotlib.figure.Figure()
            canvas = matplotlib.backends.backend_pdf.FigureCanvasPdf(figclstc) # NOTE this needs to be called to set the the canvas of figclstc
            status = viabilityToMatplotlib(clsplate.clsWells[clstcidx],figclstc,showTitle=showTitle,addWellIdsToTitle=addWellIdsToTitle)
            figclstc.savefig(pdfp,format="pdf",bbox_inches='tight')
            pdftitle+=clsplate.clsWells[clstcidx].sampleid+" "+clsplate.clsWells[clstcidx].condition    

        # cls replicate grooups
        for clstcidx in replicateGroupIndices:
            allcnt+=1
            if progressCall is not None:
                progressCall(allcnt)

            figclstc = matplotlib.figure.Figure()
            canvas = matplotlib.backends.backend_pdf.FigureCanvasPdf(figclstc) # NOTE this needs to be called to set the the canvas of figclstc
            status = viabilityToMatplotlib(clsplate.clsReplicateGroups[clstcidx],figclstc,showTitle=showTitle,addWellIdsToTitle=addWellIdsToTitle)
            figclstc.savefig(pdfp,format="pdf",bbox_inches='tight')
            pdftitle+=clsplate.clsReplicateGroups[clstcidx].sampleid+" "+clsplate.clsReplicateGroups[clstcidx].condition    

        # set some metadata
        d = pdfp.infodict()
        d['Title'] = 'Viabilities'
        if len(elementaryIndices)+len(replicateGroupIndices) == 1:
            # only one plot, put sample's full id into pdf title
            d['Title'] = pdftitle
        if creator is not None:
            d['Creator'] = creator

def survivalIntegralsToPdf(cls,pdfout,progressCall=None,conditions=None,samples=None,sampleIdToLabel=None,conditionToLabel=None):
    """
    Create barplots of survival integrals (pdf).

    :param pdfout: Filename of pdf file.
    :type pdfout: str
    :param progressCall: Function that will be called on each iteration.
    :type progressCall: @fun(int)
    :param conditions: List of conditions that shall be plotted.
    :type conditions: list(str)
    :param samples: List of samples that shall be plotted.
    :type sample: list(str)
    :param sampleIdToLabel: Dictionary of sample labels.
    :type sample: dict()
    :param conditionToLabel: Dictionary of condition labels.
    :type conditionToLabel: dict()
    """
    with contextlib.closing(matplotlib.backends.backend_pdf.PdfPages(pdfout)) as pdfp:
        fig = matplotlib.figure.Figure()
        canvas = matplotlib.backends.backend_pdf.FigureCanvasPdf(fig)
        fig.set_size_inches(18.5,10.5)
    
        survivalIntegralsToFig(cls,fig,progressCall=progressCall,sampleIdToLabel=sampleIdToLabel,conditionToLabel=conditionToLabel)
    
        fig.savefig(pdfp,format="pdf",bbox_inches='tight')

def survivalIntegralsToFig(cls,fig,progressCall=None,conditions=None,samples=None,sampleIdToLabel=None,conditionToLabel=None):
    """
    Create barplots of survival integrals (matplotlib figure).

    :param fig: matplotlib figure.
    :param progressCall: Function that will be called on each iteration.
    :type progressCall: @fun(int)
    :param conditions: List of conditions that shall be plotted.
    :type conditions: list(str)
    :param samples: List of samples that shall be plotted.
    :type sample: list(str)
    :param sampleIdToLabel: Dictionary of sample labels.
    :type sample: dict()
    :param conditionToLabel: Dictionary of condition labels.
    :type conditionToLabel: dict()
    """
    if conditions is None:
        conditions=cls.conditions()

    survival={}
    clsidx=0
    for co in conditions:
        ctcs=cls.clsReplicateGroupIdcsForCondition(co)
        survival[co]=[]
        # gather survival integrals for this condition
        for clstc in ctcs:
            if progressCall is not None:
                progressCall(clsidx)
            clsidx+=1
            if re.search('blank',clstc.sampleid,re.IGNORECASE) or re.search('background',clstc.sampleid,re.IGNORECASE):
                continue
            si, sivar, status = clstc.survivalIntegral()
            survival[co].append({
                    'si': si,
                    'sierr': math.sqrt(sivar) if sivar is not None else numpy.nan,
                    'condition': co,
                    'sampleid': clstc.sampleid,
                    })

    _survivalIntegralsToFig(survival,cls,fig,progressCall,conditions,samples,sampleIdToLabel,conditionToLabel)

    return survival

def _survivalIntegralsToFig(survival,cls,fig,progressCall=None,conditions=None,samples=None,sampleIdToLabel=None,conditionToLabel=None):
    """
    Helper function for survivalIntegralsToFig.

    For internal use only.
    """
    if samples is not None:
        newlists={}
        for co in conditions:
            newlists[co]=[]
            for clsdict in survival[co]:
                if clsdict['sampleid'] in samples:
                    newlists[co].append(clsdict)
                else:
                    print(clsdict['sampleid']+' '+co+' not in') 
        survival=newlists

    fontsize='large'
    titlefontsize='x-large'

    coloridx=0
    colors=twentysixColours()
    colordict={}
    for co in conditions:
        # sort by survival
        survival[co] = sorted(survival[co], key=lambda k: k['si'])
        # assign colors to sorted survivals (but chooses from colordict if the sample was already assigned a colour)
        for clsdict in survival[co]:        
            if clsdict['sampleid'] not in colordict:
                colordict[clsdict['sampleid']]=colors[coloridx]
                coloridx+=1
                if coloridx >= len(colors):
                    coloridx=0
            clsdict['color']=colordict[clsdict['sampleid']]

    condidx=1
    for condition in conditions:
        sortedbysurvival = survival[condition]
        sis=numpy.array([k['si'] for k in sortedbysurvival],dtype=float)
        sierrs=numpy.array([k['sierr'] for k in sortedbysurvival],dtype=float)
        names=[sampleIdToLabel[k['sampleid']] if sampleIdToLabel is not None and k['sampleid'] in sampleIdToLabel else k['sampleid']
               for k in sortedbysurvival]
        clrs=[k['color'] for k in sortedbysurvival]
        ax=fig.add_subplot(len(list(survival.keys())),1,condidx)
        ax.set_ylabel('Survival Integral',fontsize=fontsize)
        ind = numpy.arange(len(sis))  # the x locations for the groups
        width = 0.75           # the width of the bars
        ax.set_title(conditionToLabel[condition] if conditionToLabel is not None and condition in conditionToLabel else condition,
                     x=.5, y=.8, fontsize=titlefontsize)#ax1ymax*.8)
        ax.set_xlabel('Strains',fontsize=fontsize)

        # remove top and left axes
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.get_yaxis().tick_left()
        ax.get_xaxis().tick_bottom()

        # set names of samples
        ax.set_xticks(ind+width/2)
        ax.set_xticklabels(names, rotation='vertical',fontsize=fontsize)
        # adjsut fontsize of yticks
        # FIXME does not work: ax.set_yticklabels(ax.get_yticklabels(),fontsize=fontsize)
        # the final plot
        rects1 = ax.bar(ind, sis, width, color=clrs, yerr=sierrs)
        condidx+=1
