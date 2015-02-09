"""
This module implements the plot functions for GATHODE.

Growth Analysis Tool for High-throughput Optical Density Experiments
(GATHODE) plot functions.
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
import textwrap

import contextlib
import matplotlib
import matplotlib.backends.backend_pdf
import matplotlib.legend
import matplotlib.legend_handler

from platereader.statusmessage import StatusMessage, Severity

from platereader.numpytools import notNanAndGreaterEqual, notNanAndLess, nonNanSqrt, nonNanNonZeroDivide

def twentysixColours():
    """
    Return the 2010 Colour Alphabet Project.

    :return: list(str) -- List of RGB-encoded colors.
    """
    return [
        #'#FFFFFF', # (white background assumed)
        '#F0A3FF', # Amethyst
        '#0075DC', # Blue
        '#993F00', # Caramel
        '#4C005C', # Damson
        '#191919', # Ebony
        '#005C31', # Forest
        '#2BCE48', # Green
        '#FFCC99', # Honeydew
        '#808080', # Iron
        '#94FFB5', # Jade
        '#8F7C00', # Khaki
        '#9DCC00', # Lime
        '#C20088', # Mallow
        '#003380', # Navy
        '#FFA405', # Orpiment
        '#FFA8BB', # Pink
        '#426600', # Quagmire
        '#FF0010', # Red
        '#5EF1F2', # Sky
        '#00998F', # Turquoise
        '#E0FF66', # Uranium
        '#740AFF', # Violet
        '#990000', # Wine
        '#FFFF80', # Xanthin
        '#FFFF00', # Yellow
        '#FF5005', # Zinnia
        ]

def darkSubsetOfTwentysixColours():
    """
    Return only the dark colors of the 2010 Colour Alphabet Project.

    :return: list(str) -- List of RGB-encoded colors.
    """
    return [
        #'#FFFFFF', # (white background assumed)
        # too bright '#F0A3FF', # Amethyst
        '#0075DC', # Blue
        '#993F00', # Caramel
        '#4C005C', # Damson
        '#191919', # Ebony
        '#005C31', # Forest
        '#2BCE48', # Green
        # too bright '#FFCC99', # Honeydew
        '#808080', # Iron
        # too bright '#94FFB5', # Jade
        '#8F7C00', # Khaki
        '#9DCC00', # Lime
        '#C20088', # Mallow
        '#003380', # Navy
        '#FFA405', # Orpiment
        # too bright '#FFA8BB', # Pink
        '#426600', # Quagmire
        '#FF0010', # Red
        # too bright '#5EF1F2', # Sky
        '#00998F', # Turquoise
        # too bright '#E0FF66', # Uranium
        '#740AFF', # Violet
        '#990000', # Wine
        # too bright '#FFFF80', # Xanthin
        # too bright '#FFFF00', # Yellow
        '#FF5005', # Zinnia
        ]

def colorsForSampleCondition(plate):
    """
    Find a color for each sample-condition pair.

    Return a dictionary with a structure like this:
    { 'sample1': { 'condition1': '#colHex1', 'condition2': '#colHex2', ...}, 'sample2': {...}, ...}

    :return: dict() -- color for sampleid-condition
    """
    colors=darkSubsetOfTwentysixColours()
    coloridx=0
    colorForSampleCondition={}

    # the background for all conditions are colored black
    for tc in plate.backgroundWells():
        if tc.sampleid not in colorForSampleCondition:
            colorForSampleCondition[tc.sampleid]={}
        colorForSampleCondition[tc.sampleid][tc.condition]='black'
    for tc in plate.backgroundReplicateGroups():
        if tc.sampleid not in colorForSampleCondition:
            colorForSampleCondition[tc.sampleid]={}
        colorForSampleCondition[tc.sampleid][tc.condition]='black'

    # all the others are distinguished by
    for tc in plate.wells:
        if tc.sampleid not in colorForSampleCondition:
            colorForSampleCondition[tc.sampleid]={}
        if tc.condition not in colorForSampleCondition[tc.sampleid]:
            colorForSampleCondition[tc.sampleid][tc.condition]=colors[coloridx]
            coloridx += 1
            if coloridx == len(colors):
                coloridx = 0

    return colorForSampleCondition

# NOTE if you want to feed colors you might want to have a look at http://stackoverflow.com/a/4382138
# see also http://graphicdesign.stackexchange.com/questions/3682/large-color-set-for-coloring-of-many-datasets-on-a-plot
def replicateToFig(fig,replicate,addToTitle=None,
                   showTitle=True,
                   showDerivatives=None,
                   showRaw=True,showBackground=True,showSingle=True,showSmoothed=True,
                   showMaxLinearSlope=False,
                   showLogod=True,showLogodSmoothed=False,
                   showMaxGrowthrate=True, showMaxGrowthrateFromLogOdDerivative=True,
                   showDerivativeLinear=True, showSmoothedDerivativeLinear=True,
                   showLogOdDerivative=True,showLogOdDerivativeFromNonLog=True,showLogOdDerivativeFromNonLogSmoothed=True,
                   showExpFitsOd0Mu=True,showExpFitsMu=False,
                   showGrowthyield=True,
                   colors=None,
                   legendkwargs=None,
                   derlegendkwargs=None,
                   ax_tick_params=None,
                   ax2_tick_params=None,
                   derax_tick_params=None,
                   derax2_tick_params=None):
    """
    Show features (such as OD, log(OD), derivatives, maximal growth, ...) on a matplotlib figure.

    :param show*: Whether to show a certain feature of this Replicate.
    :type show*: bool
    :param showDerivatives: if None, show if derivatives defined; if False don't show them (overrides show* options)
    :type showDerivatives: bool

    :param colors: override default colors
    :type colors: dict
    :param legendkwargs: kwargs passed to legend; when None, defaults to {'loc': 0, 'prop': {'size': 8}}
    :type legendkwargs: dict
    :param ax_tick_params: kwargs passed to tick_params for first axis
    :type ax_tick_params: dict
    :param ax2_tick_params: kwargs passed to tick_params for second axis; when None, defaults to ax_tick_params
    :type ax2_tick_params: dict
    :param derax_tick_params: kwargs passed to tick_params for first axis of derivative plot; when None, defaults to ax_tick_params
    :type derax_tick_params: dict
    :param derax2_tick_params: kwargs passed to tick_params for second axis of derivative plot; when None, defaults to derax_tick_params
    :type derax2_tick_params: dict

    FIXME enhance documentation of parameters.

    :return: StatusMessage -- Statuses of this dataset.
    """
    statuslist=StatusMessage()

    if showDerivatives is None:
        # not set, if derivatives exist show them
        showDerivatives=replicate.derivative() is not None or replicate.smoothedOdDerivative() is not None

    if derax_tick_params is None and ax_tick_params is not None:
        derax_tick_params=dict(ax_tick_params)

    ## title
    title=None
    if showTitle:
        title=replicate.sampleid+" "+replicate.condition
        if addToTitle is not None:
            title+=" "+addToTitle
        # add linebreak to title if too long
        title="\n".join(textwrap.wrap(title, 100))

    if showDerivatives:
        ax=fig.add_subplot(2,1,1)
    else:
        ax=fig.add_subplot(1,1,1)

    dm, xmin, xmax= dataToMatplotlibAxes(replicate,ax,title=title,
                                         showRaw=showRaw,showBackground=showBackground,showSingle=showSingle,
                                         showSmoothed=showSmoothed,showGrowthyield=showGrowthyield,
                                         showMaxLinearSlope=showMaxLinearSlope,
                                         showLogod=showLogod,showLogodSmoothed=showLogodSmoothed,
                                         showMaxGrowthrate=showMaxGrowthrate, showMaxGrowthrateFromLogOdDerivative=showMaxGrowthrateFromLogOdDerivative,
                                         showExpFitsOd0Mu=showExpFitsOd0Mu,showExpFitsMu=showExpFitsMu,
                                         statuslist=statuslist,dontShowXlabel=showDerivatives,
                                         colors=colors,
                                         legendkwargs=legendkwargs,
                                         ax_tick_params=ax_tick_params,
                                         ax2_tick_params=ax2_tick_params)

    if showDerivatives:
        axder=fig.add_subplot(2,1,2)
        derivativesToMatplotlibAxes(replicate,axder,
                                    showMaxLinearSlope=showMaxLinearSlope,
                                    showMaxGrowthrate=showMaxGrowthrate, showMaxGrowthrateFromLogOdDerivative=showMaxGrowthrateFromLogOdDerivative,
                                    showDerivativeLinear=showDerivativeLinear, showSmoothedDerivativeLinear=showSmoothedDerivativeLinear,
                                    showLogOdDerivative=showLogOdDerivative,showLogOdDerivativeFromNonLog=showLogOdDerivativeFromNonLog,
                                    showLogOdDerivativeFromNonLogSmoothed=showLogOdDerivativeFromNonLogSmoothed,
                                    showExpFitsOd0Mu=showExpFitsOd0Mu,showExpFitsMu=showExpFitsMu,
                                    xmin=xmin,xmax=xmax,
                                    statuslist=statuslist,
                                    colors=colors,
                                    derlegendkwargs=derlegendkwargs,
                                    derax_tick_params=derax_tick_params,
                                    derax2_tick_params=derax2_tick_params)

    return statuslist

def plotcolors(colors=None):
    """
    Determine colors, based on given colors and on defaults.

    :return: dict() -- dictionary of colors (values) for features (keys).
    """
    defaultcolors={
        'od': 'blue',
        'singleod': 'lightgreen',
        'odsmoothed': 'red',
        'oderr': 'blue',
        'oderralpha': 0.2,
        'logod': 'brown',
        'logoderr': 'brown',
        'logoderralpha': 0.2,
        'logodsmoothed': 'black',
        'logodderivative': 'green',
        'rawod': 'black',
        'rawoderr': 'gray',
        'singlerawod': 'lightgreen',
        'rawoderralpha': 0.3,
        'background': 'yellow',
        'backgrounderr': 'yellow',
        'backgrounderralpha': 0.4,
        'smoothedod': 'red',
        'linearmaxslope': 'cyan',
        'growthyield': 'yellow',
        'growthyieldedge': 'black',
        'expfitsod0mu': 'red',
        'expfitsod0muerr': 'red',
        'expfitsod0muerralpha': 0.2,
        'expfitsmu': 'cyan',
        'maxgrowthrate': 'red',
        'maxgrowthratefromlogod': 'black',
        }
    colors_={}
    for key in defaultcolors:
        if colors is not None and key in colors:
            colors_[key]=colors[key]
        else:
            colors_[key]=defaultcolors[key]
    return colors_

def dataToMatplotlibAxes(replicate,ax,title=None,
                         showRaw=True,showBackground=True,showSingle=True,showSmoothed=True,
                         showOd=True,
                         showMaxLinearSlope=False,
                         showLogod=True,showLogodSmoothed=False,
                         showMaxGrowthrate=True, showMaxGrowthrateFromLogOdDerivative=True,
                         showExpFitsOd0Mu=True,showExpFitsMu=False,
                         showGrowthyield=True,
                         statuslist=None,colors=None,dontShowXlabel=False,
                         legendkwargs=None,
                         ax_tick_params=None,
                         ax2_tick_params=None):
    """
    Show features not associated with derivative (such as OD, log(OD), maximal growth, ...) on a matplotlib figure.

    :param ax: The matplotlib Axes object used for plotting.
    :type ax: matplotlib.axis.Axis

    :param show*: Whether to show a certain feature of this Replicate.
    :type show*: bool

    :param colors: override default colors
    :type colors: dict
    :param legendkwargs: kwargs passed to legend; when None, defaults to {'loc': 0, 'prop': {'size': 8}}
    :type legendkwargs: dict
    :param ax_tick_params: kwargs passed to tick_params for first axis
    :type ax_tick_params: dict
    :param ax2_tick_params: kwargs passed to tick_params for second axis; when None, defaults to ax_tick_params
    :type ax2_tick_params: dict

    FIXME enhance documentation of parameters.

    :return: StatusMessage -- Statuses of this dataset.
    """
    if legendkwargs is None:
        legendkwargs={
            'loc': 0,
            'prop': {'size': 8},
            'numpoints': 1,
            'scatterpoints': 1,
            #'markerscale': 0.8,
            #'handlelength': 3,
            }
    # make sure that error bars are not too close to symbol in legend
    matplotlib.legend.Legend.update_default_handler_map({
            matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0.8,yerr_size=0.8)
            })
    if ax_tick_params is None:
        ax_tick_params={}
    if ax2_tick_params is None:
        ax2_tick_params=ax_tick_params
    colors_=plotcolors(colors)

    if statuslist is None:
        statuslist=StatusMessage()

    if replicate.timeunit == "s":
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    if title is not None:
        ax.set_title(title)

    ax2=None

    if len(replicate.activeChildWellIndices()) == 0:
        statuslist.addStatus(StatusMessage(
                key='replicates:',shortmsg='noActiveChild',
                longmsg='none of the children is active',
                severity=Severity.failed))
        return statuslist, None, None
    if replicate.rawOd() is None:
        statuslist.addStatus(StatusMessage(
                key='replicates:',shortmsg='noRawOd',
                longmsg='raw optical density is not defined',
                severity=Severity.failed))
        return statuslist, None, None

    (ax1xmin,ax1xmax,ax1ymin,ax1ymax)=(replicate.time.min(),replicate.time.max(),None,None)
    (ax2xmin,ax2xmax,ax2ymin,ax2ymax)=(ax1xmin,ax1xmax,None,None)

    if replicate.od() is not None and showOd:
        ax.plot(replicate.time,replicate.od(),label="bckgrnd & HD corr.",color=colors_['od'])
        if replicate.odVar() is not None:
            err=numpy.sqrt(replicate.odVar())
            if showSingle:
                # show each underlying replicate
                for stc in replicate.activeChildWells():
                    if stc.od() is not None:
                        ax.plot(stc.time,stc.od(),label="",color=colors_['singleod'])
            ax.fill_between(replicate.time, replicate.od()-err, replicate.od()+err,
                            alpha=colors_['oderralpha'],edgecolor=colors_['oderr'],facecolor=colors_['oderr'])

    if ((showLogod and replicate.logOd() is not None)
        or (showLogodSmoothed and replicate.logOdSmoothed() is not None)
        or (showMaxGrowthrate and replicate.od() is not None)
        or (showMaxGrowthrateFromLogOdDerivative and replicate.od() is not None)):
        ax2 = ax.twinx()
        ax2.set_ylabel("ln(OD)")

    if replicate.logOdCutoff() is not None and replicate.logOd() is not None:
        logltidcs=notNanAndLess(replicate.logOd(),replicate.logOdCutoff())
    else:
        logltidcs=numpy.array([False for i in range(0,len(replicate.time))],dtype=bool)

    if showLogod and replicate.logOd() is not None:
        logodcopy=numpy.copy(replicate.logOd())
        logodcopy[logltidcs] = numpy.nan
        if numpy.any(~numpy.isnan(logodcopy)):
            # plotting does not work if all values are nan; it is not enough to check logltidcs, because
            # logod can also be nan if od is less than zero (or rather: less then 1e-35 cutoff)
            ax2.plot(replicate.time,logodcopy,label="natural log",color=colors_['logod'])
        if numpy.any(~numpy.isnan(logodcopy)) and replicate.odVar() is not None:
            # create errors from odVar
            logoderr=nonNanNonZeroDivide(numpy.sqrt(replicate.odVar()),numpy.absolute(replicate.od()))
            logoderr[numpy.isnan(logodcopy)] = numpy.nan
            ax2.fill_between(replicate.time,
                             numpy.ma.masked_array(logodcopy-logoderr,numpy.isnan(logoderr)),
                             numpy.ma.masked_array(logodcopy+logoderr,numpy.isnan(logoderr)),
                             alpha=colors_['logoderralpha'],edgecolor=colors_['logoderr'],facecolor=colors_['logoderr'])
        (ax2xmin,ax2xmax,ax2ymin,ax2ymax)=ax2.axis()

    if showRaw or replicate.od() is None:
        if replicate.rawOdVar() is not None:
            err=numpy.sqrt(replicate.rawOdVar())
            activeTcs=replicate.activeChildWells()
            if replicate.od() is None and len(activeTcs) > 1:
                # we only show the raw, so we can show all replicates on top of each other
                for stc in activeTcs:
                    ax.plot(stc.time,stc.rawOd(),label="",color=colors_['singlerawod'])
            ax.fill_between(replicate.time, replicate.rawOd()-err, replicate.rawOd()+err,
                            alpha=colors_['rawoderralpha'],edgecolor=colors_['rawoderr'],facecolor=colors_['rawoderr'])
        ax.plot(replicate.time,replicate.rawOd(),label="raw",color=colors_['rawod'])

    if showBackground and replicate.background is not None:
        ax.plot(replicate.time,replicate.background.rawOd(),label="background",color=colors_['background'])
        if replicate.background.rawOdVar() is not None:
            err=numpy.sqrt(replicate.background.rawOdVar())
            ax.fill_between(replicate.time, replicate.background.rawOd()-err, replicate.background.rawOd()+err,
                            alpha=colors_['backgrounderralpha'], edgecolor=colors_['backgrounderr'], facecolor=colors_['backgrounderr'])
    if showSmoothed and replicate.smoothedOd() is not None:
        ax.plot(replicate.time,replicate.smoothedOd(),color=colors_['smoothedod'],
                label="smoothed""\n(s="+str(replicate.smoothingS())+", k="+str(replicate.smoothingK())+")")
        if replicate.smoothedOdDerivative() is not None and showMaxLinearSlope:
            slopemax, slopemaxVar, interceptmax, interceptmaxVar, timemax, timemaxVar, timemaxIndices, plainSlopeStatus=replicate.odSlopemaxIntercept()
            if plainSlopeStatus is not None:
                statuslist.addStatus(plainSlopeStatus)
            if slopemax is not None:
                flin=lambda t: slopemax*t+interceptmax
                (dummyax1xmin,dummyax1xmax,ax1ymin,ax1ymax)=ax.axis()
                ax.plot(replicate.time, list(map(flin,replicate.time)),label="linear",color=colors_['linearmaxslope'])
                ax.plot([timemax],[flin(timemax)],marker='o',color=colors_['linearmaxslope'])

    if showGrowthyield and replicate.od() is not None:
        growthyield, growthyieldvar, tgrowthyield, tgrowthyieldvar, status=replicate.growthyield()
        if status is not None:
            statuslist.addStatus(status)
        if growthyield is not None:
            ax.scatter([tgrowthyield],[growthyield],s=40,alpha=1.,zorder=42,
                       marker='o',color=colors_['growthyield'],
                       edgecolor=colors_['growthyieldedge'],label='yield',
                       linewidths=.8,
                       )

    if showLogodSmoothed and replicate.logOdSmoothed() is not None:
        logodsmoothedcopy=numpy.copy(replicate.logOdSmoothed())
        logodsmoothedcopy[logltidcs] = numpy.nan
        if numpy.any(~numpy.isnan(logodsmoothedcopy)):
            ax2.plot(replicate.time,logodsmoothedcopy,color=colors_['logodsmoothed'],
                     label="smoothed""\n(s="+str(replicate.smoothingS())+", k="+str(replicate.smoothingK())+")")
            if ax2xmax is not None:
                ax2.axis((ax2xmin,ax2xmax,ax2ymin,ax2ymax))
            (ax2xmin,ax2xmax,ax2ymin,ax2ymax)=ax2.axis()

    if showMaxGrowthrate and replicate.od() is not None:
        mu, muvar, od0, od0var, maxt, maxtvar, lag, lagvar, method, expfitExpmaxStatus = replicate.maxGrowthrate()
        if mu is not None and od0 > 0:
            lagAtLogOdEquals=replicate.lagAtLogOdEquals()
            flin=lambda t: mu*t+math.log(od0)
            ax2.plot(replicate.time, list(map(flin,replicate.time)),color=colors_['maxgrowthrate'])
            y2mumax=mu*maxt+math.log(od0)
            ax2.plot([maxt],[y2mumax],label='max. $\mu$',marker='o',color=colors_['maxgrowthrate'])
            if lag is not None:
                ax2.errorbar([lag],[lagAtLogOdEquals],
                             xerr=[math.sqrt(lagvar)] if lagvar is not None else None,
                             label='lag',
                             marker='x',markeredgewidth=2,color=colors_['maxgrowthrate'])
            if ax2xmin is not None:
                ax2.axis((ax2xmin,ax2xmax,ax2ymin,ax2ymax))
            (ax2xmin,ax2xmax,ax2ymin,ax2ymax)=ax2.axis()
            if y2mumax-math.fabs(y2mumax)*.05 < ax2ymin:
                ax2ymin=y2mumax-math.fabs(y2mumax)*.05
            if lag is not None and ax2ymax >= lagAtLogOdEquals-math.fabs(lagAtLogOdEquals)*.05:
                ax2ymin=lagAtLogOdEquals-math.fabs(lagAtLogOdEquals)*.05
            ax2.axis((ax2xmin,ax2xmax,ax2ymin,ax2ymax))

    if showMaxGrowthrateFromLogOdDerivative and replicate.od() is not None:
        mu, muvar, od0, od0var, maxt, maxtvar, lag, lagvar, method, nonlogExpmaxStatus = replicate.maxGrowthrateFromLogOdDerivative()
        if mu is not None and od0 > 0:
            lagAtLogOdEquals=replicate.lagAtLogOdEquals()
            flin=lambda t: mu*t+math.log(od0)
            ax2.plot(replicate.time, list(map(flin,replicate.time)),color=colors_['maxgrowthratefromlogod'])
            y2mumax=mu*maxt+math.log(od0)
            ax2.plot([maxt],[y2mumax],label='max. $\mu$ (1/OD(t)*dOD(t)/dt smoothed)',marker='o',color=colors_['maxgrowthratefromlogod'])
            if lag is not None:
                ax2.errorbar([lag],[lagAtLogOdEquals],
                             xerr=[math.sqrt(lagvar)] if lagvar is not None else None,
                             label='lag (1/OD(t)*dOD(t)/dt smoothed)',
                             marker='x',markeredgewidth=2,color=colors_['maxgrowthratefromlogod'])
            if ax2xmin is not None:
                ax2.axis((ax2xmin,ax2xmax,ax2ymin,ax2ymax))
            (ax2xmin,ax2xmax,ax2ymin,ax2ymax)=ax2.axis()
            if y2mumax-math.fabs(y2mumax)*.05 < ax2ymin:
                ax2ymin=y2mumax-math.fabs(y2mumax)*.05
            if lag is not None and ax2ymax >= lagAtLogOdEquals-math.fabs(lagAtLogOdEquals)*.05:
                ax2ymin=lagAtLogOdEquals-math.fabs(lagAtLogOdEquals)*.05
            ax2.axis((ax2xmin,ax2xmax,ax2ymin,ax2ymax))

    if ax1xmin is not None:
        ax.axis((ax1xmin,ax1xmax,ax1ymin,ax1ymax))

    (sub1_ax1xmin,sub1_ax1xmax,sub1_ax1ymin,sub1_ax1ymax)=ax.axis()

    if ax2 is not None:
        # for legend see also http://stackoverflow.com/questions/5484922/secondary-axis-with-twinx-how-to-add-to-legend
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        maxlen=0
        for lab in labels + labels2:
            if len(lab) > maxlen:
                maxlen=len(lab)
        dummy = matplotlib.legend.Rectangle((0, 0.5), .01, .01, fc="black", alpha=.0,)
        ax2.legend(lines + [dummy] + lines2, labels + ['-' * maxlen] + labels2, **legendkwargs)
    else:
        ax.legend(**legendkwargs)

    if not dontShowXlabel:
        ax.set_xlabel("time/"+replicate.timeunit)

    ax.set_ylabel("OD")
#    ax.xaxis.set_major_locator(matplotlib.axes.mticker.MaxNLocator(10))
    ax.xaxis.set_minor_locator(matplotlib.axes.mticker.AutoMinorLocator(10))
    ax.axhline(0,ls='--',color='black')

    ax.tick_params(**ax_tick_params)
    if ax2 is not None:
        ax2.tick_params(**ax2_tick_params)
    return statuslist, sub1_ax1xmin, sub1_ax1xmax

def derivativesToMatplotlibAxes(replicate,ax,
                                showMaxLinearSlope=False,
                                showMaxGrowthrate=True, showMaxGrowthrateFromLogOdDerivative=True,
                                showDerivativeLinear=True, showSmoothedDerivativeLinear=True,
                                showLogOdDerivative=True,showLogOdDerivativeFromNonLog=True,showLogOdDerivativeFromNonLogSmoothed=True,
                                showExpFitsOd0Mu=True,showExpFitsMu=False,
                                xmin=None,xmax=None,
                                statuslist=None,colors=None,
                                derlegendkwargs=None,
                                derax_tick_params=None,
                                derax2_tick_params=None):
    """
    Show features associated with derivative (dOD/dt, maximal growth, ...) on a matplotlib figure.

    :param ax: The matplotlib Axes object used for plotting.
    :type ax: matplotlib.axis.Axis

    :param show*: Whether to show a certain feature of this Replicate.
    :type show*: bool

    :param colors: override default colors
    :type colors: dict
    :param legendkwargs: kwargs passed to legend; when None, defaults to {'loc': 0, 'prop': {'size': 8}}
    :type legendkwargs: dict
    :param derax_tick_params: kwargs passed to tick_params for first axis
    :type derax_tick_params: dict
    :param derax2_tick_params: kwargs passed to tick_params for second axis; when None, defaults to derax_tick_params
    :type derax2_tick_params: dict

    FIXME enhance documentation of parameters.

    :return: StatusMessage -- Statuses of this dataset.
    """
    if derlegendkwargs is None:
        derlegendkwargs={
            'loc': 0,
            'prop': {'size': 8},
            'numpoints': 1, 
            }
    # make sure that error bars are not too close to symbol in legend
    matplotlib.legend.Legend.update_default_handler_map({
            matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0.8,yerr_size=0.8)
            })
    if derax_tick_params is None:
        derax_tick_params={}
    if derax2_tick_params is None:
        derax2_tick_params=derax_tick_params
    colors_=plotcolors(colors)

    if statuslist is None:
        statuslist=StatusMessage()
    expfitExpmaxStatus=None
    nonlogExpmaxStatus=None

    (ax1xmin,ax1xmax,ax1ymin,ax1ymax)=(None,None,None,None)
    (ax2xmin,ax2xmax,ax2ymin,ax2ymax)=(None,None,None,None)

    if replicate.timeunit == "s":
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    if showDerivativeLinear and replicate.derivative() is not None:
        ax.plot(replicate.time[:-1],replicate.derivative(),label="bckgrnd & HD corr.",color=colors_['od'])
    if showSmoothedDerivativeLinear and replicate.smoothedOdDerivative() is not None:
        ax.plot(replicate.time[:-1],replicate.smoothedOdDerivative(),label="smoothed",color=colors_['odsmoothed'])
        if showMaxLinearSlope:
            slopemax, slopemaxVar, interceptmax, interceptmaxVar, timemax, timemaxVar, timemaxIndices, plainSlopeStatus=replicate.odSlopemaxIntercept()
            if slopemax is not None:
                ax.plot([timemax],[slopemax],marker='o',color=colors_['linearmaxslope'])

    ax.set_xlabel("time/"+replicate.timeunit)
    ax.set_ylabel("dOD(t)/dt")
    ax.axhline(0,ls='--',color='black')

    ax2=None
    if (showLogOdDerivative or showExpFitsOd0Mu or showExpFitsMu
        or showLogOdDerivativeFromNonLog
        or showLogOdDerivativeFromNonLogSmoothed
        or showMaxGrowthrate
        or showMaxGrowthrateFromLogOdDerivative):
        ax2 = ax.twinx()
    if replicate.logOdCutoff() is not None and replicate.logOd() is not None:
        tmb=int(math.floor(replicate.slidingWindowSize()/2.))
        tmt=int(math.ceil(replicate.slidingWindowSize()/2.))
        logltidcsMinWin=notNanAndLess(replicate.logOd()[tmb:-tmt],replicate.logOdCutoff())
        logltidcsMin1=notNanAndLess(replicate.logOd()[:-1],replicate.logOdCutoff())
    else:
        logltidcsMinWin=numpy.array([False for i in range(0,len(replicate.time)-replicate.slidingWindowSize())],dtype=bool)
        logltidcsMin1=numpy.array([False for i in range(0,len(replicate.time)-1)],dtype=bool)

    if showLogOdDerivative:
        logodderivative=numpy.copy(replicate.logOdDerivative())
        logodderivative[logltidcsMin1] = numpy.nan
        if numpy.any(~numpy.isnan(logodderivative)):
            ax2.plot(replicate.time[:-1],logodderivative,label="local der.",color=colors_['logodderivative'])

    tslidwin = replicate.time[int(math.floor(replicate.slidingWindowSize()/2.)):int(-math.ceil(replicate.slidingWindowSize()/2.))]
    if showExpFitsOd0Mu:
        mu_withOd0, muvar_withOd0, od0, od0var = replicate.expFitsOd0Mu()
        if mu_withOd0 is not None and len(mu_withOd0):
            muCopy_withOd0=numpy.copy(mu_withOd0)
            muCopy_withOd0[logltidcsMinWin] = numpy.nan
            if numpy.any(~numpy.isnan(muCopy_withOd0)):
                ax2.plot(tslidwin,muCopy_withOd0,label='$\mu$ (exp. fits)',color=colors_['expfitsod0mu'])

    if showExpFitsMu:
        mu_onlyMu, muvar_onlyMu = replicate.expFitsMu()
        if mu_onlyMu is not None and len(mu_onlyMu):
            muCopy_onlyMu=numpy.copy(mu_onlyMu)
            muCopy_onlyMu[logltidcsMinWin] = numpy.nan
            if numpy.any(~numpy.isnan(muCopy_onlyMu)):
                ax2.plot(tslidwin,muCopy_onlyMu,label='$\mu$ (exp. fits, fixed OD(0))',color=colors_['expfitsmu'])

    if showLogOdDerivativeFromNonLogSmoothed:
        logodderivativefnlsm=replicate.logOdDerivativeFromNonLogSmoothed()
        if logodderivativefnlsm is not None and len(logodderivativefnlsm):
            logodderivativefnlsmCopy=numpy.copy(logodderivativefnlsm)
            logodderivativefnlsm[logltidcsMin1] = numpy.nan
            if numpy.any(~numpy.isnan(logodderivativefnlsm)):
                ax2.plot(replicate.time[:-1], logodderivativefnlsm,
                         label="1/OD(t)*dOD(t)/dt smoothed",color=colors_['logodsmoothed'])

    if ax2 is not None:
        (ax2xmin,ax2xmax,ax2ymin,ax2ymax)=ax2.axis()

    if showExpFitsOd0Mu:
        # this shall not influence the y-limits, therefore it is not plotted before
        if mu_withOd0 is not None and len(mu_withOd0) and muvar_withOd0 is not None and len(muvar_withOd0):
            muerr_withOd0 = nonNanSqrt(muvar_withOd0)
            muerr_withOd0[logltidcsMinWin] = numpy.nan
            if numpy.any(~numpy.isnan(muerr_withOd0)):
                ax2.fill_between(tslidwin, 
                                 numpy.ma.masked_array(muCopy_withOd0-muerr_withOd0,numpy.isnan(muerr_withOd0)),
                                 numpy.ma.masked_array(muCopy_withOd0+muerr_withOd0,numpy.isnan(muerr_withOd0)),
                                 alpha=colors_['expfitsod0muerralpha'],edgecolor=colors_['expfitsod0muerr'],facecolor=colors_['expfitsod0muerr'])

    # we calculate these no matter whether they are shown or not, in order to set the ax2ymax limit (see below)
    mu_nl, muvar_nl, od0_nl, od0var_nl, maxt_nl, maxtvar_nl, lag_nl, lagvar_nl, method_nl, nonlogExpmaxStatus = replicate.maxGrowthrateFromLogOdDerivative()
    mu_ef, muvar_ef, od0_ef, od0var_ef, maxt_ef, maxtvar_ef, lag_ef, lagvar_ef, method_ef, expfitExpmaxStatus = replicate.maxGrowthrate()
    if showMaxGrowthrate and replicate.od() is not None and mu_ef is not None:
        ax2.errorbar([maxt_ef],[mu_ef],
                     xerr=[math.sqrt(maxtvar_ef)] if maxtvar_ef is not None else None,
                     yerr=[math.sqrt(muvar_ef)] if muvar_ef is not None else None,
                     marker='o',label='max. $\mu$',color=colors_['maxgrowthrate'])

    if showMaxGrowthrateFromLogOdDerivative and replicate.od() is not None and mu_nl is not None:
        ax2.errorbar([maxt_nl],[mu_nl],
                     xerr=[math.sqrt(maxtvar_nl)] if maxtvar_nl is not None else None,
                     yerr=[math.sqrt(muvar_nl)] if muvar_nl is not None else None,
                     marker='o',label='max. $\mu$ (1/OD(t)*dOD(t)/dt smoothed)',color=colors_['maxgrowthratefromlogod'])

    if showLogOdDerivativeFromNonLog:
        logodderivativefnl=replicate.logOdDerivativeFromNonLog()
        if logodderivativefnl is not None and len(logodderivativefnl):
            logodderivativefnlCopy=numpy.copy(logodderivativefnl)
            logodderivativefnlCopy[logltidcsMin1] = numpy.nan
            if numpy.any(~numpy.isnan(logodderivativefnlCopy)):
                ax2.plot(replicate.time[:-1], logodderivativefnlCopy,
                         label="1/OD(t)*dOD(t)/dt",color=colors_['logod'])

    if ax2 is not None:
        if ax2xmax is None:
            (ax2xmin,ax2xmax,ax2ymin,ax2ymax)=ax2.axis()
        if ax2ymin < 0 and ax2ymax > 0:
            # fix insanely low ymin
            ax2ymin=-ax2ymax/10

        # of the two mu_max, see which one would lead to a higher ymax of ax2
        mupluserr=None
        if mu_nl is not None:
            mupluserr=mu_nl+math.sqrt(muvar_nl) if muvar_nl is not None else mu_nl
        if mu_ef is not None:
            mupluserr_ef=mu_ef+math.sqrt(muvar_ef) if muvar_ef is not None else mu_ef
            if mupluserr is None or mupluserr_ef > mupluserr:
                mupluserr=mupluserr_ef
        # if the mu + error does not lead to an "unhealthy" increase of ymax, use this as ymax
        if mupluserr is not None and ax2ymax is not None and ax2ymax < mupluserr:
            ax2ymax = min((ax2ymax-ax2ymin)*1.5+ax2ymin, 1.08*mupluserr)

        ax2.axis((ax2xmin,ax2xmax,ax2ymin,ax2ymax))

    if xmin is None:
        xmin=numpy.min(replicate.time)
    if xmax is None:
        xmax=numpy.max(replicate.time)
    ax.axis(xmin=xmin, xmax=xmax)
    if ax2 is not None:
        ax2.set_ylabel("growth rate")
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        maxlen=0
        for lab in labels + labels2:
            if len(lab) > maxlen:
                maxlen=len(lab)
        dummy = matplotlib.legend.Rectangle((0, 0.5), .01, .01, fc="black", alpha=.0,)
        ax2.axis(xmin=xmin, xmax=xmax)
        ax2.legend(lines + [dummy] + lines2, labels + ['-' * maxlen] + labels2, **derlegendkwargs)
        ax2.xaxis.set_minor_locator(matplotlib.axes.mticker.AutoMinorLocator(10))

    else:
        ax.legend(**derlegendkwargs)

    if showMaxGrowthrateFromLogOdDerivative and nonlogExpmaxStatus is not None:
        statuslist.addStatus(nonlogExpmaxStatus)
    if showMaxGrowthrate and expfitExpmaxStatus is not None:
        statuslist.addStatus(expfitExpmaxStatus)

    ax.tick_params(**derax_tick_params)
    if ax2 is not None:
        ax2.tick_params(**derax2_tick_params)
    return statuslist

def odPlateOverviewToAxes(ax,plate):
    """
    Show a scheme of the plate with small growth curves inside.

    :param ax: The matplotlib Axes object used for plotting.
    :type ax: matplotlib.axis.Axis

    :param plate: The plate
    :type plate: platereader.plate
    """
    # NOTE matplotlib is too slow to generate 384 subplots, so we
    # generate the subplots ourselves by translating the arrays.

    legendprop={'size': 0}

    showPlateLabels=True
    if len(plate.wells) == 384:
        plateformat='384'
        cols=24
        rows=16
    elif len(plate.wells) == 96:
        plateformat='96'
        cols=12
        rows=8
    elif len(plate.wells) == 200:
        plateformat='200honeycomb'
        cols=20
        rows=10
    else:
        cols=int(round(math.sqrt(len(plate.wells))))
        rows=int(math.ceil(len(plate.wells)/float(cols)))
        showPlateLabels=False

    # find global x- and y-ranges
    (minxmin,maxxmax,minymin,maxymax)=(+1e32,-1e32,+1e32,-1e32)
    for tc in plate.wells:
        (xmin,xmax,ymin,ymax)=(tc.time.min(),tc.time.max(),tc.rawOd().min(),tc.rawOd().max())
        if xmin < minxmin:
            minxmin=xmin
        if xmax > maxxmax:
            maxxmax=xmax
        if ymin < minymin:
            minymin=ymin
        if ymax > maxymax:
            maxymax=ymax

    xoff=1.0*maxxmax
    yoff=1.2*maxymax
    xsep=.2*maxxmax
    ysep=.2*maxymax
    colwidth = maxxmax+xsep
    rowheight = maxymax+ysep
    width = colwidth*cols + xoff
    height = rowheight*rows + yoff
    if plateformat == '200honeycomb':
        height+=(maxymax+ysep)/2

    ax.tick_params(bottom=False,top=False,left=False,right=False,labelbottom=False,labeltop=False,labelleft=False,labelright=False)
    ax.axis((0,width,0,height))

    colorForSampleCondition=colorsForSampleCondition(plate)
    rectangles=[]
    tcidx = 0
    for tc in plate.wells:
        offsets = _overviewWellOffsets(tc,tcidx,cols,rows,xoff,yoff,xsep,ysep,colwidth,rowheight,width,height,maxxmax,maxymax,plateformat)

        ax.plot(tc.time + offsets['wellxoff'],
                tc.rawOd() + offsets['wellyoff'],
                color=colorForSampleCondition[tc.sampleid][tc.condition],
                linewidth=2)
        rectangles.append({
                'x1': offsets['x1'],
                'x2': offsets['x2'],
                'y1': offsets['y1'],
                'y2': offsets['y2'],
                'tooltip': tc.sampleid+" "+tc.condition,
                'wellidx': tcidx,
                'well': tc,
                })

        # border around subplot (NOTE for non-honeycomb wells this is done below)
        if plateformat == '200honeycomb':
            ax.axvline(offsets['x1'],ymin=offsets['rely1'],ymax=offsets['rely2'],color='gray')
            ax.axvline(offsets['x2'],ymin=offsets['rely1'],ymax=offsets['rely2'],color='gray')
            ax.axhline(offsets['y1'],xmin=offsets['relx1'],xmax=offsets['relx2'],color='gray')
            ax.axhline(offsets['y2'],xmin=offsets['relx1'],xmax=offsets['relx2'],color='gray')

        tcidx += 1

    if plateformat != '200honeycomb':
        # border around subplots (NOTE for honeycomb wells this is done above)
        xs = 1 - (colwidth*cols + .5*xsep)/width
        xe = 1 - .5*xsep/width
        ys = (rowheight*rows + .5*ysep)/height
        ye = .5*ysep/height
        for colidx in range(0,cols+1):
            ax.axvline(xoff - .5*xsep + colidx*(maxxmax+xsep),ymin=ys,ymax=ye,color='gray')
        for rowidx in range(0,rows+1):
            ax.axhline(rowidx*(maxymax+ysep) + .5*ysep,xmin=xs,xmax=xe,color='gray')
    if showPlateLabels:
        if plateformat == '200honeycomb':
            for rowidx in range(0,rows,2):
                ax.text(.5*xoff - .5*xsep,((rowidx+1)*rowheight),str(rows-rowidx),verticalalignment='center',horizontalalignment='center')
            for colidx in range(0,cols,2):
                ax.text(xoff + ((colidx+.5)*colwidth),.5*yoff + (rows+0.5)*rowheight,str(10*colidx+1),
                        verticalalignment='center',horizontalalignment='center')
            for colidx in range(1,cols,2):
                ax.text(xoff + ((colidx+.5)*colwidth),.5*yoff + (rows+0.5)*rowheight - (maxymax+ysep)/2,str(10*colidx+1),
                        verticalalignment='center',horizontalalignment='center')
        else:
            rowlabels=[chr(x) for x in range(ord('A'), ord('A') + rows)]
            rowlabels.reverse()
            for colidx in range(0,cols):
                ax.text(xoff + ((colidx+.5)*colwidth),.5*yoff + rows*rowheight,str(colidx+1),verticalalignment='center',horizontalalignment='center')
            for rowidx in range(0,rows):
                ax.text(.5*xoff - .5*xsep,((rowidx+.5)*rowheight),rowlabels[rowidx],verticalalignment='center',horizontalalignment='center')
    return rectangles

def _overviewWellOffsets(tc,tcidx,cols,rows,xoff,yoff,xsep,ysep,colwidth,rowheight,width,height,maxxmax,maxymax,plateformat='96'):
    if plateformat == '200honeycomb':
        (colidx,rowidx)=divmod(tcidx,rows)
        rowidx = rows-rowidx-1
        extrayoff = 0 if colidx%2 else (maxymax+ysep)/2

        return {
            'wellxoff':      xoff +  colidx   *(maxxmax+xsep),
            'x1': -.5*xsep + xoff +  colidx   *(maxxmax+xsep),
            'x2': -.5*xsep + xoff + (colidx+1)*(maxxmax+xsep),
            'wellyoff':   .5*ysep +  rowidx   *(maxymax+ysep) + extrayoff,
            'y1':         .5*ysep +  rowidx   *(maxymax+ysep) + extrayoff,
            'y2':         .5*ysep + (rowidx+1)*(maxymax+ysep) + extrayoff,
            'relx1': (- .5*xsep + xoff + colwidth* colidx   )/width,
            'relx2': (- .5*xsep + xoff + colwidth*(colidx+1))/width,
            'rely1': (.5*ysep + rowheight* rowidx    + extrayoff)/height,
            'rely2': (.5*ysep + rowheight*(rowidx+1) + extrayoff)/height,
            }

    else:
        (rowidx,colidx)=divmod(tcidx,cols)
        rowidx = rows-rowidx-1

        return {
            'wellxoff':      xoff +  colidx   *(maxxmax+xsep),
            'x1': -.5*xsep + xoff +  colidx   *(maxxmax+xsep),
            'x2': -.5*xsep + xoff + (colidx+1)*(maxxmax+xsep),
            'wellyoff':   .5*ysep +  rowidx   *(maxymax+ysep),
            'y1':         .5*ysep +  rowidx   *(maxymax+ysep),
            'y2':         .5*ysep + (rowidx+1)*(maxymax+ysep),
            }



def replicateNderivativePlotAsPdf(pdfout,replicate,
                                  creator=None,
                                  addToTitle=None,
                                  showTitle=True,
                                  showDerivatives=None,
                                  showRaw=True,showBackground=True,showSingle=True,showSmoothed=True,
                                  showMaxLinearSlope=False,
                                  showLogod=True,showLogodSmoothed=False,
                                  showMaxGrowthrate=True, showMaxGrowthrateFromLogOdDerivative=True,
                                  showDerivativeLinear=True, showSmoothedDerivativeLinear=True,
                                  showLogOdDerivative=True,showLogOdDerivativeFromNonLog=True,showLogOdDerivativeFromNonLogSmoothed=True,
                                  showExpFitsOd0Mu=True,showExpFitsMu=False,
                                  showGrowthyield=True):
    """
    Show features (such as OD, log(OD), derivatives, maximal growth, ...) on a matplotlib figure.

    FIXME enhance documentation of parameters.
    """
    with contextlib.closing(matplotlib.backends.backend_pdf.PdfPages(pdfout)) as pdfp:
        fig = matplotlib.figure.Figure()
        canvas = matplotlib.backends.backend_pdf.FigureCanvasPdf(fig)
    
        status=replicateToFig(fig,replicate,
                              addToTitle=addToTitle,
                              showTitle=showTitle,
                              showDerivatives=showDerivatives,
                              showRaw=showRaw,
                              showBackground=showBackground,
                              showSingle=showSingle,
                              showSmoothed=showSmoothed,
                              showLogod=showLogod,
                              showLogodSmoothed=showLogodSmoothed,
                              showMaxGrowthrate=showMaxGrowthrate,
                              showMaxGrowthrateFromLogOdDerivative=showMaxGrowthrateFromLogOdDerivative,
                              showLogOdDerivative=showLogOdDerivative,
                              showLogOdDerivativeFromNonLog=showLogOdDerivativeFromNonLog,
                              showLogOdDerivativeFromNonLogSmoothed=showLogOdDerivativeFromNonLogSmoothed,
                              showExpFitsOd0Mu=showExpFitsOd0Mu,
                              showExpFitsMu=showExpFitsMu,
                              showGrowthyield=showGrowthyield,
                              showMaxLinearSlope=showMaxLinearSlope,
                              showDerivativeLinear=showDerivativeLinear,
                              showSmoothedDerivativeLinear=showSmoothedDerivativeLinear)
        fig.savefig(pdfp,format="pdf",bbox_inches='tight')
    
        # set some metadata
        d = pdfp.infodict()
        d['Title'] = replicate.sampleid+' '+replicate.condition
        if creator is not None:
            d['Creator'] = creator

def plotReplicatesToPdfPages(plate,pdfp,
                             listOfReplicates=[],
                             progresscnt=0,
                             showTitle=True,
                             showDerivatives=None,
                             showWellString=False,
                             showRaw=True,showBackground=True,showSingle=True,showSmoothed=True,
                             showMaxLinearSlope=False,
                             showLogod=True,showLogodSmoothed=False,
                             showMaxGrowthrate=True, showMaxGrowthrateFromLogOdDerivative=True,
                             showLogOdDerivative=True,showLogOdDerivativeFromNonLog=True,showLogOdDerivativeFromNonLogSmoothed=True,
                             showDerivativeLinear=True, showSmoothedDerivativeLinear=True,
                             showExpFitsOd0Mu=True,showExpFitsMu=False,
                             showGrowthyield=True,
                             progressCall=None):
    """
    Create a multi-page pdf with many properties in plots.

    :param pdfout: Filename.
    :type pdfout: string
    :param progressCall: Function that will be called on each iteration.
    :type progressCall: @fun(int)

    FIXME enhance documentation of parameters.
    """

    for tc in listOfReplicates:
        if progressCall is not None:
            progressCall(progresscnt)
        progresscnt+=1
        addlabel=None
        if showWellString:
            addlabel=tc.activeChildWellIdStr()

        figtc = matplotlib.figure.Figure()
        canvas = matplotlib.backends.backend_pdf.FigureCanvasPdf(figtc)
        status=replicateToFig(figtc,tc,addToTitle=addlabel,
                              showTitle=showTitle,
                              showDerivatives=showDerivatives,
                              showRaw=showRaw,showBackground=showBackground,showSingle=showSingle,showSmoothed=showSmoothed,
                              showLogod=showLogod,showLogodSmoothed=showLogodSmoothed,
                              showMaxGrowthrate=showMaxGrowthrate,showMaxGrowthrateFromLogOdDerivative=showMaxGrowthrateFromLogOdDerivative,
                              showLogOdDerivative=showLogOdDerivative,
                              showLogOdDerivativeFromNonLog=showLogOdDerivativeFromNonLog,
                              showLogOdDerivativeFromNonLogSmoothed=showLogOdDerivativeFromNonLogSmoothed,
                              showExpFitsOd0Mu=showExpFitsOd0Mu,showExpFitsMu=showExpFitsMu,
                              showGrowthyield=showGrowthyield,
                              showMaxLinearSlope=showMaxLinearSlope,
                              showDerivativeLinear=showDerivativeLinear,
                              showSmoothedDerivativeLinear=showSmoothedDerivativeLinear)
        figtc.savefig(pdfp,format="pdf",bbox_inches='tight')

    return progresscnt

def plotFullOdPlate(plate,pdfout,creator=None,
                    showDerivatives=None,
                    showWellString=False,
                    showReplicateGroups=True,
                    includeBackground=False,
                    showTitle=True,
                    showRaw=True,showBackground=True,showSingle=True,showSmoothed=True,
                    showMaxLinearSlope=False,
                    showLogod=True,showLogodSmoothed=False,
                    showMaxGrowthrate=True, showMaxGrowthrateFromLogOdDerivative=True,
                    showLogOdDerivative=True,showLogOdDerivativeFromNonLog=True,showLogOdDerivativeFromNonLogSmoothed=True,
                    showDerivativeLinear=True, showSmoothedDerivativeLinear=True,
                    showExpFitsOd0Mu=True,showExpFitsMu=False,
                    showGrowthyield=True,
                    progressCall=None):
    """
    Create a multi-page pdf with many properties in plots.

    :param pdfout: Filename.
    :type pdfout: string
    :param progressCall: Function that will be called on each iteration.
    :type progressCall: @fun(int)

    FIXME enhance documentation of parameters.
    """

    with contextlib.closing(matplotlib.backends.backend_pdf.PdfPages(pdfout)) as pdfp:
        figall = matplotlib.figure.Figure()
        canvas = matplotlib.backends.backend_pdf.FigureCanvasPdf(figall)
        ax=figall.add_subplot(111)
    
        if showReplicateGroups:
            if includeBackground:
                listOfReplicates=plate.replicateGroups
            else:
                listOfReplicates=plate.nonBackgroundReplicates()
        else:
            if includeBackground:
                listOfReplicates=plate.wells
            else:
                listOfReplicates=plate.nonBackgroundWells()
    
        plotReplicatesToPdfPages(plate,pdfp,
                                 listOfReplicates=listOfReplicates,
                                 showTitle=showTitle,
                                 showDerivatives=showDerivatives,
                                 showWellString=showWellString,
                                 showRaw=showRaw,
                                 showBackground=showBackground,
                                 showSingle=showSingle,
                                 showSmoothed=showSmoothed,
                                 showMaxLinearSlope=showMaxLinearSlope,
                                 showLogod=showLogod,
                                 showLogodSmoothed=showLogodSmoothed,
                                 showMaxGrowthrate=showMaxGrowthrate,
                                 showMaxGrowthrateFromLogOdDerivative=showMaxGrowthrateFromLogOdDerivative,
                                 showLogOdDerivative=showLogOdDerivative,
                                 showLogOdDerivativeFromNonLog=showLogOdDerivativeFromNonLog,
                                 showLogOdDerivativeFromNonLogSmoothed=showLogOdDerivativeFromNonLogSmoothed,
                                 showDerivativeLinear=showDerivativeLinear,
                                 showSmoothedDerivativeLinear=showSmoothedDerivativeLinear,
                                 showExpFitsOd0Mu=showExpFitsOd0Mu,
                                 showExpFitsMu=showExpFitsMu,
                                 showGrowthyield=showGrowthyield,
                                 progressCall=progressCall)
    
        # set some metadata
        d = pdfp.infodict()
        d['Title'] = plate.plateId
        if creator is not None:
            d['Creator'] = creator
