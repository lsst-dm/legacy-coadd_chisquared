#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

from __future__ import with_statement
"""Plot a histogram for a chi squared coadd and overlay a chi squared distribution
"""
import os
import sys
import math

import numpy
import pyfits
import matplotlib.pyplot as pyplot

NBins = 500
UseLogForY = False
UseSqrtForX = False

ChiSqOffsets = (0.0, 0.25, 0.50, 0.75, 1.0)


def clipOutliers(arr):
    """Clip values 3 sigma outside the median
    where sigma is estimated as the inter-quartile range * 0.741
    """
    arr.sort()
    arrLen = len(arr)
    iqr = arr[arrLen * 3 / 4] - arr[arrLen / 4]
    threeSigma = 4 * iqr * 0.741
    median = arr[arrLen / 2]
    minGood = median - threeSigma
    maxGood = median + threeSigma
    return numpy.extract((arr >= minGood) & (arr <= maxGood), arr)


def plotHistogram(coaddName, weightMapName):
    """Plot a histogram given paths to the coadd and weight map
    """
    coadd = pyfits.open(coaddName)
    weightMap = pyfits.open(weightMapName)
    weightMapData = weightMap[0].data
    chiSqOrder = weightMapData.max()
    coaddData = coadd[0].data
    if coaddData == None:  # handle MEF
        coaddData = coadd[1].data
    if coaddData.shape != weightMapData.shape:
        raise RuntimeError("Image shape = %s != %s = weight map shape" %
                           (coaddData.shape, weightMapData.shape))
    goodData = numpy.extract(weightMapData.flat == chiSqOrder, coaddData.flat)
    numWrongOrder = len(coaddData.flat) - len(goodData)
    tempLen = len(goodData)
    goodData = numpy.extract(numpy.isfinite(goodData), goodData)
    numNotFinite = tempLen - len(goodData)

    # undo normalization
    goodData *= float(chiSqOrder)
    # get rid of large values -- these are clearly not noise
    tempLen = len(goodData)
    goodData = numpy.extract(goodData < 50, goodData)
    numBig = tempLen - len(goodData)
    numTotal = len(coaddData.flat)
    print "ChiSquared order = %d; %d good pixels; %0.1f%% had wrong order; %0.1f%% were not finite; %0.1f%% >= 50" % \
        (chiSqOrder, len(goodData), numWrongOrder * 100.0 / numTotal,
         numNotFinite * 100.0 / numTotal, numBig * 100.0 / numTotal)

    hist, binEdges = numpy.histogram(goodData, bins=NBins)
    hist = numpy.array(hist, dtype=float)
    hist /= hist.sum()

    if UseLogForY:
        dataY = numpy.log10(hist)
    else:
        dataY = hist

    dataX = binEdges[0:-1]
    if UseSqrtForX:
        plotDataX = numpy.sqrt(dataX)
    else:
        plotDataX = dataX

    plotNameSet = []

    # plot histogram: log10(frequency) vs. value
    plotNameSet.append((pyplot.plot(plotDataX, dataY, drawstyle="steps"), "Data %0d" % (chiSqOrder,)))
    pyplot.suptitle("Histogram for Chi Squared Coadd %s" % (os.path.basename(coaddName,)))
    if UseLogForY:
        pyplot.ylabel('log10 frequency')
    else:
        pyplot.ylabel('frequency')

    if UseSqrtForX:
        pyplot.xlabel('sqrt of sum of (counts/noise)^2')
    else:
        pyplot.xlabel('sum of (counts/noise)^2')

    maxYInd = None
    for chiSqFudge in ChiSqOffsets:
        fudgedOrder = chiSqOrder + chiSqFudge

        # plot chiSq probability distribution
        chiSqX = dataX
        chiSqDist = numpy.power(chiSqX, (fudgedOrder / 2.0) - 1) * numpy.exp(-chiSqX / 2.0)
        chiSqDist /= chiSqDist.sum()

        # normalize chiSqDist to match data
        endInd = chiSqDist.argmax()  # index to peak of chi squared distribution
        if chiSqFudge == 0.0:
            maxYInd = endInd
        startInd = endInd / 2
        scaleArr = hist[startInd:endInd] / chiSqDist[startInd:endInd]
        chiSqDist *= scaleArr.mean()

        if UseLogForY:
            chiSqDistY = numpy.log10(chiSqDist)
        else:
            chiSqDistY = chiSqDist
        plotNameSet.append((pyplot.plot(plotDataX, chiSqDistY), "ChiSq %0.1f" % (fudgedOrder,)))

    plots, plotNames = zip(*plotNameSet)
    pyplot.legend(plots, plotNames, loc=0)

    # set plot limits
    # compute min and max, ignoring non-finite values
    finiteYValues = numpy.extract(numpy.isfinite(dataY), dataY)
    minY = finiteYValues.min()
    maxY = finiteYValues.max()
    # compute min of tail (portion after maximum), ignoring non-finite values
    tailMinY = numpy.extract(numpy.isfinite(dataY[maxYInd:]), dataY[maxYInd:]).min()
    yRange = maxY - tailMinY
    # plot out to where tail falls to 1% of max value
    yEndVal = tailMinY + (yRange * 0.01)
    endInd = numpy.where(dataY[maxYInd:] <= yEndVal)[0][0] + maxYInd
    pyplot.xlim((0, plotDataX[endInd]))
    yMargin = yRange * 0.05
    pyplot.ylim((minY, maxY + yMargin))

    pyplot.show()

if __name__ == "__main__":
    helpStr = """Usage: plotHistogram.py coaddfile

where:
- coaddfile is the path of the coadd image
"""
    if len(sys.argv) != 2:
        print helpStr
        sys.exit(0)

    coaddName = sys.argv[1]
    for suffix in ("_img.fits", ".fits"):
        if coaddName.lower().endswith(suffix):
            weightMapName = coaddName[0:-len(suffix)] + "_weight.fits"
            break
    else:
        raise RuntimeError("Cannot find weight map")

    plotHistogram(coaddName, weightMapName)
