#!/usr/bin/env python
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
UseLogForY = True
UseSqrtForX = False


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
    if coaddData.shape != weightMapData.shape:
        raise RuntimeError("Image shape = %s != %s = weight map shape" % \
            (coaddData.shape, weightMapData.shape))
    goodData = numpy.extract(weightMapData.flat == chiSqOrder, coaddData.flat)
    numWrongOrder = len(coaddData.flat) - len(goodData)
    tempLen = len(goodData)
    goodData = numpy.extract(numpy.isfinite(goodData), goodData)
    numNotFinite = tempLen - len(goodData)

    # undo normalization
    goodData *= float(chiSqOrder)
    # get rid of large values -- clearly not noise
    tempLen = len(goodData)
    goodData = numpy.extract(goodData < 50, goodData)
    numBig = tempLen - len(goodData)
    print "ChiSquared order = %d; %d usable pixels; %d have wrong order; %d are not finite; %d >= 50" % \
        (chiSqOrder, len(goodData), numWrongOrder, numNotFinite, numBig)
    
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

    # plot histogram: log10(frequency) vs. value
    pyplot.plot(plotDataX, dataY, drawstyle="steps")
    if UseLogForY:
        pyplot.ylabel('log10 frequency')
    else:
        pyplot.ylabel('frequency')

    if UseSqrtForX:
        pyplot.xlabel('sqrt of sum of (counts/noise)^2')
    else:
        pyplot.xlabel('sum of (counts/noise)^2')

    # plot chiSq probability distribution
    chiSqX = dataX
    chiSqDist = numpy.power(chiSqX, (chiSqOrder / 2.0) - 1) * numpy.exp(-chiSqX / 2.0)
    chiSqDist /= chiSqDist.sum()
    if UseLogForY:
        chiSqDistY = numpy.log10(chiSqDist)
    else:
        chiSqDistY = chiSqDist
    pyplot.plot(plotDataX, chiSqDistY)

    # set plot limits
    goodY = numpy.extract(numpy.isfinite(dataY), dataY)
    minY = goodY.min()
    maxY = goodY.max()
    maxYInd = goodY.argmax()
    tailMinY = goodY[maxYInd:].min()
    yRange = maxY - tailMinY
    # plot out to where tail falls to 1% of max value
    yEndVal = tailMinY + (yRange * 0.01)
    endInd = numpy.where(goodY[maxYInd:] <= yEndVal)[0][0] + maxYInd
    endInd = len(goodY)-1
    pyplot.xlim((0, plotDataX[endInd]))
    yMargin = yRange * 0.05
    pyplot.ylim((minY, maxY + yMargin))
    
    pyplot.show()


if __name__ == "__main__":
    helpStr = """Usage: plotHistogram.py coaddfile chiSqOrder

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
