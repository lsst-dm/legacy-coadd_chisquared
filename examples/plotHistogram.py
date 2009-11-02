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
UseLogForY = False
UseSqrtForX = False

if __name__ == "__main__":
    helpStr = """Usage: plotHistogram.py coaddfile nImages

where:
- coaddfile is the path of the coadd
"""
    if len(sys.argv) != 3:
        print helpStr
        sys.exit(0)
    
    coaddName = sys.argv[1]
    nImages = int(sys.argv[2])
    
coadd = pyfits.open(coaddName)
coaddData = coadd[0].data
# undo normalization
coaddData *= float(nImages)
# get rid of nans and take square root
goodData = numpy.extract(numpy.isfinite(coaddData.flat), coaddData.flat)
hist, binEdges = numpy.histogram(goodData, bins=NBins)
if UseLogForY:
    dataY = numpy.log10(numpy.where(hist > 1.0, hist, 1.0))
else:
    dataY = hist

dataX = binEdges[0:-1]
if UseSqrtForX:
    plotDataX = numpy.sqrt(dataX)
else:
    plotDataX = dataX

# plot log10(frequency) vs. sqrt of value
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
chiSqDist = numpy.power(chiSqX, (nImages / 2.0) - 1) * numpy.exp(-chiSqX / 2.0)

# fit scale by fitting points up to peak in histogram (crude but it's a start)
maxPos = hist.argmax()
goodVals = hist[0:maxPos] > 1
scaleFactors = numpy.extract(goodVals, hist[0:maxPos] / chiSqDist[0:maxPos])
weights = numpy.extract(goodVals, hist[0:maxPos])
meanScaleFactor = numpy.average(scaleFactors, weights=weights)
chiSqDist *= meanScaleFactor
if UseLogForY:
    chiSqDistY = numpy.log10(numpy.where(chiSqDist > 1.0, chiSqDist, 1.0))
else:
    chiSqDistY = chiSqDist
pyplot.plot(plotDataX, chiSqDistY)

pyplot.show()
