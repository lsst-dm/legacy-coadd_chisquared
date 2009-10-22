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

NBins = 200
UseLogForY = False

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
goodData = numpy.sqrt(numpy.ma.array(coaddData.flat, mask=numpy.isnan(coaddData.flat)).compressed())
hist, binEdges = numpy.histogram(goodData, bins=NBins)
if UseLogForY:
    dataY = numpy.log10(numpy.where(hist > 1.0, hist, 1.0))
else:
    dataY = hist

dataX = binEdges[0:-1]

# plot log10(frequency) vs. sqrt of value
pyplot.plot(dataX, dataY, drawstyle="steps")
if UseLogForY:
    pyplot.ylabel('log10 frequency')
else:
    pyplot.ylabel('frequency')
pyplot.xlabel('sqrt of sum of (counts/noise)^2')

# plot chiSq probability distribution
x = dataX
chiSqDist = numpy.power(x, (nImages / 2.0) - 1) * numpy.exp(-x / 2.0)
# need some way to scale it properly, but for now just chuck in a factor
# that will get it in the ballpark
chiSqDist *= 50
if UseLogForY:
    chiSqDistY = numpy.log10(numpy.where(chiSqDist > 1.0, chiSqDist, 1.0))
else:
    chiSqDistY = chiSqDist
pyplot.plot(dataX, chiSqDistY)

pyplot.show()
