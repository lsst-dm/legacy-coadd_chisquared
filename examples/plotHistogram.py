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

NBins = 1000

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
# get rid of nans and take square root
goodData = numpy.sqrt(numpy.ma.array(coaddData.flat, mask=numpy.isnan(coaddData.flat)).compressed())
hist, binEdges = numpy.histogram(goodData, bins=NBins)
logHist = numpy.log10(hist)
logHist = numpy.where(logHist >= 0, logHist, 0)

degFree = len(goodData)

sqrtValue = binEdges[0:-1]

# plot log10(frequency) vs. sqrt of value
pyplot.plot(sqrtValue, logHist, drawstyle="steps")
pyplot.ylabel('log10 frequency')
pyplot.xlabel('sqrt of sum of (counts/noise)^2')

# plot chiSq probability distribution
chiSqDist = numpy.power(sqrtValue, (nImages / 2.0) - 1) * numpy.exp(-sqrtValue / 2.0)
# need some way to scale it properly, but for now just chuck in a factor
# that will get it in the ballpark
logChiSqDist = numpy.log10(chiSqDist * 1000)
logChiSqDist = numpy.where(logChiSqDist >= 0, logChiSqDist, 0)

pyplot.plot(sqrtValue, logChiSqDist)

pyplot.show()
