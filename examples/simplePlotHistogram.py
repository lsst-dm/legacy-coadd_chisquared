#!/usr/bin/env python
from __future__ import with_statement
"""Plot a histogram of the counts in an image
"""
import os
import sys
import math

import numpy
import pyfits
import matplotlib.pyplot as pyplot

NBins = 300

if __name__ == "__main__":
    helpStr = """Usage: simplePlotHistogram.py image

where:
- image is the path of the image
"""
    if len(sys.argv) != 2:
        print helpStr
        sys.exit(0)
    
    imagePath = sys.argv[1]
    
image = pyfits.open(imagePath)
imageData = image[0].data
# get rid of nans and take square root
goodData = numpy.ma.array(imageData.flat, mask=numpy.isnan(imageData.flat)).compressed()
hist, binEdges = numpy.histogram(goodData, bins=NBins)
counts = binEdges[0:-1]

# plot log10(frequency) vs. sqrt of value
pyplot.plot(counts, hist, drawstyle="steps")
pyplot.ylabel('frequency')
pyplot.xlabel('counts')

pyplot.show()
