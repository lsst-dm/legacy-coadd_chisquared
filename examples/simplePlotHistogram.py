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
if imageData == None:  # handle MEF
    imageData = image[1].data
# get rid of nans and take square root
goodData = numpy.ma.array(imageData.flat, mask=numpy.isnan(imageData.flat)).compressed()
hist, binEdges = numpy.histogram(goodData, bins=NBins)
counts = binEdges[0:-1]

# plot log10(frequency) vs. sqrt of value
pyplot.plot(counts, hist, drawstyle="steps")
pyplot.ylabel('frequency')
pyplot.xlabel('counts')

pyplot.show()
