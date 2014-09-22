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
"""Test Coadd class
"""
import unittest

import numpy

import lsst.utils.tests as utilsTests
import lsst.pex.policy as pexPolicy
import lsst.afw.image as afwImage
import lsst.afw.image.utils as imageUtils
import lsst.afw.image.testUtils as afwTestUtils
import lsst.coadd.chisquared as coaddChiSq

doPlot = False

if doPlot:
    import matplotlib.pyplot as pyplot

def makeHistogram(coadd, numBins, numImages):
    """Generate a histogram for a given coadd maskedImage

    Inputs:
    - coadd: a chiSquared coadd MaskedImage
    - numBins: number of bins for histogram
    - numImages: number of images that went into the coadd
    
    Returns:
    - histX: x values for histogram of coadd data (counts)
    - histY: y values for histogram of coadd data (number of pixels)
    - chiSqY: chi squared distribution values corresponding to histX
    """
    coaddData = coadd.getImage().getArray()
    # undo normalization
    coaddData *= float(numImages)
    # get rid of nans and infs
    goodData = numpy.extract(numpy.isfinite(coaddData.flat), coaddData.flat)
    goodData = numpy.extract(goodData < 50, goodData)
    
    # compute histogram
    histY, binEdges = numpy.histogram(goodData, bins=numBins)
    histX = binEdges[0:-1]
    histY = numpy.array(histY, dtype=float) # convert from int to float
    histY /= histY.sum()

    # compute chiSq probability distribution; chi squared order = numImages
    chiSqY = numpy.power(histX, (numImages / 2.0) - 1) * numpy.exp(-histX / 2.0)
    chiSqY /= chiSqY.sum()
    
    return (histX, histY, chiSqY)

class CoaddTestCase(unittest.TestCase):
    def testNoiseCoadd(self):
        """Build a coadd from noise images and compare the histogram to a chi squared distribution
        """
        numImages = 4
        imShape = (150, 150)
        imSigma = 1.0
        imVariance = 1.0
        numBins = 200
        maxStdDevErr = 0.2
        maxMeanErr = 1.0e-12
        
        badMaskPlanes = ["EDGE"]
    
        numpy.random.seed(0)
        
        coadd = None
        wcs = None
        for imInd in range(numImages):
            maskedImage = afwTestUtils.makeGaussianNoiseMaskedImage(
                dimensions=imShape, sigma=imSigma, variance=imVariance)
            exposure = afwImage.ExposureF(maskedImage, wcs)
            
            if not coadd:
                coadd = coaddChiSq.Coadd(
                    bbox = exposure.getBBox(),
                    wcs = exposure.getWcs(),
                    badMaskPlanes = badMaskPlanes)
    
            coadd.addExposure(exposure)
    
        coadd.getWeightMap()
        coaddExposure = coadd.getCoadd()
        
        histX, histY, chiSqY = makeHistogram(coaddExposure.getMaskedImage(), numBins=numBins, numImages=numImages)

        if doPlot:
            pyplot.plot(histX, histY, drawstyle="steps")
            pyplot.plot(histX, chiSqY)
            pyplot.ylabel('frequency')
            pyplot.xlabel('sum of (counts/noise)^2')
            pyplot.show()
        
        errArr = (histY - chiSqY) * numBins # crude scaling
#         print "Mean error =", errArr.mean()
#         print "Std dev error = ", errArr.std()
        if errArr.std() > maxStdDevErr:
            self.fail("Standard deviation of error = %s > %s limit" % (errArr.std(), maxStdDevErr))
        if errArr.mean() > maxMeanErr:
            self.fail("Mean of error = %s > %s limit" % (errArr.mean(), maxMeanErr))

    def testFilters(self):
        """Test that the coadd filter is set correctly
        """
        filterPolicyFile = pexPolicy.DefaultPolicyFile("afw", "SdssFilters.paf", "tests")
        filterPolicy = pexPolicy.Policy.createPolicy(filterPolicyFile, filterPolicyFile.getRepositoryPath(), True)
        imageUtils.defineFiltersFromPolicy(filterPolicy, reset=True)
        
        unkFilter = afwImage.Filter()
        gFilter = afwImage.Filter("g")
        rFilter = afwImage.Filter("r")

        imShape = (150, 150)
        imSigma = 1.0
        imVariance = 1.0
        
        badMaskPlanes = ["EDGE"]
    
        numpy.random.seed(0)
        
        coadd = None
        wcs = None
        maskedImage = afwTestUtils.makeGaussianNoiseMaskedImage(
            dimensions=imShape, sigma=imSigma, variance=imVariance)
        inExp = afwImage.ExposureF(maskedImage, wcs)

        coadd = coaddChiSq.Coadd(
            bbox = inExp.getBBox(),
            wcs = inExp.getWcs(),
            badMaskPlanes = badMaskPlanes,
        )

        inExp.setFilter(gFilter)
        coadd.addExposure(inExp)
        self.assertEqualFilters(coadd.getCoadd().getFilter(), gFilter)
        self.assertEqualFilterSets(coadd.getFilters(), (gFilter,))
        coadd.addExposure(inExp)
        self.assertEqualFilters(coadd.getCoadd().getFilter(), gFilter)
        self.assertEqualFilterSets(coadd.getFilters(), (gFilter,))
        
        inExp.setFilter(rFilter)
        coadd.addExposure(inExp)
        self.assertEqualFilters(coadd.getCoadd().getFilter(), unkFilter)
        self.assertEqualFilterSets(coadd.getFilters(), (gFilter, rFilter))
    
    def assertEqualFilters(self, f1, f2):
        """Compare two filters
        
        Right now compares only the name, but if == ever works for Filters (ticket #1744)
        then use == instead
        """
        self.assertEquals(f1.getName(), f2.getName())
    
    def assertEqualFilterSets(self, fs1, fs2):
        """Assert that two collections of filters are equal, ignoring order
        """
        self.assertEquals(set(f.getName() for f in fs1), set(f.getName() for f in fs2))

def suite():
    """Return a suite containing all the test cases in this module.
    """
    utilsTests.init()

    suites = [
        unittest.makeSuite(CoaddTestCase),
        unittest.makeSuite(utilsTests.MemoryTestCase),
    ]

    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
