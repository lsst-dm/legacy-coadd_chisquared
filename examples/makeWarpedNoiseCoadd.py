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
"""Make a coadd from warped gaussian noise images
"""
import os
import sys
import math
import optparse

import numpy

import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.afw.image as afwImage
import lsst.afw.image.testUtils as imTestUtils
import lsst.coadd.chisquared as coaddChiSq
import lsst.coadd.utils as coaddUtils

BaseDir = os.path.dirname(__file__)
PolicyPath = os.path.join(BaseDir, "MakeWarpedNoiseCoaddDictionary.paf")

def makeNoiseMaskedImage(shape, sigma, variance=1.0):
    """Make a gaussian noise MaskedImageF
    
    Inputs:
    - shape: shape of output array (cols, rows)
    - sigma; sigma of image distribution
    - variance: constant value for variance plane
    """
    image = numpy.random.normal(loc=0.0, scale=sigma, size=shape)
    mask = numpy.zeros(shape, dtype=int)
    variance = numpy.zeros(shape, dtype=float) + variance
    
    return imTestUtils.maskedImageFromArrays((image, mask, variance), afwImage.MaskedImageF)
    

if __name__ == "__main__":
    pexLog.Trace.setVerbosity('lsst.coadd', 5)
    helpStr = """Usage: makeWarpedNoiseCoadd.py coaddPath numImages

where:
- coaddfile is the desired name or path of the output coadd
- indata is a file containing a list of:
    pathToExposure
  where:
  - pathToExposure is the path to an Exposure;
    only the WCS and image size are used from the exposure; all other data is ignored.
  - the first exposure listed is taken to be the reference exposure;
    all other images are warped to match its WCS.
  - empty lines and lines that start with # are ignored.

Make a chi-squared coadd from a set of warped Gaussian noise images.
The WCSs for the Gaussian images are taken from a list of input images
(which are only used for their WCS).
The intent is to see how correlated noise affects the statistics of a pure noise coadd.

The policy controlling the parameters is %s
""" % (PolicyPath,)
    if len(sys.argv) != 3:
        print helpStr
        sys.exit(0)
    
    coaddPath = sys.argv[1]
    if os.path.exists(coaddPath):
        print "Coadd file %s already exists" % (coaddPath,)
        print helpStr
        sys.exit(1)
    weightPath = os.path.splitext(coaddPath)[0] + "_weight.fits"
    
    indata = sys.argv[2]

    policy = pexPolicy.Policy.createPolicy(PolicyPath)

    saveDebugImages = policy.getBool("saveDebugImages")
    imageSigma = policy.getDouble("imageSigma")
    variance = policy.getDouble("variance")
    warpPolicy = policy.getPolicy("warpPolicy")
    allowedMaskPlanes = policy.getPolicy("coaddPolicy").get("allowedMaskPlanes")
    
    sys.stdout.write("""
coaddPath  = %s
imageSigma = %0.1f
variance   = %0.1f
saveDebugImages = %s
""" % (coaddPath, imageSigma, variance, saveDebugImages))
    
    numpy.random.seed(0)

    # process exposures
    coadd = None
    with file(indata, "rU") as infile:
        for lineNum, line in enumerate(infile):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            filePath = line
            fileName = os.path.basename(filePath)
            
            print "Processing exposure %s" % (filePath,)
            try:
                inputExposure = afwImage.ExposureF(filePath)
            except Exception, e:
                print "Skipping %s: %s" % (filePath, e)
                continue
            imageShape = tuple(inputExposure.getMaskedImage().getDimensions())
            wcs = inputExposure.getWcs()

            print "Create Gaussian noise exposure"
            maskedImage = makeNoiseMaskedImage(shape=imageShape, sigma=imageSigma, variance=variance)
            exposure = afwImage.ExposureF(maskedImage, wcs)
            
            if not coadd:
                print "Create warper and coadd with size and WCS matching the first exposure"
                warper = coaddUtils.Warp.fromPolicy(warpPolicy)
                coadd = coaddChiSq.Coadd(
                    dimensions = maskedImage.getDimensions(),
                    xy0 = exposure.getXY0(),
                    wcs = exposure.getWcs(),
                    allowedMaskPlanes = allowedMaskPlanes)

                if saveDebugImages:
                    exposure.writeFits("warped%s" % (fileName,))

                coadd.addExposure(exposure)
            else:
                warpedExposure = warper.warpExposure(
                    dimensions = coadd.getDimensions(),
                    xy0 = coadd.getXY0(),
                    wcs = coadd.getWcs(),
                    exposure = exposure)
                
                coadd.addExposure(warpedExposure)
                
                if saveDebugImages:
                    warpedExposure.writeFits("warped%s" % (fileName,))

    print "Save resulting coadd and weight map"
    weightMap = coadd.getWeightMap()
    weightMap.writeFits(weightPath)
    coaddExposure = coadd.getCoadd()
    coaddExposure.writeFits(coaddPath)
