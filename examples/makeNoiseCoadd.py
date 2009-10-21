#!/usr/bin/env python
from __future__ import with_statement
"""Make a coadd from gaussian noise images
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

BaseDir = os.path.dirname(__file__)
PolicyPath = os.path.join(BaseDir, "makeNoiseCoadd_policy.paf")

def makeNoiseExposure(shape, sigma, variance=1.0):
    """Make a gaussian noise Exposure
    
    Inputs:
    - shape: shape of output array (cols, rows)
    - sigma; sigma of image distribution
    - variance: constant value for variance plane
    """
    image = numpy.random.normal(loc=0.0, scale=sigma, size=shape)
    mask = numpy.zeros(shape, dtype=int)
    variance = numpy.zeros(shape, dtype=float) + variance
    
    maskedImage = imTestUtils.maskedImageFromArrays((image, mask, variance))
    exposure = afwImage.ExposureF(maskedImage)
    return exposure
    

if __name__ == "__main__":
#    pexLog.Trace.setVerbosity('lsst.coadd', 5)
    helpStr = """Usage: makeCoadd.py coaddPath numImages

where:
- coaddPath is the desired name or path of the output coadd
- numImages is the desired number of images

The policy controlling the parameters is %s
""" % (PolicyPath,)
    if len(sys.argv) != 3:
        print helpStr
        sys.exit(0)
    
    coaddPath = sys.argv[1]
    coaddBaseName, coaddExt = os.path.splitext(coaddPath)
    weightOutName = "%s_weight.fits" % (coaddBaseName,)
    
    numImages = int(sys.argv[2])

    policy = pexPolicy.Policy.createPolicy(PolicyPath)

    saveDebugImages = policy.getBool("saveDebugImages")
    imageShape = (256, 256)
    imageShape = policy.getArray("imageShape")
    print "imageShape=%r" % (imageShape,)
    imageSigma = policy.getDouble("imageSigma")
    variance = policy.getDouble("variance")
    
    sys.stdout.write("""
coaddPath  = %s
numImages  = %d
imageShape = %s
imageSigma = %0.1f
variance   = %0.1f
""" % (coaddPath, numImages, imageShape, imageSigma, variance))
    
    numpy.random.seed(0)
    
    coadd = None
    for imInd in range(numImages):
        print "Create exposure %d" % (imInd,)
        exposure = makeNoiseExposure(shape=imageShape, sigma=imageSigma, variance=variance)
        
        if not coadd:
            print "Create coadd with exposure %d" % (imInd,)
            coadd = coaddChiSq.Coadd(exposure, policy)
        else:
            print "Add exposure %d to coadd" % (imInd,)
            coadd.addExposure(exposure)
        if saveDebugImages:
            expName = "%d_%s" % (imInd, coaddPath)
            print "Save intermediate exposure %s" % (expName,)
            exposure.writeFits(expName)

    print "Save weight map as %s" % (weightOutName,)
    weightMap = coadd.getWeightMap()
    weightMap.writeFits(weightOutName)
    coaddExposure = coadd.getCoadd()
    print "Save coadd as %s" % (coaddPath,)
    coaddExposure.writeFits(coaddPath)
