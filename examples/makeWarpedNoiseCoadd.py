#!/usr/bin/env python
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

BaseDir = os.path.dirname(__file__)
PolicyPath = os.path.join(BaseDir, "makeWarpedNoiseCoadd_policy.paf")

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
  - pathToExposure is the path to an Exposure (without the final _img.fits);
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
    if os.path.exists(coaddPath + "_img.fits"):
        print "Coadd file %s already exists" % (coaddPath,)
        print helpStr
        sys.exit(1)
    coaddPath = sys.argv[1]
    coaddBaseName, coaddExt = os.path.splitext(coaddPath)
    weightOutName = "%s_weight.fits" % (coaddBaseName,)
    
    indata = sys.argv[2]

    policy = pexPolicy.Policy.createPolicy(PolicyPath)

    saveDebugImages = policy.getBool("saveDebugImages")
    imageSigma = policy.getDouble("imageSigma")
    variance = policy.getDouble("variance")
    
    sys.stdout.write("""
coaddPath  = %s
imageSigma = %0.1f
variance   = %0.1f
saveDebugImages = %s
""" % (coaddPath, imageSigma, variance, saveDebugImages))
    
    numpy.random.seed(0)

    # process exposures
    ImageSuffix = "_img.fits"
    coadd = None
    with file(indata, "rU") as infile:
        for lineNum, line in enumerate(infile):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            filePath = line
            fileName = os.path.basename(filePath)
            if not os.path.isfile(filePath + ImageSuffix):
                print "Skipping exposure %s; image file %s not found" % (fileName, filePath + ImageSuffix,)
                continue
            
            print "Processing exposure %s" % (filePath,)
            inputExposure = afwImage.ExposureF(filePath)
            imageShape = tuple(inputExposure.getMaskedImage().getDimensions())
            wcs = inputExposure.getWcs()
            print "Create Gaussian noise exposure"
            maskedImage = makeNoiseMaskedImage(shape=imageShape, sigma=imageSigma, variance=variance)
            exposure = afwImage.ExposureF(maskedImage, wcs)
            
            if not coadd:
                print "Create coadd using this exposure as a reference"
                coadd = coaddChiSq.Coadd(exposure, policy)
                if saveDebugImages:
                    warpedExposure = coadd.getWarpedReferenceExposure()
                    warpedExposure.writeFits("warped%s" % (fileName,))
            else:
                print "Warp and add this exposure to the coadd"
                # note: we are not psf-matching so the returned psfMatchedExposure
                # is identical to the returned warpedExposure
                warpedExposure, psfMatchedExposure = coadd.addExposure(exposure)[0:2]
                if saveDebugImages:
                    warpedExposure.writeFits("warped%s" % (fileName,))

    print "Save resulting coadd and weight map"
    weightMap = coadd.getWeightMap()
    weightMap.writeFits(weightOutName)
    coaddExposure = coadd.getCoadd()
    coaddExposure.writeFits(coaddPath)
