#!/usr/bin/env python
from __future__ import with_statement
"""
This example requires:
- A set of science exposures
- A file containing the paths to each, as:
  exposure1
  exposure2
  ...
The first exposure's WCS and size are used for the coadd.

@todo: modify to use ChiSquaredCoaddDictionary as the base dictionary
"""
import os
import sys
import math

import numpy

import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.coadd.utils as coaddUtils
import lsst.coadd.chisquared as coaddChiSq

BaseDir = os.path.dirname(__file__)
PolicyPath = os.path.join(BaseDir, "makeCoadd_policy.paf")

BackgroundCells = 256

def subtractBackground(maskedImage, doDisplay = False):
    """Subtract the background from a MaskedImage
    
    Note: at present the mask and variance are ignored, but they might used be someday.
    
    Returns the background object returned by afwMath.makeBackground.
    """
    if doDisplay:
        ds9.mtv(maskedImage)
    bkgControl = afwMath.BackgroundControl(afwMath.NATURAL_SPLINE)
    bkgControl.setNxSample(int(maskedImage.getWidth() // BackgroundCells) + 1)
    bkgControl.setNySample(int(maskedImage.getHeight() // BackgroundCells) + 1)
    bkgControl.sctrl.setNumSigmaClip(3)
    bkgControl.sctrl.setNumIter(3)

    image = maskedImage.getImage()
    bkgObj = afwMath.makeBackground(image, bkgControl)
    image -= bkgObj.getImageF()
    if doDisplay:
        ds9.mtv(image)
    return bkgObj

if __name__ == "__main__":
    pexLog.Trace.setVerbosity('lsst.coadd', 5)
    helpStr = """Usage: makeCoadd.py coaddPath indata

where:
- coaddPath is the desired name or path of the output coadd
- indata is a file containing a list of:
    pathToExposure
  where:
  - pathToExposure is the path to an Exposure (without the final _img.fits)
  - the first exposure listed is taken to be the reference exposure,
    which determines the size and WCS of the coadd
  - empty lines and lines that start with # are ignored.

The policy controlling the parameters is %s
""" % (PolicyPath,)
    if len(sys.argv) != 3:
        print helpStr
        sys.exit(0)
    
    outName = sys.argv[1]
    if os.path.exists(outName + "_img.fits"):
        print "Coadd file %s already exists" % (outName,)
        print helpStr
        sys.exit(1)
    weightOutName = outName + "_weight.fits"
    
    indata = sys.argv[2]

    policy = pexPolicy.Policy.createPolicy(PolicyPath)

    saveDebugImages = policy.getBool("saveDebugImages")

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
            exposure = afwImage.ExposureF(filePath)
            
            if not coadd:
                print "Create coadd"
                coadd = coaddChiSq.Coadd(exposure.getMaskedImage().getDimensions(), exposure.getWcs(), policy)
            
            print "Subtract background"
            subtractBackground(exposure.getMaskedImage())
            if saveDebugImages:
                exposure.writeFits("bgsub%s" % (fileName,))

            print "Warp and add to coadd"
            warpedExposure = coadd.addExposure(exposure)
            if saveDebugImages:
                warpedExposure.writeFits("warped%s" % (fileName,))

    weightMap = coadd.getWeightMap()
    weightMap.writeFits(weightOutName)
    coaddExposure = coadd.getCoadd()
    coaddExposure.writeFits(outName)
