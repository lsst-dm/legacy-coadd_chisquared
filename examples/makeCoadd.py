#!/usr/bin/env python
from __future__ import with_statement
"""
This example requires:
- A set of science exposures
- A file containing the paths to each, as:
  exposure1 (the reference exposure)
  exposure2
  ...
The first exposure is the reference exposure: this should have the worst seeing
(because all other exposures are PSF-matched to this one). It is also the exposure
whose WCS and size are used for the coadd.

@todo: modify to use ip_diffim dictionary as the base policy for psf-matching once Policy supports this.
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
import lsst.coadd.chisquared as coaddChiSquared

BaseDir = os.path.dirname(__file__)
DefPolicyPath = os.path.join(BaseDir, "makeCoadd_policy.paf")

def subtractBackground(maskedImage, doDisplay = False):
    """Subtract the background from a MaskedImage
    
    Note: at present the mask and variance are ignored, but they might used be someday.
    
    Returns the background object returned by afwMath.makeBackground.
    """
    if doDisplay:
        ds9.mtv(maskedImage)
    bkgControl = afwMath.BackgroundControl(afwMath.NATURAL_SPLINE)
    bkgControl.setNxSample(max(2, int(maskedImage.getWidth()/256) + 1))
    bkgControl.setNySample(max(2, int(maskedImage.getHeight()/256) + 1))
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
    helpStr = """Usage: makeCoadd.py coaddfile indata

where:
- coaddfile is the desired name or path of the output coadd
- indata is a file containing a list of:
    pathToExposure
  where:
  - pathToExposure is the path to an Exposure (without the final _img.fits)
  - the first exposure listed is taken to be the reference exposure;
    this one should have the worst PSF of the data set.
  - empty lines and lines that start with # are ignored.

The policy controlling the parameters is makeCoadd_policy.paf
See makeCoadd_policy.paf for documentation
"""
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

    makeCoaddPolicyPath = DefPolicyPath
    makeCoaddPolicy = pexPolicy.Policy.createPolicy(makeCoaddPolicyPath)

    saveDebugImages = makeCoaddPolicy.getBool("saveDebugImages")

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
            
            print "Subtract background"
            subtractBackground(exposure.getMaskedImage())
            
            if not coadd:
                print "First exposure is the reference: warp but do not psf-match"
                coadd = coaddChiSquared.Coadd(exposure, makeCoaddPolicy)
                if saveDebugImages:
                    warpedExposure = coadd.getWarpedReferenceExposure()
                    warpedExposure.writeFits("warped%s" % (fileName,))
            else:
                print "Warp, psf-match and add to coadd"
                warpedExposure = coadd.addExposure(exposure)
                if saveDebugImages:
                    warpedExposure.writeFits("warped%s" % (fileName,))

    weightMap = coadd.getWeightMap()
    weightMap.writeFits(weightOutName)
    coaddExposure = coadd.getCoadd()
    coaddExposure.writeFits(outName)
