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
"""
import os
import sys

import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.coadd.utils as coaddUtils
import lsst.coadd.chisquared as coaddChiSq

SaveDebugImages = False

BackgroundCells = 256

PackageName = "coadd_chisquared"
PolicyDictName = "chiSquaredCoadd_dict.paf"

def subtractBackground(maskedImage, doDisplay = False):
    """Subtract the background from a MaskedImage
    
    Note: at present the mask and variance are ignored, but they might used be someday.
    
    Returns the background object returned by afwMath.makeBackground.
    """
    if doDisplay:
        ds9.mtv(maskedImage)
    bkgControl = afwMath.BackgroundControl(afwMath.Interpolate.NATURAL_SPLINE)
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
    helpStr = """Usage: makeCoadd.py coaddPath indata [policy]

where:
- coaddPath is the desired name or path of the output coadd
- indata is a file containing a list of:
    pathToExposure
  where:
  - pathToExposure is the path to an Exposure (without the final _img.fits)
  - the first exposure listed is taken to be the reference exposure,
    which determines the size and WCS of the coadd
  - empty lines and lines that start with # are ignored.
- policy: path to a policy file

The policy dictionary is: policy/%s
""" % (PolicyDictName,)
    if len(sys.argv) not in (3, 4):
        print helpStr
        sys.exit(0)
    
    outName = sys.argv[1]
    if os.path.exists(outName + "_img.fits"):
        print "Coadd file %s already exists" % (outName,)
        print helpStr
        sys.exit(1)
    weightOutName = outName + "_weight.fits"
    
    indata = sys.argv[2]
    
    if len(sys.argv) > 3:
        policyPath = sys.argv[3]
        policy = pexPolicy.Policy(policyPath)
    else:
        policy = pexPolicy.Policy()

    policyFile = pexPolicy.DefaultPolicyFile(PackageName, PolicyDictName, "policy")
    defPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)
    policy.mergeDefaults(defPolicy.getDictionary())

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
            if SaveDebugImages:
                exposure.writeFits("bgsub%s" % (fileName,))

            print "Warp and add to coadd"
            warpedExposure = coadd.addExposure(exposure)
            if SaveDebugImages:
                warpedExposure.writeFits("warped%s" % (fileName,))

    weightMap = coadd.getWeightMap()
    weightMap.writeFits(weightOutName)
    coaddExposure = coadd.getCoadd()
    coaddExposure.writeFits(outName)
