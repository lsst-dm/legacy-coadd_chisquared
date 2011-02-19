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
import traceback

import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9
import lsst.coadd.utils as coaddUtils
import lsst.coadd.chisquared as coaddChiSq

PolicyPackageName = "coadd_utils"
PolicyDictName = "WarpAndCoaddDictionary.paf"

if __name__ == "__main__":
    pexLog.Trace.setVerbosity('lsst.coadd', 5)
    helpStr = """Usage: warpAndCoadd.py coaddPath indata [policy]

where:
- coaddPath is the desired name or path of the output coadd
- indata is a file containing a list of:
    pathToExposure
  where:
  - pathToExposure is the path to an Exposure
  - the first exposure listed is taken to be the reference exposure,
    which determines the size and WCS of the coadd
  - empty lines and lines that start with # are ignored.
- policy: path to a policy file

The policy dictionary is: policy/%s
""" % (PolicyDictName,)
    if len(sys.argv) not in (3, 4):
        print helpStr
        sys.exit(0)
    
    coaddPath = sys.argv[1]
    if os.path.exists(coaddPath):
        print >> sys.stderr, "Coadd file %s already exists" % (coaddPath,)
        sys.exit(1)
    weightPath = os.path.splitext(coaddPath)[0] + "_weight.fits"
    
    indata = sys.argv[2]
    
    if len(sys.argv) > 3:
        policyPath = sys.argv[3]
        policy = pexPolicy.Policy(policyPath)
    else:
        policy = pexPolicy.Policy()

    policyFile = pexPolicy.DefaultPolicyFile(PolicyPackageName, PolicyDictName, "policy")
    defPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)
    policy.mergeDefaults(defPolicy.getDictionary())
    warpPolicy = policy.getPolicy("warpPolicy")
    coaddPolicy = policy.getPolicy("coaddPolicy")
    saveDebugImages = policy.get("saveDebugImages")
    bboxMin = policy.getArray("bboxMin")
    bboxSize = policy.getArray("bboxSize")
    bbox = afwImage.BBox(afwImage.PointI(bboxMin[0], bboxMin[1]), bboxSize[0], bboxSize[1])
    print >> sys.stderr, "saveDebugImages =", saveDebugImages
    print >> sys.stderr, "BBox =", bbox

    # process exposures
    coadd = None
    expNum = 0
    numExposuresInCoadd = 0
    numExposuresFailed = 0
    with file(indata, "rU") as infile:
        for exposurePath in infile:
            exposurePath = exposurePath.strip()
            if not exposurePath or exposurePath.startswith("#"):
                continue
            expNum += 1

            try:
                print >> sys.stderr, "Processing exposure %s" % (exposurePath,)
                try:
                    exposure = afwImage.ExposureF(exposurePath, 0, bbox)
                except Exception, e:
                    print >> sys.stderr, "Skipping %s: %s" % (exposurePath, e)
                    continue
                if saveDebugImages:
                    exposure.writeFits("exposure%s.fits" % (expNum,))
                
                if not coadd:
                    print >> sys.stderr, "Create warper and coadd with size and WCS matching the first exposure"
                    maskedImage = exposure.getMaskedImage()
                    warper = coaddUtils.Warp.fromPolicy(warpPolicy)
                    coadd = coaddChiSq.Coadd.fromPolicy(
                        bbox = coaddUtils.bboxFromImage(exposure),
                        wcs = exposure.getWcs(),
                        policy = coaddPolicy)
                    print >> sys.stderr, "badPixelMask=", coadd.getBadPixelMask()
                
                    coadd.addExposure(exposure)
                else:
                    warpedExposure = warper.warpExposure(
                        wcs = coadd.getWcs(),
                        exposure = exposure,
                        maxBBox = coadd.getBBox())
                    if saveDebugImages:
                        warpedExposure.writeFits("warped%s.fits" % (expNum,))
                        
                    coadd.addExposure(warpedExposure)

                numExposuresInCoadd += 1
            except Exception, e:
                print >> sys.stderr, "Exposure %s failed: %s" % (exposurePath, e)
                if os.path.exists(exposurePath):
                    traceback.print_exc(file=sys.stderr)
                numExposuresFailed += 1
                continue

    print >> sys.stderr, "Coadded %d exposures and failed %d" % (numExposuresInCoadd, numExposuresFailed)
    weightMap = coadd.getWeightMap()
    weightMap.writeFits(weightPath)
    coaddExposure = coadd.getCoadd()
    coaddExposure.writeFits(coaddPath)
