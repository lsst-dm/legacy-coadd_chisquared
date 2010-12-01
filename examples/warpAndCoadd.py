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

import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9
import lsst.coadd.utils as coaddUtils
import lsst.coadd.chisquared as coaddChiSq

SaveDebugImages = False

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
        print "Coadd file %r already exists" % (coaddPath,)
        print helpStr
        sys.exit(1)
    weightPath = os.splitext(coaddPath)[0] + "_weight.fits"
    
    indata = sys.argv[2]
    
    if len(sys.argv) > 3:
        policyPath = sys.argv[3]
        policy = pexPolicy.Policy(policyPath)
    else:
        policy = pexPolicy.Policy()

    policyFile = pexPolicy.DefaultPolicyFile(PolicyPackageName, PolicyDictName, "policy")
    defPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)
    policy.mergeDefaults(defPolicy.getDictionary())
    warpingKernelName = policy.getPolicy("warpPolicy").get("warpingKernelName")
    allowedMaskPlanes = policy.getPolicy("coaddPolicy").get("allowedMaskPlanes")

    # process exposures
    coadd = None
    with file(indata, "rU") as infile:
        for lineNum, line in enumerate(infile):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            filePath = line
            fileName = os.path.basename(filePath)
            if not os.path.isfile(filePath):
                print "Skipping exposure %s; file %r not found" % (fileName, filePath)
                continue
            
            print "Processing exposure %s" % (filePath,)
            exposure = afwImage.ExposureF(filePath)
            
            if not coadd:
                print "First exposure is the reference; create warper and coadd"
                maskedImage = exposure.getMaskedImage()
                warper = coaddUtils.Warp(warpingKernelName)
                coadd = coaddChiSq.Coadd(maskedImage.getDimensions(), exposure.getWcs(), allowedMaskPlanes)
            
                print "Add reference exposure to coadd (without warping)"
                coadd.addExposure(exposure)
            else:
                print "Warp exposure"
                warpedExposure = warper.warpExposure(coadd.getDimensions(), coadd.getWcs(), exposure)
                if SaveDebugImages:
                    warpedExposure.writeFits("warped%s" % (fileName,))
                    
                print "Add warped exposure to coadd"
                coadd.addExposure(warpedExposure)

    weightMap = coadd.getWeightMap()
    weightMap.writeFits(weightPath)
    coaddExposure = coadd.getCoadd()
    coaddExposure.writeFits(coaddPath)
