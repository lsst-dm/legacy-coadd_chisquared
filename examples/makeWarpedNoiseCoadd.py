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

"""Make a coadd from warped gaussian noise images
"""
from __future__ import print_function
import os
import sys
import traceback

import numpy as np

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.image.testUtils as afwTestUtils
import lsst.coadd.chisquared as coaddChiSq
from lsst.log import Log
from noiseCoaddConfig import NoiseCoaddConfig

if __name__ == "__main__":
    Log.getLogger('coadd').setLevel(Log.DEBUG)
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
"""
    if len(sys.argv) != 3:
        print(helpStr)
        sys.exit(0)

    coaddPath = sys.argv[1]
    if os.path.exists(coaddPath):
        print("Coadd file %s already exists" % (coaddPath,), file=sys.stderr)
        sys.exit(1)
    weightPath = os.path.splitext(coaddPath)[0] + "_weight.fits"

    indata = sys.argv[2]

    config = NoiseCoaddConfig()

    sys.stderr.write("""
coaddPath  = %s
imageSigma = %0.1f
variance   = %0.1f
saveDebugImages = %s
""" % (coaddPath, config.imageSigma, config.variance, config.saveDebugImages))

    np.random.seed(0)

    # process exposures
    coadd = None
    numExposuresInCoadd = 0
    numExposuresFailed = 0
    expNum = 0
    with file(indata, "rU") as infile:
        for lineNum, line in enumerate(infile):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            exposurePath = line
            expNum += 1

            try:
                print("Processing exposure %s" % (exposurePath,), file=sys.stderr)
                inputExposure = afwImage.ExposureF(exposurePath)

                print("Create Gaussian noise exposure", file=sys.stderr)
                maskedImage = afwTestUtils.makeGaussianNoiseMaskedImage(
                    dimensions=inputExposure.getDimensions(),
                    sigma=config.imageSigma,
                    variance=config.variance,
                )
                exposure = afwImage.ExposureF(maskedImage, inputExposure.getWcs())

                if config.saveDebugImages:
                    exposure.writeFits("exposure%d.fits" % (expNum,))

                if not coadd:
                    print("Create warper and coadd with size and WCS matching the first exposure",
                          file=sys.stderr)
                    warper = afwMath.Warper.fromConfig(config.warp)
                    coadd = coaddChiSq.Coadd.fromConfig(
                        bbox=exposure.getBBox(),
                        wcs=exposure.getWcs(),
                        config=config.coadd,
                    )
                    print("badPixelMask=", coadd.getBadPixelMask(), file=sys.stderr)

                    coadd.addExposure(exposure)
                else:
                    print("Warp exposure", file=sys.stderr)
                    warpedExposure = warper.warpExposure(
                        destWcs=coadd.getWcs(),
                        srcExposure=exposure,
                        maxBBox=coadd.getBBox(),
                    )

                    coadd.addExposure(warpedExposure)

                    if config.saveDebugImages:
                        warpedExposure.writeFits("warped%d.fits" % (expNum,))

                numExposuresInCoadd += 1
            except Exception as e:
                print("Exposure %s failed: %s" % (exposurePath, e), file=sys.stderr)
                if os.path.exists(exposurePath):
                    traceback.print_exc(file=sys.stderr)
                numExposuresFailed += 1
                continue

    print("Coadded %d exposures and failed %d" % (numExposuresInCoadd, numExposuresFailed), file=sys.stderr)
    weightMap = coadd.getWeightMap()
    weightMap.writeFits(weightPath)
    coaddExposure = coadd.getCoadd()
    coaddExposure.writeFits(coaddPath)
