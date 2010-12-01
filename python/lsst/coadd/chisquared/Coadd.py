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

"""chi-squared coadd
"""
import sys
import lsst.pex.logging as pexLog
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDetect
import lsst.coadd.utils as coaddUtils
import chisquaredLib

__all__ = ["Coadd"]

class Coadd(coaddUtils.Coadd):
    def __init__(self, dimensions, wcs, policy):
        """Create a chi-squared coadd
        
        Inputs:
        - dimensions: dimensions of coadd
        - wcs: WCS of coadd
        - policy: see policy/chiSquaredCoadd_dict.paf
        """
        coaddUtils.Coadd.__init__(self,
            dimensions = dimensions,
            wcs = wcs,
            policy = policy,
            logName = "coadd.chisquared.Coadd",
        )

    def addExposure(self, exposure, weightFactor=1.0):
        """Add a an exposure to the coadd; it is assumed to have the same WCS as the coadd

        Inputs:
        - exposure: Exposure to add to coadd; must be background-subtracted and warped to match the coadd.
        - weight: weight of good pixels for the weight map
        """
        self._log.log(pexLog.Log.INFO, "add exposure to coadd")
        chisquaredLib.addToCoadd(self._coadd.getMaskedImage(), self._weightMap,
            exposure.getMaskedImage(), self._badPixelMask, weightFactor)
