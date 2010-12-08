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
import lsst.pex.logging as pexLog
import lsst.coadd.utils as coaddUtils
import chisquaredLib

__all__ = ["Coadd"]

class Coadd(coaddUtils.Coadd):
    def __init__(self, bbox, wcs, allowedMaskPlanes):
        """Create a chi-squared coadd
        
        Inputs:
        - bbox: bounding box of coadd Exposure with respect to parent (lsst.afw.geom.BoxI):
            coadd dimensions = bbox.getDimensions(); xy0 = bbox.getMin()
        - wcs: WCS of coadd exposure (lsst.afw.math.Wcs)
        - allowedMaskPlanes: mask planes to allow (ignore) when rejecting masked pixels.
            Specify as a single string containing space-separated names
        """
        coaddUtils.Coadd.__init__(self,
            bbox = bbox,
            wcs = wcs,
            allowedMaskPlanes = allowedMaskPlanes,
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
