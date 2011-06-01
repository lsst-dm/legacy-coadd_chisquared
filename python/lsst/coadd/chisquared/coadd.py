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
    def __init__(self, bbox, wcs, badMaskPlanes):
        """Create a chi-squared coadd
        
        Inputs:
        @param[in] bbox: bounding box of coadd Exposure with respect to parent (afwGeom.Box2I):
            coadd dimensions = bbox.getDimensions(); xy0 = bbox.getMin()
        @param[in] wcs: WCS of coadd exposure (lsst.afw.math.Wcs)
        @param[in] badMaskPlanes: mask planes to pay attention to when rejecting masked pixels.
            Specify as a collection of names.
            badMaskPlanes should always include "EDGE".
        """
        coaddUtils.Coadd.__init__(self,
            bbox = bbox,
            wcs = wcs,
            badMaskPlanes = badMaskPlanes,
            coaddZeroPoint = 27.0, # this value is ignored
            logName = "coadd.chisquared.Coadd",
        )
    
    @classmethod
    def fromPolicy(cls, bbox, wcs, policy):
        """Create a coadd
        
        Inputs:
        @param[in] bbox: bounding box of coadd Exposure with respect to parent (afwGeom.Box2I):
            coadd dimensions = bbox.getDimensions(); xy0 = bbox.getMin()
        @param[in] wcs: WCS of coadd exposure (lsst.afw.math.Wcs)
        @param[in] policy: coadd policy; see policy/CoaddPolicyDictionary.paf
        """
        badMaskPlanes = policy.getArray("badMaskPlanes")
        return cls(bbox, wcs, badMaskPlanes)

    def addExposure(self, exposure, weightFactor=1.0):
        """Add a an exposure to the coadd; it is assumed to have the same WCS as the coadd

        Inputs:
        @param[in] exposure: Exposure to add to coadd; must be background-subtracted and warped to match the coadd.
        @param[in] weightFactor: weight with which to add exposure to coadd

        @return
        - overlapBBox: region of overlap between exposure and coadd in parent coordinates (afwGeom.Box2I)
        - weight: weight with which exposure was added to coadd; weight = weightFactor for this kind of coadd
        """
        self._log.log(pexLog.Log.INFO, "add exposure to coadd")
        overlapBBox = chisquaredLib.addToCoadd(self._coadd.getMaskedImage(), self._weightMap,
            exposure.getMaskedImage(), self._badPixelMask, weightFactor)

        return overlapBBox, weightFactor
    
    def _setCalib(self, *dumArgs):
        """Set _coaddZeroPoint and add a Calib object to self._coadd
        
        Chi-squared coadds have no calibration so make this a null operation
        """
        pass
