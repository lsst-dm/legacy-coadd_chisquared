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

class Coadd(object):
    def __init__(self, dimensions, wcs, policy):
        """Create a chi-squared coadd
        
        Inputs:
        - dimensions: dimensions of coadd
        - wcs: WCS of coadd
        - policy: Parameters include:
            - doWarpExposures: if True then warp each exposure to match the reference exposure;
                otherwise assume this has already been done.
            - allowedMaskPlanes: a list of space-separated mask names of bits that are allowed in the coadd
            
            The following is only read if doWarpExposures is True:
            - warpingKernelName: name of warping kernel (see lsst.afw.math.makeWarpingKernel)
        """
        self._log = pexLog.Log(pexLog.Log.getDefaultLog(), "coadd.chisquared.Coadd")
        self._policy = policy
        self._doWarpExposures = policy.get("doWarpExposures")

        allowedMaskPlanes = policy.get("allowedMaskPlanes").split()
        allowedMask = 0
        for maskPlaneName in allowedMaskPlanes:
            allowedMask |= 1 << afwImage.MaskU.getMaskPlane(maskPlaneName)
        self._badPixelMask = 0xFFFF - allowedMask

        self._wcs = wcs # for convenience
        blankMaskedImage = afwImage.MaskedImageF(dimensions)
        self._coadd = afwImage.ExposureF(blankMaskedImage, wcs)

        if self._doWarpExposures:
            self._warpingKernel = afwMath.makeWarpingKernel(policy.get("warpingKernelName"))
        self._wcs = self._coadd.getWcs() # merely a convenience
        self._weightMap = afwImage.ImageF(self._coadd.getMaskedImage().getDimensions(), 0)

    def addExposure(self, exposure):
        """Add an Exposure to the coadd
        
        If doWarpExposures is True then first warp the exposure to match the WCS of the coadd.
        Otherwise the exposure is assumed to already have been warped.
        
        Inputs:
        - exposure: Exposure to add to coadd; must have the background subtracted.
            
        Returns:
        - warpedExposure: exposure warped to match the WCS of the coadd,
            or the original exposure if doWarpExposure false
        """
        if self._doWarpExposures:
            # warp exposure
            warpedExposure = self._warpExposure(exposure)
        else:
            warpedExposure = exposure
            
        self._basicAddExposure(warpedExposure)
        return warpedExposure

    def getCoadd(self):
        """Get the coadd Exposure, as computed so far
        
        Return the coadd Exposure consisting of the reference exposure and all exposures
        you have added so far. You may call addExposure and getCoadd as often as you like.
        """
        # make a deep copy so I can scale it
        coaddMaskedImage = self._coadd.getMaskedImage()
        scaledMaskedImage = coaddMaskedImage.__class__(coaddMaskedImage, True)

        # set the edge pixels
        coaddUtils.setCoaddEdgeBits(scaledMaskedImage.getMask(), self._weightMap)
        
        # scale non-edge pixels by weight map
        coaddUtils.divide(scaledMaskedImage, self._weightMap)
        
        return afwImage.makeExposure(scaledMaskedImage, self._wcs)
        
    def getWeightMap(self):
        """Get the weight map
        
        The weight map is a float Image of the same dimensions as the coadd;
        the value of each pixel is the number of input images
        that contributed to the associated pixel of the coadd.
        """
        return self._weightMap
    
    def _basicAddExposure(self, warpedExposure):
        """Add a an exposure to the coadd; it is assumed to have the same WCS as the coadd

        Inputs:
        - warpedExposure: Exposure to add to coadd; must have the background subtracted
            and have been warped to match the WCS of the coadd.
        - weight: weight of good pixels for the weight map
        """
        self._log.log(pexLog.Log.INFO, "add masked image to coadd")
        chisquaredLib.addToCoadd(self._coadd.getMaskedImage(), self._weightMap,
            warpedExposure.getMaskedImage(), self._badPixelMask, 1.0)

    def _warpExposure(self, exposure):
        """Warp an exposure to match the WCS of the coadd
        
        _coadd and _warpingKernel must have been set.
        """
        self._log.log(pexLog.Log.INFO, "warp exposure")
        warpedExposure = coaddUtils.makeBlankExposure(self._coadd)
        afwMath.warpExposure(warpedExposure, exposure, self._warpingKernel)
        return warpedExposure
