"""chi-squared coadd
"""
import sys
import lsst.pex.logging as pexLog
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDetect
import lsst.coadd.utils as coaddUtils
import lsst.ip.diffim as ipDiffim
import chisquaredLib
import lsst.coadd.psfmatched as coaddPsfMatch

__all__ = ["Coadd"]

class Coadd(object):
    def __init__(self, referenceExposure, policy):
        """Create a chi-squared coadd
        
        Inputs:
        - referenceExposure: the first exposure to add to the reference.
            If matchPolicy is None then the coadd has the identical wcs and size.
            Otherwise the coadd size is determined by resolutionFactor.
        - policy: Parameters include:
            - doWarpExposures: if True then warp and psf-match each exposure to match the reference exposure;
                otherwise assume this has already been done.
            - allowedMaskPlanes: a list of space-separated mask names of bits that are allowed in the coadd
            
            The following are only read if doWarpExposures is True:
            - warpingKernelName: name of warping kernel (see lsst.afw.math.makeWarpingKernel)
            - resolutionFactor: resolution of coadd (along x or y) relative to referenceMaskedImage
            - psfMatchPolicy: a sub-policy containing:
                - kernelCols/Rows: size of PSF-matching kernel
                - fpNpixMin
                - fpNpixMax
                - fpGrowKsize
                - minCleanFp
                - detThreshold
                - detThresholdScaling
                - detThresholdMin
                - detThresholdType
                ...
        """
        self._log = pexLog.Log(pexLog.Log.getDefaultLog(), "coadd.chisquared.Coadd")
        self._referenceExposure = referenceExposure
        self._policy = policy
        self._doWarpExposures = policy.get("doWarpExposures")

        allowedMaskPlanes = policy.get("allowedMaskPlanes").split()
        allowedMask = 0
        for maskPlaneName in allowedMaskPlanes:
            allowedMask |= 1 << afwImage.MaskU.getMaskPlane(maskPlaneName)
        self._badPixelMask = 0xFFFF - allowedMask

        if self._doWarpExposures:
            self._psfMatchPolicy = policy.get("psfMatchPolicy")
            self._warpingKernel = afwMath.makeWarpingKernel(policy.get("warpingKernelName"))
            resolutionFactor = policy.get("resolutionFactor")
            self._coadd = coaddUtils.makeBlankCoadd(referenceExposure, resolutionFactor)
            self._warpedReferenceExposure = self._warpExposure(referenceExposure)
        else:
            self._coadd = coaddUtils.makeBlankExposure(referenceExposure)
            self._warpedReferenceExposure = referenceExposure
        self._wcs = self._coadd.getWcs() # merely a convenience
        self._weightMap = afwImage.ImageF(self._coadd.getMaskedImage().getDimensions(), 0)
        self._basicAddExposure(self._warpedReferenceExposure)

    def addExposure(self, exposure):
        """Add an Exposure to the coadd
        
        If doWarpExposures is True then first warp the exposure to match the WCS of the coadd,
        and psf-match that to the warped reference Exposure before adding to the result to the coadd.
        Otherwise the exposure is assumed to already have been warped and PSF-matched.
        
        Inputs:
        - exposure: Exposure to add to coadd; must have the background subtracted,
            and if warpExposure is True then it should have a larger PSF than the reference exposure.
            
        Returns:
        - warpedExposure: exposure warped to match the WCS of the coadd;
            the original exposure if warpExposure false
        - psfMatchedExposure: warped Exposure psf-mathed to warped reference Exposure;
            the original exposure if warpExposure false
        """
        if self._doWarpExposures:
            # warp exposure
            warpedExposure = self._warpExposure(exposure)
    
            # psf-match warped exposure to reference exposure
            self._log.log(pexLog.Log.INFO, "psf-match masked image")
            psfMatchedMaskedImage, kernelSum, backgroundModel = coaddPsfMatch.psfMatchMaskedImage(
                self._warpedReferenceExposure.getMaskedImage(), warpedExposure.getMaskedImage(),
                self._psfMatchPolicy)
            addMaskedImage = psfMatchedMaskedImage
            psfMatchedExposure = afwImage.makeExposure(psfMatchedMaskedImage, self._wcs)
            weight = 1.0 / kernelSum
        else:
            warpedExposure = exposure
            psfMatchedExposure = exposure
            weight = 1.0
            
        self._basicAddExposure(psfMatchedExposure, weight)
        return (warpedExposure, psfMatchedExposure)

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
    
    def getWarpedReferenceExposure(self):
        """Get the warped reference exposure (for debugging)
        """
        return self._warpedReferenceExposure

    def _basicAddExposure(self, warpedExposure, weight=1.0):
        """Add a warped and psf-matched exposure to the coadd

        Inputs:
        - warpedExposure: Exposure to add to coadd; must have the background subtracted
            and have been warped and psf-matched to warped reference exposure.
        - weight: weight of good pixels for the weight map
        """
        self._log.log(pexLog.Log.INFO, "add masked image to coadd; weight=%0.2f" % (weight,))
        chisquaredLib.addToCoadd(self._coadd.getMaskedImage(), self._weightMap,
            warpedExposure.getMaskedImage(), self._badPixelMask, weight)

    def _warpExposure(self, exposure):
        """Warp an exposure to match the coadd
        
        _coadd and _warpingKernel must have been set.
        """
        self._log.log(pexLog.Log.INFO, "warp exposure")
        warpedExposure = coaddUtils.makeBlankExposure(self._coadd)
        afwMath.warpExposure(warpedExposure, exposure, self._warpingKernel)
        return warpedExposure
