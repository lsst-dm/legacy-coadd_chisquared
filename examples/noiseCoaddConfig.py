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
"""Config for making coadds of images of pure noise
"""
import lsst.pex.config as pexConfig
import lsst.afw.math as afwMath
import lsst.coadd.chisquared as coaddChiSq

class NoiseCoaddConfig(pexConfig.Config):
    saveDebugImages = pexConfig.Field(
        dtype = bool,
        doc = "Save warped intermediate images?",
        default = False,
    )
    imageShape = pexConfig.ListField(
        dtype = int,
        doc = "Constant value of variance pixels",
        length = 2,
        default = (256, 256),
    )
    imageSigma = pexConfig.Field(
        dtype = float,
        doc = "Sigma of Gaussian noise for image pixels",
        default = 1.0,
    )
    variance = pexConfig.Field(
        dtype = float,
        doc = "Constant value of variance pixels",
        default = 1.0,
    )
    warp = pexConfig.ConfigField(
        dtype = afwMath.Warper.ConfigClass,
        doc = "Policy to control warping.",
    )
    coadd = pexConfig.ConfigField(
        dtype = coaddChiSq.Coadd.ConfigClass,
        doc = "Policy to control coadd.",
    )
