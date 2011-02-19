// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
/**
* @file
*
* @author Russell Owen
*/
#include <cmath>

#include "lsst/pex/exceptions.h"
#include "lsst/coadd/chisquared/addToCoadd.h"

namespace pexExcept = lsst::pex::exceptions;
namespace afwImage = lsst::afw::image;
namespace afwGeom = lsst::afw::geom;
namespace coaddChiSq = lsst::coadd::chisquared;

/**
* @brief add good pixels from a masked image to a coadd and associated weight map
* using the chi squared algorithm
*
* For good pixels (image.mask & badPixelMask == 0), coadd and weightMap are altered as follows:
* coadd.image += (image.image / sqrt(image.variance))**2
* coadd.mask |= image.mask
* weightMap += weight
* For bad pixels, coadd and weightMap are not altered.
*
* Note that coadd.variance is not altered.
*
* @return overlapBBox: overlapping bounding box, relative to parent image (hence xy0 is taken into account)
*
* @throw pexExcept::InvalidParameterException if coadd and weightMap dimensions or xy0 do not match.
*/
template <typename CoaddPixelT, typename WeightPixelT>
afwGeom::BoxI coaddChiSq::addToCoadd(
    // spell out lsst:afw::image to make Doxygen happy
    lsst::afw::image::MaskedImage<CoaddPixelT, lsst::afw::image::MaskPixel,
        lsst::afw::image::VariancePixel> &coadd,        ///< [in,out] coadd to be modified
    lsst::afw::image::Image<WeightPixelT> &weightMap,   ///< [in,out] weight map to be modified
    lsst::afw::image::MaskedImage<CoaddPixelT, lsst::afw::image::MaskPixel,
        lsst::afw::image::VariancePixel> const &image,  ///< masked image to add to coadd
    lsst::afw::image::MaskPixel const badPixelMask, ///< skip input pixel if input mask | badPixelMask != 0
    WeightPixelT weight ///< relative weight of this image
) {
    typedef typename afwImage::MaskedImage<CoaddPixelT, afwImage::MaskPixel, afwImage::VariancePixel> Coadd;
    typedef typename afwImage::Image<WeightPixelT> WeightMap;

    if (coadd.getDimensions() != weightMap.getDimensions()) {
        throw LSST_EXCEPT(pexExcept::InvalidParameterException,
            (boost::format("coadd and weightMap dimensions differ: %dx%d != %dx%d") %
            coadd.getWidth() % coadd.getHeight() % weightMap.getWidth() % weightMap.getHeight()).str());
    }
    if (coadd.getXY0() != weightMap.getXY0()) {
        throw LSST_EXCEPT(pexExcept::InvalidParameterException,
            (boost::format("coadd and weightMap xy0 differ: %d,%d != %d,%d") %
            coadd.getX0() % coadd.getY0() % weightMap.getX0() % weightMap.getY0()).str());
    }

    afwGeom::BoxI overlapBox = afwGeom::BoxI(
        afwGeom::makePointI(coadd.getX0(), coadd.getY0()),
        afwGeom::makeExtentI(coadd.getWidth(), coadd.getHeight()));
    overlapBox.clip(afwGeom::BoxI(
        afwGeom::makePointI(image.getX0(), image.getY0()),
        afwGeom::makeExtentI(image.getWidth(), image.getHeight())));
    if (overlapBox.isEmpty()) {
        return overlapBox;
    }

    afwImage::BBox coaddSubregion(
        afwImage::PointI(overlapBox.getMinX() - coadd.getX0(), overlapBox.getMinY() - coadd.getY0()),
        overlapBox.getWidth(), overlapBox.getHeight());
    Coadd coaddView(coadd, coaddSubregion, false);
    WeightMap weightMapView(weightMap, coaddSubregion, false);

    afwImage::BBox imageSubregion(
        afwImage::PointI(overlapBox.getMinX() - image.getX0(), overlapBox.getMinY() - image.getY0()),
        overlapBox.getWidth(), overlapBox.getHeight());
    Coadd imageView(image, imageSubregion, false);

    for (int y = 0, endY = imageView.getHeight(); y != endY; ++y) {
        typename Coadd::const_x_iterator imageIter = imageView.row_begin(y);
        typename Coadd::const_x_iterator const imageEndIter = imageView.row_end(y);
        typename Coadd::x_iterator coaddIter = coaddView.row_begin(y);
        typename WeightMap::x_iterator weightMapIter = weightMapView.row_begin(y);
        for (; imageIter != imageEndIter; ++imageIter, ++coaddIter, ++weightMapIter) {
            if ((imageIter.mask() & badPixelMask) == 0) {
                CoaddPixelT value = imageIter.image() / std::sqrt(imageIter.variance());
                coaddIter.image() += value * value;
                coaddIter.mask() |= imageIter.mask();
                *weightMapIter += weight;
            }
        }
    }
    return overlapBox;
}

//
// Explicit instantiations
//
#define MASKEDIMAGE(IMAGEPIXEL) afwImage::MaskedImage<IMAGEPIXEL, \
    afwImage::MaskPixel, afwImage::VariancePixel>
#define INSTANTIATE(COADDPIXEL, WEIGHTPIXEL) \
    template afwGeom::BoxI coaddChiSq::addToCoadd<COADDPIXEL, WEIGHTPIXEL>( \
        MASKEDIMAGE(COADDPIXEL) &coadd, \
        afwImage::Image<WEIGHTPIXEL> &weightMap, \
        MASKEDIMAGE(COADDPIXEL) const &image, \
        afwImage::MaskPixel const badPixelMask, \
        WEIGHTPIXEL weight \
    );

INSTANTIATE(double, double);
INSTANTIATE(double, float);
INSTANTIATE(double, int);
INSTANTIATE(double, boost::uint16_t);
INSTANTIATE(float, double);
INSTANTIATE(float, float);
INSTANTIATE(float, int);
INSTANTIATE(float, boost::uint16_t);
