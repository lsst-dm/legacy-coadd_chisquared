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

template <typename CoaddPixelT, typename WeightPixelT>
afwGeom::Box2I coaddChiSq::addToCoadd(
    // spell out lsst:afw::image to make Doxygen happy
    lsst::afw::image::MaskedImage<CoaddPixelT, lsst::afw::image::MaskPixel,
        lsst::afw::image::VariancePixel> &coadd,
    lsst::afw::image::Image<WeightPixelT> &weightMap,
    lsst::afw::image::MaskedImage<CoaddPixelT, lsst::afw::image::MaskPixel,
        lsst::afw::image::VariancePixel> const &image,
    lsst::afw::image::MaskPixel const badPixelMask,
    WeightPixelT weight
) {
    typedef typename afwImage::MaskedImage<CoaddPixelT, afwImage::MaskPixel, afwImage::VariancePixel> Coadd;
    typedef typename afwImage::Image<WeightPixelT> WeightMap;

    if (coadd.getBBox() != weightMap.getBBox()) {
        throw LSST_EXCEPT(pexExcept::InvalidParameterError,
            (boost::format("coadd and weightMap parent bboxes differ: %s != %s") %
            coadd.getBBox() % weightMap.getBBox()).str());
    }

    afwGeom::Box2I overlapBBox = coadd.getBBox();
    overlapBBox.clip(image.getBBox());
    if (overlapBBox.isEmpty()) {
        return overlapBBox;
    }

    Coadd coaddView(coadd, overlapBBox, afwImage::PARENT, false);
    WeightMap weightMapView(weightMap, overlapBBox, afwImage::PARENT, false);
    Coadd imageView(image, overlapBBox, afwImage::PARENT, false);

    for (int y = 0, endY = imageView.getHeight(); y != endY; ++y) {
        typename Coadd::const_x_iterator imageIter = imageView.row_begin(y);
        typename Coadd::const_x_iterator const imageEndIter = imageView.row_end(y);
        typename Coadd::x_iterator coaddIter = coaddView.row_begin(y);
        typename WeightMap::x_iterator weightMapIter = weightMapView.row_begin(y);
        for (; imageIter != imageEndIter; ++imageIter, ++coaddIter, ++weightMapIter) {
            if ((imageIter.mask() & badPixelMask) == 0) {
                CoaddPixelT value = imageIter.image() * imageIter.image() / imageIter.variance();
                coaddIter.image() += value;
                coaddIter.mask() |= imageIter.mask();
                *weightMapIter += weight;
            }
        }
    }
    return overlapBBox;
}

//
// Explicit instantiations
//
/// \cond
#define MASKEDIMAGE(IMAGEPIXEL) afwImage::MaskedImage<IMAGEPIXEL, \
    afwImage::MaskPixel, afwImage::VariancePixel>
#define INSTANTIATE(COADDPIXEL, WEIGHTPIXEL) \
    template afwGeom::Box2I coaddChiSq::addToCoadd<COADDPIXEL, WEIGHTPIXEL>( \
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
/// \endcond
