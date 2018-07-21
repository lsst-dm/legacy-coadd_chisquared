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

#ifndef LSST_COADD_CHISQUARED_ADDTOCOADD_H
#define LSST_COADD_CHISQUARED_ADDTOCOADD_H
/**
 * @file
 *
 * @author Russell Owen
 */
#include "lsst/afw/geom.h"
#include "lsst/afw/image.h"

namespace lsst {
namespace coadd {
namespace chisquared {

/**
 * @brief add good pixels from a masked image to a coadd and associated weight map
 * using the chi squared algorithm
 *
 * For good pixels (image.mask & badPixelMask == 0), coadd and weightMap are altered as follows:
 * coadd.image += image.image**2 / image.variance
 * coadd.mask |= image.mask
 * weightMap += weight
 * For bad pixels, coadd and weightMap are not altered.
 *
 * Note that coadd.variance is not altered.
 *
 * @return overlapBBox: overlapping bounding box, relative to parent image (hence xy0 is taken into account)
 *
 * @throw pexExcept::InvalidParameterError if coadd and weightMap dimensions or xy0 do not match.
 */
template <typename CoaddPixelT, typename WeightPixelT>
lsst::afw::geom::Box2I addToCoadd(
        lsst::afw::image::MaskedImage<CoaddPixelT, lsst::afw::image::MaskPixel,
                                      lsst::afw::image::VariancePixel>
                &coadd,                                    ///< [in,out] coadd to be modified
        lsst::afw::image::Image<WeightPixelT> &weightMap,  ///< [in,out] weight map to be modified
        lsst::afw::image::MaskedImage<CoaddPixelT, lsst::afw::image::MaskPixel,
                                      lsst::afw::image::VariancePixel> const
                &maskedImage,  ///< masked image to add to coadd
        lsst::afw::image::MaskPixel const
                badPixelMask,  ///< skip input pixel if input mask & badPixelMask !=0
        WeightPixelT weight    ///< relative weight of this image
);

}  // namespace chisquared
}  // namespace coadd
}  // namespace lsst

#endif  // !defined(LSST_COADD_CHISQUARED_ADDTOCOADD_H)
