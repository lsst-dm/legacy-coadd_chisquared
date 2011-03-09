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
#include "boost/cstdint.hpp"

#include "lsst/afw/geom.h"
#include "lsst/afw/image.h"

namespace lsst {
namespace coadd {
namespace chisquared {

    template<typename CoaddPixelT, typename WeightPixelT>
    lsst::afw::geom::Box2I addToCoadd(
        lsst::afw::image::MaskedImage<CoaddPixelT, lsst::afw::image::MaskPixel,
            lsst::afw::image::VariancePixel> &coadd,
        lsst::afw::image::Image<WeightPixelT> &weightMap,
        lsst::afw::image::MaskedImage<CoaddPixelT, lsst::afw::image::MaskPixel,
            lsst::afw::image::VariancePixel> const &maskedImage,
        lsst::afw::image::MaskPixel const badPixelMask,
        WeightPixelT weight
    );

}}} // lsst::coadd::chisquared

#endif // !defined(LSST_COADD_CHISQUARED_ADDTOCOADD_H)
