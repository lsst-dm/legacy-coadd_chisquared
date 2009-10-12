// -*- LSST-C++ -*-
#ifndef LSST_COADD_CHISQUARED_ADDTOCOADD_H
#define LSST_COADD_CHISQUARED_ADDTOCOADD_H
/**
* @file
*
* @author Russell Owen
*/
#include "boost/cstdint.hpp"

#include "lsst/afw/image.h"

namespace lsst {
namespace coadd {
namespace chisquared {

    template<typename CoaddPixelT, typename WeightPixelT>
    void addToCoadd(
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
