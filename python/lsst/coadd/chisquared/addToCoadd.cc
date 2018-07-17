/*
 * LSST Data Management System
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 * See the COPYRIGHT file
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */
#include "pybind11/pybind11.h"

#include <cstdint>

#include "lsst/coadd/chisquared/addToCoadd.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace coadd {
namespace chisquared {

namespace {

/**
 * Wrap addToCoadd
 *
 * @tparam CoaddPixelT  Pixel type of image plane of coadd and masked image
 * @tparam WeightPixelT  Pixel type of weight map and weight scalar
 * @param mod  pybind11 module
 */
template <typename CoaddPixelT, typename WeightPixelT>
void declareAddToCoadd(py::module& mod) {
    mod.def("addToCoadd", &addToCoadd<CoaddPixelT, WeightPixelT>, "coadd"_a, "weightMap"_a, "maskedImage"_a,
            "badPixelMask"_a, "weight"_a);
}

}  // namespace lsst::coadd::chisquared::<anonymous>

PYBIND11_MODULE(addToCoadd, mod) {
    py::module::import("lsst.afw.geom");
    py::module::import("lsst.afw.image");

    declareAddToCoadd<double, double>(mod);
    declareAddToCoadd<double, float>(mod);
    declareAddToCoadd<double, int>(mod);
    declareAddToCoadd<double, std::uint16_t>(mod);
    declareAddToCoadd<float, double>(mod);
    declareAddToCoadd<float, float>(mod);
    declareAddToCoadd<float, int>(mod);
    declareAddToCoadd<float, std::uint16_t>(mod);
}

}  // chisquared
}  // coadd
}  // lsst
