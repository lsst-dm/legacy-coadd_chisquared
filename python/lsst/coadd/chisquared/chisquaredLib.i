// -*- lsst-c++ -*-

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
 
%define chisquaredLib_DOCSTRING
"
Python interface to lsst::coadd::chisquared functions and classes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.coadd.chisquared", docstring=chisquaredLib_DOCSTRING) chisquaredLib

%{
#include "boost/cstdint.hpp"
#include "lsst/pex/logging.h"
#include "lsst/afw/math.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/coadd/chisquared.h"
%}

%include "lsst/p_lsstSwig.i"
%import  "lsst/afw/image/imageLib.i" 
%import  "lsst/afw/math/mathLib.i" 

%lsst_exceptions()

%include "lsst/coadd/chisquared/addToCoadd.h"
%define %ADDTOCOADD(COADDPIXEL, WEIGHTPIXEL)
    %template(addToCoadd) lsst::coadd::chisquared::addToCoadd<COADDPIXEL, WEIGHTPIXEL>;
%enddef
%ADDTOCOADD(double, double);
%ADDTOCOADD(double, float);
%ADDTOCOADD(double, int);
%ADDTOCOADD(double, boost::uint16_t);
%ADDTOCOADD(float, double);
%ADDTOCOADD(float, float);
%ADDTOCOADD(float, int);
%ADDTOCOADD(float, boost::uint16_t);
