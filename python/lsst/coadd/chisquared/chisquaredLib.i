// -*- lsst-c++ -*-
%define chisquaredLib_DOCSTRING
"
Python interface to lsst::coadd::chisquared functions and classes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.coadd.chisquared", docstring=chisquaredLib_DOCSTRING) chisquaredLib

%{
#include "boost/cstdint.hpp"
#include "lsst/coadd/chisquared.h"
#include "lsst/afw/geom.h" /* why is this needed?
Without it I see the following when compiling utilsLib_wrap:
python/lsst/coadd/utils/utilsLib_wrap.cc: In function 'void* _p_lsst__afw__geom__ellipses__QuadrupoleTo_p_lsst__afw__geom__ellipses__BaseCore(void*, int*)':
python/lsst/coadd/utils/utilsLib_wrap.cc:13321: error: 'lsst::afw::geom::ellipses' has not been declared
...
*/
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
