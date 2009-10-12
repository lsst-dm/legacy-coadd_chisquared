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
