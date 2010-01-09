# -*- python -*-
#
# Setup our environment
#
import os.path
import lsst.SConsUtils as scons

env = scons.makeEnv(
    "coadd_chisquared",
    r"$HeadURL: svn+ssh://svn.lsstcorp.org/DMS/coadd/chisquared/trunk/SConstruct $",
    [
        ["boost", "boost/version.hpp", "boost_system:C++"],
        ["boost", "boost/version.hpp", "boost_filesystem:C++"],
        ["boost", "boost/regex.hpp", "boost_regex:C++"],
        ["boost", "boost/serialization/base_object.hpp", "boost_serialization:C++"],
        ["python", "Python.h"],
        ["cfitsio", "fitsio.h", "m cfitsio", "ffopen"],
        ["wcslib", "wcslib/wcs.h", "m wcs"],
        ["xpa", "xpa.h", "xpa", "XPAPuts"],
        ["minuit2", "Minuit2/FCNBase.h", "Minuit2:C++"],
        ["gsl", "gsl/gsl_rng.h", "gslcblas gsl"],
        ["pex_exceptions", "lsst/pex/exceptions.h", "pex_exceptions:C++"],
        ["utils", "lsst/utils/Utils.h", "utils:C++"],
        ["daf_base", "lsst/daf/base.h", "daf_base:C++"],
        ["pex_logging", "lsst/pex/logging/Trace.h", "pex_logging:C++"],
        ["security", "lsst/security/Security.h", "security:C++"],
        ["pex_policy", "lsst/pex/policy/Policy.h", "pex_policy:C++"],
        ["daf_persistence", "lsst/daf/persistence.h", "daf_persistence:C++"],
        ["daf_data", "lsst/daf/data.h", "daf_data:C++"],
        ["eigen", "Eigen/Core.h"],
        ["afw", "lsst/afw/image.h", "afw:C++"],
        ["coadd_utils", "lsst/coadd/utils.h", "coadd_utils:C++"],
    ],
)
env.libs["coadd_chisquared"] += env.getlibs("boost wcslib cfitsio minuit2 gsl utils daf_base daf_data daf_persistence pex_exceptions pex_logging pex_policy security afw coadd_utils")

#
# Build/install things
#
for d in Split("doc examples lib python/lsst/coadd/chisquared tests"):
    SConscript(os.path.join(d, "SConscript"))

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

Alias("install", [
    env.InstallAs(os.path.join(env['prefix'], "doc", "doxygen"), os.path.join("doc", "htmlDir")),
    env.Install(env['prefix'], "examples"),
    env.Install(env['prefix'], "include"),
    env.Install(env['prefix'], "lib"),
    env.Install(env['prefix'], "policy"),
    env.Install(env['prefix'], "python"),
    env.Install(env['prefix'], "src"),
    env.Install(env['prefix'], "tests"),
    env.InstallEups(env['prefix'] + "/ups"),
])

scons.CleanTree(r"*~ core *.so *.os *.o *.pyc")

files = scons.filesToTag()
if files:
    env.Command("TAGS", files, "etags -o $TARGET $SOURCES")

env.Declare()
env.Help("""
LSST implementation of PSF-matching coaddition algorithm
""")
