"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documenation builds.
"""

from documenteer.sphinxconfig.stackconf import build_package_configs
import lsst.coadd.chisquared


_g = globals()
_g.update(build_package_configs(
    project_name='coadd_chisquared',
    version=lsst.coadd.chisquared.version.__version__))
