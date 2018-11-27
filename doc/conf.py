"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documenation builds.
"""

from documenteer.sphinxconfig.stackconf import build_package_configs
import lsst.dia.pipe


_g = globals()
_g.update(build_package_configs(
    project_name='dia_pipe',
    version=lsst.dia.pipe.version.__version__))
