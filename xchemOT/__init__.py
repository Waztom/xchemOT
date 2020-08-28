"""
xchemOT
 Automate generation of synthetic workflows for the OpenTrons robotics platform
"""

# Add imports here
from .xchemOT import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
