# This Source Code Form is subject to the terms of the Oxford Nanopore
# Technologies PLC. Public License Version 1.0. See the LICENSE file.
"""Multi-enzyme digestion for Pore-C / CiFi concatemers.

A drop-in extension of pore-c-py's ``digest`` step that accepts several
restriction enzymes, for protocols that used more than one.
"""
__version__ = "1.0.0"
#: version of pore-c-py this tool mirrors byte-for-byte in single-enzyme mode
UPSTREAM_VERSION = "2.0.6"
