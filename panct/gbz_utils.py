"""
Utilities for dealing with GBZ files
"""

import logging
import os
from shutil import which
import subprocess
import tempfile

def ExtractRegionFromGBZ(gbz_file, region, reference):
    tmpfile = tempfile.NamedTemporaryFile(delete=False)
    cmd = ["query", "--sample", reference, \
        "--contig", region[0],
        "--interval", str(region[1])+".."+str(region[2]), \
        gbz_file]
    proc = subprocess.run(cmd, stdout=tmpfile)
    if proc.returncode != 0:
        return None
    else:
        return tmpfile
    
def CheckGBZBaseInstalled(log):
    """
    Check that gbz2db and query from
    gbz-base are installed

    Returns
    -------
    passed : bool
       True if both tools are found
    """
    if which("gbz2db") is None:
        log.warning("Could not find gbz2db")
        return False
    if which("query") is None:
        log.critical("Could not find query")
        return False
    return True

def IndexGBZ(gbz_file: str):
    """
    Index the GBZ file with gbz2db

    Parameters
    ----------
    gbz_file : str
        Path to the GBZ file

    Returns
    -------
    passed : bool
        True if we were able to create the .gbz.db file
    """
    cmd = ["gbz2db", gbz_file]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE)
    return proc.returncode == 0

def CheckGBZFile(gbz_file: str, log: logging.Logger):
    """
    Check if the GBZ file exists and is
    indexed by GBZ-Base

    Parameters
    ----------
    gbz_file : str
        Path to the GBZ file
    log : logging.Logger

    Returns
    -------
    passed : bool
        True if GBZ file and GBZ-base database exist
    """
    if not os.path.exists(gbz_file):
        log.critical(f"{gbz_file} does not exist\n")
        return False
    if not os.path.exists(gbz_file + ".db"):
        log.info(f"{gbz_file}.db does not exist. Attempting to create")
        if not IndexGBZ(gbz_file):
            log.critical("Failed to create GBZ index")
            return False
    return True
