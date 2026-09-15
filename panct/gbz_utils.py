"""
Utilities for dealing with GBZ files
"""

import os
import logging
import tempfile
import subprocess
from shutil import which
from pathlib import Path
from typing import Optional
import psutil
import time
import multiprocessing
import resource

from .data import Region
from . import graph_utils as gutils

_log = logging.getLogger(__name__)

class RegionExtractionFailed(Exception):
    """Raised when NodeTable/LinkTable construction failed or exceeded a configured memory cap."""
    pass

def _table_worker(result_queue, table_cls, kwargs, memory_limit_bytes):
    """Runs in a child process. Builds a NodeTable or LinkTable under a memory cap."""
    if memory_limit_bytes is not None:
        _set_mem_limit(memory_limit_bytes)
    try:
        table = table_cls(**kwargs)
        result_queue.put(("ok", table))
    except MemoryError:
        result_queue.put(("memory_error", None))
    except Exception as e:
        result_queue.put(("error", str(e)))
def _pool_worker_init(memory_limit_bytes):
    """Runs once when a persistent worker process starts."""
    if memory_limit_bytes is not None:
        resource.setrlimit(resource.RLIMIT_AS, (memory_limit_bytes, memory_limit_bytes))


def _pool_table_builder(args):
    """Runs inside a persistent worker, once per region call, no re-import cost."""
    table_cls, table_kwargs = args
    try:
        table = table_cls(**table_kwargs)
        return ("ok", table)
    except MemoryError:
        return ("memory_error", None)
    except Exception as e:
        return ("error", str(e))


def make_table_pool(memory_limit_gb: Optional[float]):
    """
    Create a single persistent worker process for building NodeTable/LinkTable
    objects under a memory cap, reused across many regions to avoid paying
    interpreter-startup/import cost per region. Returns None if no cap is set
    (callers should build tables in-process directly in that case).
    """
    if memory_limit_gb is None:
        return None
    ctx = multiprocessing.get_context("spawn")
    return ctx.Pool(
        processes=1,
        initializer=_pool_worker_init,
        initargs=(int(memory_limit_gb * 1024**3),),
    )

def build_table_with_limit(
    table_cls, table_kwargs: dict,
    memory_limit_gb: Optional[float], log: logging.Logger,
    region_label: str = "",
    skipped_bed_file: Optional[Path] = None,
    bed_fields: Optional[tuple] = None,
    pool=None,
):
    """
    Build a NodeTable or LinkTable, optionally capping total memory.
    If `pool` is given (see make_table_pool), reuses that persistent worker
    process instead of spawning a fresh one per call. On failure, if
    skipped_bed_file is given, appends bed_fields to it before raising.
    """
    if memory_limit_gb is None:
        return table_cls(**table_kwargs)

    if pool is not None:
        status, payload = pool.apply(_pool_table_builder, ((table_cls, table_kwargs),))
        exitcode_ok = (status == "ok")
    else:
        # Fallback: spawn a fresh process per call (slow — prefer passing `pool`)
        ctx = multiprocessing.get_context("spawn")
        q = ctx.Queue()
        p = ctx.Process(
            target=_table_worker,
            args=(q, table_cls, table_kwargs, int(memory_limit_gb * 1024**3)),
        )
        p.start()
        p.join()
        status, payload = q.get() if not q.empty() else ("crashed", None)
        exitcode_ok = (p.exitcode == 0 and status == "ok")

    if not exitcode_ok:
        log.warning(
            f"{table_cls.__name__} construction failed for {region_label} "
            f"(status={status}); likely exceeded {memory_limit_gb}GB cap."
        )
        if skipped_bed_file is not None and bed_fields is not None:
            with open(skipped_bed_file, "a") as skipf:
                skipf.write("\t".join(str(x) for x in bed_fields) + "\n")
        raise RegionExtractionFailed(f"{table_cls.__name__} for {region_label} exceeded {memory_limit_gb}GB cap")

    return payload

def run_query_with_memory_limit(cmd, tmpfile, memory_limit_gb, log, poll_interval=0.5):
    """
    Run a subprocess, polling its actual RSS memory usage, and kill it if it
    exceeds memory_limit_gb. Returns a returncode: the process's real
    returncode on normal completion, or -9 (mimicking SIGKILL) if we killed
    it for exceeding the cap.
    """
    proc = subprocess.Popen(cmd, stdout=tmpfile, stderr=subprocess.PIPE)
    limit_bytes = memory_limit_gb * 1024**3
    try:
        ps_proc = psutil.Process(proc.pid)
        while proc.poll() is None:
            try:
                rss = ps_proc.memory_info().rss
            except psutil.NoSuchProcess:
                break
            if rss > limit_bytes:
                log.warning(f"'query' exceeded {memory_limit_gb}GB (RSS={rss/1024**3:.2f}GB); terminating.")
                proc.kill()
                proc.wait()
                return -9, b""
            time.sleep(poll_interval)
        _, stderr = proc.communicate()
        return proc.returncode, stderr
    finally:
        pass


def extract_region_from_gbz(
    gbz_file: Path, region: Region, reference: str,
    memory_limit_gb: Optional[float] = None,
    log: logging.Logger = None,
    skipped_bed_file: Optional[Path] = None,
) -> Optional[Path]:
    if log is None:
        log = _log
    tmpfile = tempfile.NamedTemporaryFile(delete=False)
    cmd = [
        "query",
        "--sample", reference,
        "--contig", region.chrom,
        "--interval", str(region.start) + ".." + str(region.end),
        str(gbz_file) + ".db",
    ]

    if memory_limit_gb is not None:
        returncode, stderr = run_query_with_memory_limit(cmd, tmpfile, memory_limit_gb, log)
    else:
        proc = subprocess.run(cmd, stdout=tmpfile, stderr=subprocess.PIPE)
        returncode, stderr = proc.returncode, proc.stderr

    if returncode != 0:
        region_label = f"{region.chrom}:{region.start}-{region.end}"
        if returncode < 0:
            log.error(
                f"'query' was killed by signal {-returncode} "
                f"(likely OOM-killed{', exceeded ' + str(memory_limit_gb) + 'GB cap' if memory_limit_gb else ''}) "
                f"for {region_label}."
            )
        else:
            log.error(
                f"'query' failed (code {returncode}) for {region_label}. "
                f"stderr: {stderr.decode(errors='replace').strip()}"
            )
        if skipped_bed_file is not None:
            with open(skipped_bed_file, "a") as skipf:
                skipf.write(f"{region.chrom}\t{region.start}\t{region.end}\n")
        try:
            os.unlink(tmpfile.name)
        except OSError as cleanup_err:
            log.debug(f"Could not remove leftover temp file {tmpfile.name}: {cleanup_err}")
        return None

    return Path(tmpfile.name)
    
def check_gbzbase_installed(log: logging.Logger = None):
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


def index_gbz(gbz_file: Path, log: logging.Logger = None):
    """
    Index the GBZ file with gbz2db

    Parameters
    ----------
    gbz_file : Path
        Path to the GBZ file

    Returns
    -------
    passed : bool
        True if we were able to create the .gbz.db file
    """
    cmd = ["gbz2db", gbz_file]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE)
    return proc.returncode == 0


def check_gbzfile(gbz_file: Path, log: logging.Logger = None):
    """
    Check if the GBZ file exists and is
    indexed by GBZ-Base

    Parameters
    ----------
    gbz_file : Path
        Path to the GBZ file
    log : logging.Logger

    Returns
    -------
    passed : bool
        True if GBZ file and GBZ-base database exist
    """
    if not gbz_file.exists():
        log.critical(f"{gbz_file} does not exist\n")
        return False
    if not os.path.exists(str(gbz_file) + ".db"):
        log.info(f"{gbz_file}.db does not exist. Attempting to create")
        if not index_gbz(gbz_file):
            log.critical("Failed to create GBZ index")
            return False
    return True

def _set_mem_limit(max_bytes: int):
    resource.setrlimit(resource.RLIMIT_AS, (max_bytes, max_bytes))

def load_link_table_from_gbz(
    gbz_file: Path, region: Region, reference: str, exclude_samples=[], walk_file=None,
    log: logging.Logger = None, memory_limit_gb: Optional[float] = None,
    skipped_bed_file: Optional[Path] = None, pool=None,
) -> gutils.LinkTable:
    if log is None:
        log = _log
    gfa_file = extract_region_from_gbz(gbz_file, region, reference, memory_limit_gb=memory_limit_gb, log=log, skipped_bed_file=skipped_bed_file)
    if gfa_file is None:
        return gutils.LinkTable()

    region_label = f"{region.chrom}:{region.start}-{region.end}"
    return build_table_with_limit(
        gutils.LinkTable,
        dict(gfa_file=gfa_file, exclude_samples=exclude_samples),
        memory_limit_gb, log,
        region_label=f"{region_label} (link table)",
        skipped_bed_file=skipped_bed_file,
        bed_fields=(region.chrom, region.start, region.end),
        pool=pool,
    )