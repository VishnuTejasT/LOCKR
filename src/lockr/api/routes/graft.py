"""GET /graft/status, POST /graft, GET /graft/download/{job_id}.

The threading scan takes 1-3 minutes (a real PyRosetta run over ~17
positions), so /graft runs it off the event loop via run_in_threadpool
rather than blocking the whole API while it works.

Grafted PDBs are too large for the JSON response, so they're written to a
temp file server-side, keyed by job_id, and cleaned up lazily on next access
past their expiry rather than on a background timer, this is a single-user
local tool, not a hosted service that needs a real job queue.
"""

from __future__ import annotations

import os
import time
import uuid
from importlib.metadata import PackageNotFoundError, version as pkg_version
from importlib.resources import files

from fastapi import APIRouter
from fastapi.concurrency import run_in_threadpool
from fastapi.responses import FileResponse, JSONResponse

from lockr.engine import graft

from ..errors import ApiError
from ..schemas.graft import GraftRequest, GraftResponse, GraftStatusResponse, ScoredPosition

router = APIRouter()

_PYROSETTA_INSTALL_URL = "https://west.rosettacommons.org/pyrosetta/quarterly/release"

_JOB_EXPIRY_SECONDS = 60 * 60
# job_id -> (pdb_path, expiry_unix_time)
_jobs: dict[str, tuple[str, float]] = {}

_TEMPLATE_NAME = "lucCage_template_clean.pdb"


def _template_path() -> str:
    return str(files("lockr.data").joinpath(_TEMPLATE_NAME))


def _template_bundled() -> bool:
    return os.path.exists(_template_path())


def _pyrosetta_version() -> str | None:
    if not graft.PYROSETTA_AVAILABLE:
        return None
    try:
        return pkg_version("pyrosetta")
    except PackageNotFoundError:
        return None


def _sweep_expired_jobs() -> None:
    now = time.time()
    expired = [job_id for job_id, (_, expiry) in _jobs.items() if expiry < now]
    for job_id in expired:
        path, _ = _jobs.pop(job_id)
        if os.path.exists(path):
            os.remove(path)


def _store_job(pdb_path: str) -> str:
    _sweep_expired_jobs()
    job_id = str(uuid.uuid4())
    _jobs[job_id] = (pdb_path, time.time() + _JOB_EXPIRY_SECONDS)
    return job_id


@router.get("/graft/status", response_model=GraftStatusResponse)
def graft_status() -> GraftStatusResponse:
    return GraftStatusResponse(
        available=graft.PYROSETTA_AVAILABLE,
        version=_pyrosetta_version(),
        template_bundled=_template_bundled(),
        calibration_warning=graft._calibration_mismatch_warning() if graft.PYROSETTA_AVAILABLE else None,
    )


@router.post("/graft", response_model=GraftResponse)
async def run_graft(request: GraftRequest) -> GraftResponse:
    if not graft.PYROSETTA_AVAILABLE:
        return JSONResponse(status_code=503, content={"error": {
            "code": "PYROSETTA_UNAVAILABLE",
            "message": "PyRosetta is not installed. See INSTALL.md for setup.",
            "install_guide": _PYROSETTA_INSTALL_URL,
        }})

    template_path = _template_path()

    try:
        if request.scan_all or request.specific_position is None:
            result = await run_in_threadpool(
                graft.find_best_graft_position,
                request.sequence, template_path,
                latch_start=request.latch_start, latch_end=request.latch_end,
            )
            all_scores = [ScoredPosition(position=p, score=s) for p, s in result.all_scores]
            job_id = _store_job(result.grafted_pdb_path)
            return GraftResponse(
                best_position=result.best_position, best_score=result.best_score,
                verdict=result.verdict, all_scores=all_scores,
                grafted_sequence=result.grafted_sequence, job_id=job_id,
                runtime_seconds=result.runtime_seconds, binder_length=result.binder_length,
                calibration_warning=result.calibration_warning,
            )
        else:
            start = time.time()
            at_result = await run_in_threadpool(
                graft.graft_at_position, request.sequence, template_path, request.specific_position,
            )
            job_id = _store_job(at_result.grafted_pdb_path)
            return GraftResponse(
                best_position=at_result.position, best_score=at_result.score,
                verdict=at_result.verdict, all_scores=[ScoredPosition(position=at_result.position, score=at_result.score)],
                grafted_sequence=at_result.grafted_sequence, job_id=job_id,
                runtime_seconds=time.time() - start, binder_length=len(request.sequence),
                calibration_warning=at_result.calibration_warning,
            )
    except ValueError as e:
        raise ApiError("NO_VALID_POSITIONS", str(e))


@router.get("/graft/download/{job_id}")
def download_graft(job_id: str):
    _sweep_expired_jobs()
    if job_id not in _jobs:
        raise ApiError("JOB_NOT_FOUND", "This grafted PDB has expired or never existed, run the graft again.",
                       status_code=404)
    path, _ = _jobs[job_id]
    return FileResponse(path, media_type="chemical/x-pdb", filename=f"lucCage_grafted_{job_id[:8]}.pdb")
