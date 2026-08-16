import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import NamedTuple

from rdkit import Chem


_DIMORPHITE_BRIDGE = Path(__file__).with_name("dimorphite_bridge.py")


class ParentFormReference(NamedTuple):
    smiles: str
    backend_requested: str
    backend_used: str
    candidates: list[str]
    comment: str


def _canonicalize_smiles(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _dedupe_smiles(values: list[str]) -> list[str]:
    unique_values = []
    seen = set()
    for value in values:
        if value not in seen:
            seen.add(value)
            unique_values.append(value)
    return unique_values


def _fallback_parent_form_reference(
    parent_smiles: str,
    backend_requested: str,
    comment: str,
) -> ParentFormReference:
    canonical_parent = _canonicalize_smiles(parent_smiles) or parent_smiles
    return ParentFormReference(
        smiles=canonical_parent,
        backend_requested=backend_requested,
        backend_used="fallback",
        candidates=[canonical_parent],
        comment=comment,
    )


def _dimorphite_reference_from_candidates(
    parent_smiles: str,
    candidates: list[str],
    comment: str,
) -> ParentFormReference:
    canonical_parent = _canonicalize_smiles(parent_smiles) or parent_smiles
    canonical_candidates = []
    for candidate in candidates:
        canonical_candidate = _canonicalize_smiles(candidate)
        if canonical_candidate is not None:
            canonical_candidates.append(canonical_candidate)

    canonical_candidates = _dedupe_smiles(canonical_candidates)
    if not canonical_candidates:
        canonical_candidates = [canonical_parent]

    return ParentFormReference(
        smiles=canonical_candidates[0],
        backend_requested="dimorphite_dl",
        backend_used="dimorphite_dl",
        candidates=canonical_candidates,
        comment=comment,
    )


def _conda_env_python(env_name: str | None) -> str | None:
    if not env_name:
        return None

    conda_exe = os.environ.get("CONDA_EXE") or shutil.which("conda")
    if not conda_exe:
        return None

    conda_root = Path(conda_exe).resolve().parent.parent
    if os.name == "nt":
        candidate = conda_root / "envs" / env_name / "python.exe"
    else:
        candidate = conda_root / "envs" / env_name / "bin" / "python"

    if candidate.exists():
        return str(candidate)
    return None


def _resolve_dimorphite_command(
    dimorphite_python: str | None = None,
    dimorphite_conda_env: str | None = "dimorphite",
) -> tuple[list[str] | None, str]:
    explicit_python = dimorphite_python or os.environ.get("MOLJAM_DIMORPHITE_PYTHON")
    if explicit_python:
        return [explicit_python], f"external python '{explicit_python}'"

    env_python = _conda_env_python(dimorphite_conda_env)
    if env_python:
        return [env_python], f"conda env python '{env_python}'"

    conda_exe = os.environ.get("CONDA_EXE") or shutil.which("conda")
    if conda_exe and dimorphite_conda_env:
        return (
            [conda_exe, "run", "--no-capture-output", "-n", dimorphite_conda_env, "python"],
            f"conda run env '{dimorphite_conda_env}'",
        )

    return None, "no external Dimorphite-DL command configured"


def _run_dimorphite_bridge(
    parent_smiles_list: list[str],
    ph: float,
    dimorphite_python: str | None = None,
    dimorphite_conda_env: str | None = "dimorphite",
) -> tuple[dict[str, list[str]] | None, str]:
    command, command_comment = _resolve_dimorphite_command(
        dimorphite_python=dimorphite_python,
        dimorphite_conda_env=dimorphite_conda_env,
    )
    if command is None:
        return None, command_comment

    payload = {
        "ph": float(ph),
        "smiles": list(parent_smiles_list),
    }
    try:
        completed = subprocess.run(
            [*command, str(_DIMORPHITE_BRIDGE)],
            check=False,
            capture_output=True,
            text=True,
            input=json.dumps(payload),
            timeout=120,
        )
    except Exception as exc:  # pragma: no cover - defensive fallback
        return None, f"Dimorphite-DL bridge failed via {command_comment} ({exc})"

    if completed.returncode != 0:
        stderr = completed.stderr.strip() or "unknown error"
        return None, f"Dimorphite-DL bridge failed via {command_comment} ({stderr})"

    try:
        raw_results = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:  # pragma: no cover - defensive fallback
        return None, f"Dimorphite-DL bridge returned invalid JSON via {command_comment} ({exc})"

    results = {}
    for parent_smiles in parent_smiles_list:
        value = raw_results.get(parent_smiles, [])
        if isinstance(value, list):
            results[parent_smiles] = value
        else:
            results[parent_smiles] = []
    return results, f"Parent form derived with Dimorphite-DL at pH {ph:.1f} via {command_comment}"


def _run_local_dimorphite(
    parent_smiles_list: list[str],
    ph: float,
) -> tuple[dict[str, list[str]] | None, str]:
    try:
        from dimorphite_dl import protonate_smiles as local_protonate_smiles
    except ImportError:
        return None, "Dimorphite-DL unavailable in the current environment"

    results = {}
    try:
        for parent_smiles in parent_smiles_list:
            results[parent_smiles] = local_protonate_smiles(
                parent_smiles,
                ph_min=ph,
                ph_max=ph,
                precision=0.0,
                max_variants=32,
            )
    except Exception as exc:  # pragma: no cover - defensive fallback
        return None, f"Local Dimorphite-DL failed ({exc})"

    return results, f"Parent form derived with Dimorphite-DL at pH {ph:.1f} in-process"


def _select_dimorphite_parent_forms(
    parent_smiles_list: list[str],
    ph: float,
    dimorphite_python: str | None = None,
    dimorphite_conda_env: str | None = "dimorphite",
) -> dict[str, ParentFormReference]:
    unique_parent_smiles = _dedupe_smiles(parent_smiles_list)
    if not unique_parent_smiles:
        return {}

    batch_results, batch_comment = _run_dimorphite_bridge(
        unique_parent_smiles,
        ph,
        dimorphite_python=dimorphite_python,
        dimorphite_conda_env=dimorphite_conda_env,
    )
    if batch_results is None:
        batch_results, batch_comment = _run_local_dimorphite(unique_parent_smiles, ph)

    references = {}
    for parent_smiles in unique_parent_smiles:
        if batch_results is None:
            references[parent_smiles] = _fallback_parent_form_reference(
                parent_smiles,
                "dimorphite_dl",
                f"{batch_comment}; using structural parent as parent form",
            )
            continue

        references[parent_smiles] = _dimorphite_reference_from_candidates(
            parent_smiles,
            batch_results.get(parent_smiles, []),
            batch_comment,
        )

    return references


def _select_chemaxon_parent_form(
    parent_smiles: str,
    ph: float,
    executable: str,
) -> ParentFormReference:
    canonical_parent = _canonicalize_smiles(parent_smiles) or parent_smiles
    if shutil.which(executable) is None:
        return ParentFormReference(
            smiles=canonical_parent,
            backend_requested="chemaxon",
            backend_used="fallback",
            candidates=[canonical_parent],
            comment=f"ChemAxon executable '{executable}' not found; using structural parent as parent form",
        )

    command = [
        executable,
        "-N",
        "ih",
        canonical_parent,
        "majormicrospecies",
        "-H",
        f"{ph:.1f}",
        "-f",
        "smiles",
    ]
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except Exception as exc:  # pragma: no cover - defensive fallback
        return ParentFormReference(
            smiles=canonical_parent,
            backend_requested="chemaxon",
            backend_used="fallback",
            candidates=[canonical_parent],
            comment=f"ChemAxon execution failed ({exc}); using structural parent as parent form",
        )

    if completed.returncode != 0:
        stderr = completed.stderr.strip() or "unknown error"
        return ParentFormReference(
            smiles=canonical_parent,
            backend_requested="chemaxon",
            backend_used="fallback",
            candidates=[canonical_parent],
            comment=f"ChemAxon majormicrospecies failed ({stderr}); using structural parent as parent form",
        )

    output_lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    canonical_candidates = []
    for line in output_lines:
        first_field = line.split("\t", 1)[0].strip()
        canonical_candidate = _canonicalize_smiles(first_field)
        if canonical_candidate is not None:
            canonical_candidates.append(canonical_candidate)

    canonical_candidates = _dedupe_smiles(canonical_candidates)
    if not canonical_candidates:
        canonical_candidates = [canonical_parent]

    return ParentFormReference(
        smiles=canonical_candidates[0],
        backend_requested="chemaxon",
        backend_used="chemaxon",
        candidates=canonical_candidates,
        comment=f"Parent form derived with ChemAxon at pH {ph:.1f}",
    )


def derive_parent_form_references(
    parent_smiles_list: list[str],
    backend: str = "dimorphite_dl",
    ph: float = 7.4,
    chemaxon_executable: str = "cxcalc",
    dimorphite_python: str | None = None,
    dimorphite_conda_env: str | None = "dimorphite",
) -> dict[str, ParentFormReference]:
    backend_key = (backend or "dimorphite_dl").strip().lower()
    unique_parent_smiles = _dedupe_smiles(parent_smiles_list)

    if backend_key == "chemaxon":
        return {
            parent_smiles: _select_chemaxon_parent_form(parent_smiles, ph, chemaxon_executable)
            for parent_smiles in unique_parent_smiles
        }
    if backend_key == "dimorphite_dl":
        return _select_dimorphite_parent_forms(
            unique_parent_smiles,
            ph,
            dimorphite_python=dimorphite_python,
            dimorphite_conda_env=dimorphite_conda_env,
        )

    return {
        parent_smiles: _fallback_parent_form_reference(
            parent_smiles,
            backend_key,
            f"Unknown parent form backend '{backend}'; using structural parent as parent form",
        )
        for parent_smiles in unique_parent_smiles
    }


def derive_parent_form_reference(
    parent_smiles: str,
    backend: str = "dimorphite_dl",
    ph: float = 7.4,
    chemaxon_executable: str = "cxcalc",
    dimorphite_python: str | None = None,
    dimorphite_conda_env: str | None = "dimorphite",
) -> ParentFormReference:
    references = derive_parent_form_references(
        [parent_smiles],
        backend=backend,
        ph=ph,
        chemaxon_executable=chemaxon_executable,
        dimorphite_python=dimorphite_python,
        dimorphite_conda_env=dimorphite_conda_env,
    )
    return references[parent_smiles]
