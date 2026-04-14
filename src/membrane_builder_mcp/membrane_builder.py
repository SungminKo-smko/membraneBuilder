"""membrane_builder.py — packmol-memgen CLI wrapper.

Wraps AmberTools25's ``packmol-memgen`` command so it can be called from the
MCP server as plain async Python functions.  All heavy work is offloaded to a
thread-pool via ``asyncio.to_thread`` so the event loop never blocks.
"""

import asyncio
import logging
import os
import re
import shutil
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Conda / environment configuration
# ---------------------------------------------------------------------------

CONDA_ENV = os.environ.get("MEMBRANE_CONDA_ENV", "AmberTools25")
CONDA_BASE = Path(os.environ.get("CONDA_BASE", str(Path.home() / "miniconda3")))


ENV_BIN = CONDA_BASE / "envs" / CONDA_ENV / "bin"


def _get_env() -> dict[str, str]:
    """Return an environment dict with the conda env's bin/ prepended to PATH.

    This avoids ``conda run`` which sometimes fails to propagate PATH correctly,
    causing tools like ``reduce`` or ``packmol`` to not be found.
    """
    env = os.environ.copy()
    env["PATH"] = f"{ENV_BIN}:{env.get('PATH', '')}"
    # AmberTools programs also look for AMBERHOME
    amber_home = CONDA_BASE / "envs" / CONDA_ENV
    env["AMBERHOME"] = str(amber_home)
    return env


def _get_executable(name: str) -> str:
    """Return the absolute path to a binary inside the conda env."""
    exe = ENV_BIN / name
    if exe.is_file() and os.access(exe, os.X_OK):
        return str(exe)
    raise FileNotFoundError(f"{name} not found in {ENV_BIN}")


# ---------------------------------------------------------------------------
# list_available_lipids
# ---------------------------------------------------------------------------

def _parse_lipid_line(line: str) -> dict | None:
    """Parse a single output line from ``packmol-memgen --available_lipids_all``.

    Expected formats (examples from the tool):
        POPC     charge:  0   1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine
        DOPE     charge: -1   ...
    Returns None when the line does not match.
    """
    # Strip ANSI escape codes that some terminal-aware programs emit
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    clean = ansi_escape.sub("", line).strip()

    # Actual format: NAME    CHARGE    FULL_NAME    [COMMENT]
    # e.g.: AHPA         -1               1-arachidonoyl-2-...
    m = re.match(
        r"^(?P<name>[A-Z][A-Z0-9]{2,5})\s+"
        r"(?P<charge>-?\d+)\s+"
        r"(?P<rest>.+)?$",
        clean,
    )
    if not m:
        return None

    rest = (m.group("rest") or "").strip()
    parts = rest.rsplit("  ", 1)
    full_name = parts[0].strip() if parts else rest
    comment = parts[1].strip() if len(parts) > 1 else ""
    return {
        "name": m.group("name"),
        "charge": int(m.group("charge")),
        "full_name": full_name,
        "comment": comment,
    }


def _run_list_lipids() -> list[dict]:
    """Blocking helper — run ``packmol-memgen --available_lipids_all``."""
    cmd = [_get_executable("packmol-memgen"), "--available_lipids_all"]
    logger.debug("list_available_lipids cmd: %s", " ".join(cmd))

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
            env=_get_env(),
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"packmol-memgen not found in {ENV_BIN}. "
            "Make sure miniconda3/anaconda3 is installed."
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError("packmol-memgen --available_lipids_all timed out.") from exc

    output = result.stdout + result.stderr
    lipids: list[dict] = []
    for line in output.splitlines():
        parsed = _parse_lipid_line(line)
        if parsed:
            lipids.append(parsed)

    return lipids


async def list_available_lipids() -> list[dict]:
    """Return all lipids supported by ``packmol-memgen``.

    Each entry is a dict::

        {
            "name":      str,   # e.g. "POPC"
            "charge":    int,   # net charge (0, -1, …)
            "full_name": str,   # long name from the tool output
            "comment":   str,   # additional comment (may be empty)
        }
    """
    return await asyncio.to_thread(_run_list_lipids)


# ---------------------------------------------------------------------------
# build_membrane
# ---------------------------------------------------------------------------

def _run_build_membrane(
    *,
    output_dir: Path,
    pdb_filename: str,
    cmd: list[str],
    timeout: int,
) -> dict:
    """Blocking helper that actually executes ``packmol-memgen``."""
    logger.info("build_membrane cwd=%s cmd=%s", output_dir, " ".join(cmd))

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(output_dir),
            capture_output=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=timeout,
            env=_get_env(),
        )
        full_log = proc.stdout or ""
        returncode = proc.returncode
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"packmol-memgen not found in {ENV_BIN}. "
            f"Searched in {CONDA_BASE}. "
            "Ensure AmberTools25 conda env exists."
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"packmol-memgen exceeded the {timeout}s timeout."
        ) from exc

    # ------------------------------------------------------------------
    # Collect output files
    # ------------------------------------------------------------------
    stem = Path(pdb_filename).stem  # filename without extension

    # packmol-memgen names outputs like <stem>_lipid.pdb / .top / .crd
    # but sometimes just <stem>.pdb when no protein is provided.
    # Search flexibly.
    pdb_candidates = sorted(output_dir.glob("*_lipid*.pdb")) + sorted(
        output_dir.glob(f"{stem}*.pdb")
    )
    top_candidates = sorted(output_dir.glob("*.top")) + sorted(
        output_dir.glob("*.prmtop")
    )
    crd_candidates = (
        sorted(output_dir.glob("*.crd"))
        + sorted(output_dir.glob("*.rst7"))
        + sorted(output_dir.glob("*.inpcrd"))
    )

    output_pdb = str(pdb_candidates[0]) if pdb_candidates else None
    output_top = str(top_candidates[0]) if top_candidates else None
    output_crd = str(crd_candidates[0]) if crd_candidates else None

    # Last 30 lines of the log for the caller
    log_tail = "\n".join(full_log.splitlines()[-30:])

    return {
        "success": returncode == 0,
        "returncode": returncode,
        "output_pdb": output_pdb,
        "output_top": output_top,
        "output_crd": output_crd,
        "log": log_tail,
        "command": " ".join(cmd),
    }


async def build_membrane(
    protein_pdb_path: str,
    output_dir: str = "./output",
    lipids: str = "POPC",
    ratio: str = "1",
    distxy_fix: float | None = None,
    dist: float = 15.0,
    dist_wat: float = 17.5,
    salt: bool = True,
    salt_concentration: float = 0.15,
    parametrize: bool = True,
    ffprot: str = "ff19SB",
    ffwat: str = "tip3p",
    fflip: str = "lipid21",
    preoriented: bool = False,
    keep_ligands: bool = True,
    keep_files: bool = True,
    timeout: int = 7200,
) -> dict:
    """Run ``packmol-memgen`` and return paths to the generated files.

    Parameters
    ----------
    protein_pdb_path:
        Absolute or relative path to the input protein PDB file.
    output_dir:
        Directory where output files will be placed.  Created if absent.
    lipids:
        Colon-separated lipid names, e.g. ``"POPC"`` or ``"DOPE:DOPG"``.
    ratio:
        Colon-separated molar ratios matching *lipids*, e.g. ``"3:1"``.
    distxy_fix:
        Fixed XY box size in Å.  ``None`` → packmol-memgen chooses automatically.
    dist:
        Minimum distance (Å) between the protein and the membrane edge.
    dist_wat:
        Water layer thickness in Å (default 17.5).
    salt:
        Whether to add KCl ions.
    salt_concentration:
        Ion concentration in mol/L (default 0.15 M).
    parametrize:
        Run ``tleap`` to generate AMBER topology/coordinate files.
    ffprot:
        Protein force field (``"ff14SB"`` or ``"ff19SB"``).
    ffwat:
        Water model (``"tip3p"``, ``"opc"``, ``"spce"``).
    fflip:
        Lipid force field (``"lipid21"`` or ``"lipid17"``).
    preoriented:
        Protein is already oriented along the membrane normal — skip OPM lookup.
    keep_ligands:
        Pass ``--keepligs`` to preserve HETATM ligand records.
    keep_files:
        Pass ``--keep`` to retain intermediate files.
    timeout:
        Wall-clock timeout in seconds (default 3600 = 1 hour).

    Returns
    -------
    dict
        ``success`` (bool), ``output_pdb``, ``output_top``, ``output_crd``,
        ``log`` (last 30 lines), ``command`` (the shell command executed).
    """
    # ------------------------------------------------------------------
    # 1. Prepare output directory
    # ------------------------------------------------------------------
    out_dir = Path(output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 2. Copy protein PDB into output_dir
    #    packmol-memgen writes outputs in the *current working directory*,
    #    so we run it with cwd=output_dir and hand it a bare filename.
    # ------------------------------------------------------------------
    src_pdb = Path(protein_pdb_path).resolve()
    if not src_pdb.is_file():
        raise FileNotFoundError(f"Protein PDB not found: {src_pdb}")

    dest_pdb = out_dir / src_pdb.name
    if src_pdb != dest_pdb:
        shutil.copy2(src_pdb, dest_pdb)
    pdb_filename = src_pdb.name  # relative, used as the --pdb argument

    # ------------------------------------------------------------------
    # 3. Build the command
    # ------------------------------------------------------------------
    cmd: list[str] = [_get_executable("packmol-memgen")]
    cmd += ["--pdb", pdb_filename]
    cmd += ["--lipids", lipids]
    cmd += ["--ratio", ratio]

    if distxy_fix is not None:
        cmd += ["--distxy_fix", str(distxy_fix)]

    cmd += ["--dist", str(dist)]
    cmd += ["--dist_wat", str(dist_wat)]

    if salt:
        cmd += ["--salt", "--saltcon", str(salt_concentration)]

    if parametrize:
        cmd += ["--parametrize"]

    cmd += ["--ffprot", ffprot]
    cmd += ["--ffwat", ffwat]
    cmd += ["--fflip", fflip]

    if preoriented:
        cmd += ["--preoriented"]
    if keep_ligands:
        cmd += ["--keepligs"]
    if keep_files:
        cmd += ["--keep"]

    cmd += ["--overwrite"]

    # ------------------------------------------------------------------
    # 4. Execute (non-blocking)
    # ------------------------------------------------------------------
    return await asyncio.to_thread(
        _run_build_membrane,
        output_dir=out_dir,
        pdb_filename=pdb_filename,
        cmd=cmd,
        timeout=timeout,
    )


# ---------------------------------------------------------------------------
# analyze_protein
# ---------------------------------------------------------------------------

def _parse_pdb_atoms(pdb_path: Path) -> list[dict]:
    """Return a list of atom dicts from ATOM/HETATM lines in a PDB file."""
    atoms: list[dict] = []
    try:
        with pdb_path.open("r", errors="replace") as fh:
            for line in fh:
                record = line[:6].strip()
                if record not in ("ATOM", "HETATM"):
                    continue
                try:
                    atoms.append(
                        {
                            "record": record,
                            "serial": int(line[6:11]),
                            "name": line[12:16].strip(),
                            "resname": line[17:20].strip(),
                            "chain": line[21].strip(),
                            "resseq": int(line[22:26]),
                            "x": float(line[30:38]),
                            "y": float(line[38:46]),
                            "z": float(line[46:54]),
                        }
                    )
                except (ValueError, IndexError):
                    continue
    except OSError as exc:
        raise FileNotFoundError(f"Cannot open PDB: {pdb_path}") from exc
    return atoms


def _run_analyze_protein(protein_pdb_path: str) -> dict:
    """Blocking analysis of a PDB file."""
    pdb_path = Path(protein_pdb_path).resolve()
    if not pdb_path.is_file():
        raise FileNotFoundError(f"Protein PDB not found: {pdb_path}")

    atoms = _parse_pdb_atoms(pdb_path)

    if not atoms:
        return {
            "file": str(pdb_path),
            "atom_count": 0,
            "hetatm_count": 0,
            "residue_count": 0,
            "chains": [],
            "bounding_box": None,
            "error": "No ATOM/HETATM records found.",
        }

    # Counts
    atom_records = [a for a in atoms if a["record"] == "ATOM"]
    hetatm_records = [a for a in atoms if a["record"] == "HETATM"]

    # Unique residues: (chain, resseq, resname) tuples
    residues: set[tuple] = set()
    for a in atom_records:
        residues.add((a["chain"], a["resseq"], a["resname"]))

    # Chains (from ATOM records only, preserving insertion order)
    chains: list[str] = []
    seen_chains: set[str] = set()
    for a in atom_records:
        c = a["chain"]
        if c not in seen_chains:
            chains.append(c)
            seen_chains.add(c)

    # Bounding box from all atoms
    xs = [a["x"] for a in atoms]
    ys = [a["y"] for a in atoms]
    zs = [a["z"] for a in atoms]
    bounding_box = {
        "x_min": min(xs), "x_max": max(xs),
        "y_min": min(ys), "y_max": max(ys),
        "z_min": min(zs), "z_max": max(zs),
        "x_size": max(xs) - min(xs),
        "y_size": max(ys) - min(ys),
        "z_size": max(zs) - min(zs),
    }

    return {
        "file": str(pdb_path),
        "atom_count": len(atom_records),
        "hetatm_count": len(hetatm_records),
        "residue_count": len(residues),
        "chains": chains,
        "bounding_box": bounding_box,
    }


async def analyze_protein(protein_pdb_path: str) -> dict:
    """Analyze a PDB file and return basic structural information.

    Returns
    -------
    dict
        ``file``, ``atom_count``, ``hetatm_count``, ``residue_count``,
        ``chains`` (list of chain IDs from ATOM records),
        ``bounding_box`` (min/max/size for x, y, z in Å).
    """
    return await asyncio.to_thread(_run_analyze_protein, protein_pdb_path)
