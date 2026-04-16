"""pdb_preprocessor.py — PDB chain/ligand extractor.

Parses a PDB file line-by-line (fixed-width format) and writes a filtered
copy that retains only the requested chains and ligands.  All blocking I/O
is offloaded to a thread-pool via ``asyncio.to_thread`` so the event loop
never blocks.

Typical usage for 7CFN — keep GPCR19 (chain R) + ligands FX0, CLR, PLM,
drop G-protein chains A/B/G/N and crystallographic waters (HOH):

    result = await extract_chains(
        pdb_path="7CFN.pdb",
        chains=["R"],
        keep_ligands=["FX0", "CLR", "PLM"],
    )
"""

import asyncio
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# PDB column offsets (0-based, following the PDB fixed-width spec)
# ---------------------------------------------------------------------------
#   cols  0- 5  record type  (ATOM / HETATM / ...)
#   cols  6-10  atom serial number
#   cols 11-11  space
#   cols 12-15  atom name
#   col  16     alternate location indicator
#   cols 17-19  residue name
#   col  20     space
#   col  21     chain ID
#   cols 22-25  residue sequence number
#   col  26     code for insertion of residues
#   cols 30-37  x orthogonal coordinate (Å)
#   cols 38-45  y orthogonal coordinate (Å)
#   cols 46-53  z orthogonal coordinate (Å)
#   cols 54-59  occupancy
#   cols 60-65  temperature factor
#   cols 76-77  element symbol
#   cols 78-79  charge

_ATOM_RECORDS = {"ATOM", "HETATM"}


def _pad_to_80(line: str) -> str:
    """Ensure a PDB line is exactly 80 chars (pad / trim trailing newline)."""
    stripped = line.rstrip("\n\r")
    if len(stripped) < 80:
        return stripped.ljust(80)
    return stripped[:80]


def _should_keep(
    record: str,
    chain: str,
    resname: str,
    keep_chains: set[str],
    keep_ligands: set[str] | None,
    remove_waters: bool,
) -> bool:
    """Return True if this ATOM/HETATM record should be included in output."""
    if record == "ATOM":
        # Standard protein/nucleic-acid atoms — keep if chain matches
        return chain in keep_chains

    # HETATM — ligand, solvent, or modified residue
    if remove_waters and resname.strip() == "HOH":
        return False

    # If the HETATM belongs to one of the target chains, keep it
    if chain in keep_chains:
        if keep_ligands is None:
            # None means "keep all non-water HETATMs in the chain"
            return True
        return resname.strip() in keep_ligands

    return False


def _renumber_serial(line: str, new_serial: int) -> str:
    """Replace the atom serial field (cols 6-10) with *new_serial*."""
    # Serial is right-justified in a 5-character field
    serial_str = f"{new_serial:5d}"
    return line[:6] + serial_str + line[11:]


def _run_extract_chains(
    pdb_path: Path,
    keep_chains: set[str],
    keep_ligands: set[str] | None,
    remove_waters: bool,
    output_path: Path,
) -> dict:
    """Blocking implementation of the chain/ligand extraction."""
    logger.info(
        "extract_chains: input=%s chains=%s ligands=%s output=%s",
        pdb_path,
        sorted(keep_chains),
        sorted(keep_ligands) if keep_ligands is not None else "all",
        output_path,
    )

    # ------------------------------------------------------------------
    # First pass: collect the set of atom serials we keep so that CONECT
    # records referencing only kept atoms are preserved.
    # ------------------------------------------------------------------
    kept_serials: set[int] = set()
    try:
        with pdb_path.open("r", errors="replace") as fh:
            for raw in fh:
                line = raw.rstrip("\n\r")
                record = line[:6].strip()
                if record not in _ATOM_RECORDS:
                    continue
                try:
                    chain = line[21]
                    resname = line[17:20]
                except IndexError:
                    continue
                if _should_keep(
                    record, chain, resname, keep_chains, keep_ligands, remove_waters
                ):
                    try:
                        kept_serials.add(int(line[6:11]))
                    except ValueError:
                        pass
    except OSError as exc:
        raise FileNotFoundError(f"Cannot open PDB: {pdb_path}") from exc

    # ------------------------------------------------------------------
    # Second pass: write filtered output with renumbered atom serials.
    # ------------------------------------------------------------------
    atom_count = 0
    hetatm_count = 0
    chains_found: set[str] = set()
    ligands_found: set[str] = set()

    # Map old serial → new serial for CONECT rewriting
    serial_map: dict[int, int] = {}
    new_serial = 0

    output_lines: list[str] = []

    try:
        with pdb_path.open("r", errors="replace") as fh:
            for raw in fh:
                line = raw.rstrip("\n\r")
                if not line:
                    continue

                record = line[:6].strip()

                if record in _ATOM_RECORDS:
                    try:
                        chain = line[21]
                        resname = line[17:20]
                    except IndexError:
                        continue

                    if not _should_keep(
                        record,
                        chain,
                        resname,
                        keep_chains,
                        keep_ligands,
                        remove_waters,
                    ):
                        continue

                    try:
                        old_serial = int(line[6:11])
                    except ValueError:
                        old_serial = None

                    new_serial += 1
                    if old_serial is not None:
                        serial_map[old_serial] = new_serial

                    padded = _pad_to_80(line)
                    padded = _renumber_serial(padded, new_serial)
                    output_lines.append(padded)

                    # Tally
                    if record == "ATOM":
                        atom_count += 1
                        chains_found.add(chain)
                    else:
                        hetatm_count += 1
                        ligands_found.add(resname.strip())

                elif record == "TER":
                    # TER may have a chain ID at col 21; only keep if in scope
                    try:
                        chain = line[21]
                        if chain in keep_chains:
                            output_lines.append(_pad_to_80(line))
                    except IndexError:
                        # Bare TER with no chain info — keep it
                        output_lines.append(_pad_to_80(line))

                elif record == "CONECT":
                    # CONECT cols 6-10, 11-15, 16-20, 21-25, 26-30 (1-indexed)
                    # i.e. 0-based: 6:11, 11:16, 16:21, 21:26, 26:31
                    try:
                        serials_in_line: list[int] = []
                        for start in (6, 11, 16, 21, 26):
                            field = line[start : start + 5].strip()
                            if field:
                                serials_in_line.append(int(field))
                        # Keep the CONECT record only if ALL referenced atoms
                        # are in the kept set (avoids dangling references)
                        if serials_in_line and all(
                            s in kept_serials for s in serials_in_line
                        ):
                            # Rewrite serial numbers
                            new_line = line[:6]
                            for start in (6, 11, 16, 21, 26):
                                field = line[start : start + 5].strip()
                                if field:
                                    orig = int(field)
                                    mapped = serial_map.get(orig, orig)
                                    new_line += f"{mapped:5d}"
                                else:
                                    new_line += "     "
                            output_lines.append(_pad_to_80(new_line))
                    except (ValueError, IndexError):
                        pass  # Malformed CONECT — skip silently

                elif record == "END":
                    output_lines.append(_pad_to_80(line))
                    break  # Stop after END

    except OSError as exc:
        raise FileNotFoundError(f"Cannot open PDB: {pdb_path}") from exc

    # Ensure file ends with END if it wasn't already present
    if not output_lines or output_lines[-1].strip() != "END":
        output_lines.append("END" + " " * 77)

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as out:
        out.write("\n".join(output_lines) + "\n")

    logger.info(
        "extract_chains: wrote %d ATOM + %d HETATM records to %s",
        atom_count,
        hetatm_count,
        output_path,
    )

    return {
        "success": True,
        "output_path": str(output_path),
        "atom_count": atom_count,
        "hetatm_count": hetatm_count,
        "chains_found": sorted(chains_found),
        "ligands_found": sorted(ligands_found),
    }


# ---------------------------------------------------------------------------
# Public async API
# ---------------------------------------------------------------------------


async def extract_chains(
    pdb_path: str,
    chains: list[str],
    keep_ligands: list[str] | None = None,
    remove_waters: bool = True,
    output_path: str | None = None,
) -> dict:
    """Extract specific chains and ligands from a PDB file.

    Reads the PDB file line-by-line using the fixed-width PDB format and
    writes a new PDB containing only the requested chains and (optionally)
    a whitelist of HETATM residue names.  Atom serial numbers are renumbered
    sequentially starting at 1 in the output.  CONECT records are preserved
    when all atoms they reference are retained; otherwise they are dropped.

    Parameters
    ----------
    pdb_path:
        Path to the source PDB file.
    chains:
        List of chain IDs to keep (e.g. ``["R"]``).  Case-sensitive.
    keep_ligands:
        Whitelist of HETATM residue names to keep (e.g. ``["FX0", "CLR",
        "PLM"]``).  ``None`` means keep *all* non-water HETATM records that
        belong to the selected chains.
    remove_waters:
        Drop HOH (water) HETATM records even when they are in the selected
        chains.  Defaults to ``True``.
    output_path:
        Destination file path.  When ``None`` the output is written next to
        the source PDB with ``_extracted`` appended to the stem, e.g.
        ``7CFN_extracted.pdb``.

    Returns
    -------
    dict
        ``success`` (bool), ``output_path`` (str), ``atom_count`` (int),
        ``hetatm_count`` (int), ``chains_found`` (list[str]),
        ``ligands_found`` (list[str]).
    """
    src = Path(pdb_path).resolve()
    if not src.is_file():
        raise FileNotFoundError(f"PDB file not found: {src}")

    if output_path is None:
        dest = src.with_name(f"{src.stem}_extracted{src.suffix}")
    else:
        dest = Path(output_path).resolve()

    keep_chains: set[str] = set(chains)
    keep_lig_set: set[str] | None = (
        set(keep_ligands) if keep_ligands is not None else None
    )

    return await asyncio.to_thread(
        _run_extract_chains,
        src,
        keep_chains,
        keep_lig_set,
        remove_waters,
        dest,
    )
