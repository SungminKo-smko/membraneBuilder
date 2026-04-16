"""openmm_runner.py — OpenMM MD simulation runner for membrane protein systems.

Runs multi-stage MD simulations (energy minimization → NVT equilibration →
NPT equilibration → production) using OpenMM with GROMACS input files.
Designed for GPCR-in-lipid-bilayer systems; uses MonteCarloMembraneBarostat
for semiisotropic pressure coupling during NPT and production stages.

All heavy work is offloaded to a thread-pool via ``asyncio.to_thread`` so
the event loop never blocks.
"""

import asyncio
import logging
import os
from pathlib import Path

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Optional OpenMM import — fail gracefully if not installed
# ---------------------------------------------------------------------------

try:
    import openmm
    from openmm import (
        LangevinMiddleIntegrator,
        MonteCarloMembraneBarostat,
        CustomExternalForce,
        Platform,
        unit,
    )
    from openmm.app import (
        DCDReporter,
        GromacsGroFile,
        GromacsTopFile,
        HBonds,
        PME,
        Simulation,
        StateDataReporter,
    )
    from openmm.unit import (
        bar,
        femtosecond,
        kelvin,
        kilojoule_per_mole,
        nanometer,
        picosecond,
    )

    _OPENMM_AVAILABLE = True
except ImportError:
    _OPENMM_AVAILABLE = False

# ---------------------------------------------------------------------------
# Water / ion residue names to skip when adding position restraints
# ---------------------------------------------------------------------------

_WATER_ION_RESIDUES: frozenset[str] = frozenset(
    {
        # Water models
        "HOH", "WAT", "TIP3", "TIP4", "TIP5", "SOL", "SPC",
        # Common ions
        "NA", "CL", "K", "MG", "CA", "ZN", "FE", "MN",
        "Na+", "Cl-", "K+",
        # GROMACS-style ion names
        "SOD", "CLA", "POT",
    }
)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _select_platform(requested: str, precision: str) -> "Platform":
    """Return an OpenMM Platform, falling back to CPU if *requested* is unavailable."""
    properties: dict[str, str] = {}

    # Try the requested platform first
    for platform_name in (requested, "OpenCL", "CPU"):
        try:
            platform = Platform.getPlatformByName(platform_name)
            if platform_name in ("CUDA", "OpenCL"):
                properties = {"Precision": precision}
            logger.info("Using OpenMM platform: %s (precision=%s)", platform_name, precision if properties else "default")
            return platform, properties
        except Exception:
            logger.warning("Platform %s not available, trying next.", platform_name)

    raise RuntimeError("No usable OpenMM platform found (tried CUDA, OpenCL, CPU).")


def _add_position_restraints(
    system: "openmm.System",
    topology: "openmm.app.Topology",
    positions,
    force_constant: float = 1000.0,
) -> "CustomExternalForce":
    """Add harmonic position restraints to all non-hydrogen protein heavy atoms.

    Skips atoms belonging to water/ion residues.  Returns the
    ``CustomExternalForce`` object so the caller can later zero the force
    constant to release restraints.

    Parameters
    ----------
    force_constant:
        Spring constant in kJ/mol/nm².
    """
    restraint = CustomExternalForce(
        "k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
    )
    restraint.addGlobalParameter("k", force_constant * kilojoule_per_mole / nanometer**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    restrained_count = 0
    for atom in topology.atoms():
        # Skip water and ions
        if atom.residue.name in _WATER_ION_RESIDUES:
            continue
        # Skip hydrogen atoms
        if atom.element is not None and atom.element.symbol == "H":
            continue
        pos = positions[atom.index]
        restraint.addParticle(
            atom.index,
            [pos.x, pos.y, pos.z],
        )
        restrained_count += 1

    system.addForce(restraint)
    logger.debug("Added position restraints to %d heavy atoms.", restrained_count)
    return restraint


def _run_simulation(
    *,
    gro_path: str,
    top_path: str,
    output_dir: Path,
    minimize: bool,
    min_max_iterations: int,
    nvt_steps: int,
    npt_steps: int,
    production_steps: int,
    temperature: float,
    pressure: float,
    timestep: float,
    nonbonded_cutoff: float,
    report_interval: int,
    checkpoint_interval: int,
    platform: str,
    precision: str,
    gromacs_include_dir: str | None,
) -> dict:
    """Blocking simulation runner — called via asyncio.to_thread."""
    if not _OPENMM_AVAILABLE:
        raise ImportError(
            "OpenMM is not installed in the current Python environment. "
            "Install it with:\n"
            "  conda install -n AmberTools25 -c conda-forge openmm\n"
            "or activate the AmberTools25 environment before running."
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    gro_file = Path(gro_path).resolve()
    top_file = Path(top_path).resolve()
    if not gro_file.is_file():
        raise FileNotFoundError(f".gro file not found: {gro_file}")
    if not top_file.is_file():
        raise FileNotFoundError(f".top file not found: {top_file}")

    # ------------------------------------------------------------------
    # 1. Load GROMACS files
    # ------------------------------------------------------------------
    logger.info("Loading GROMACS files: gro=%s  top=%s", gro_file, top_file)

    gro = GromacsGroFile(str(gro_file))

    top_kwargs: dict = {"periodicBoxVectors": gro.getPeriodicBoxVectors()}
    if gromacs_include_dir:
        top_kwargs["includeDir"] = gromacs_include_dir

    top = GromacsTopFile(str(top_file), **top_kwargs)

    # ------------------------------------------------------------------
    # 2. Create system
    # ------------------------------------------------------------------
    logger.info("Creating system (PME, cutoff=%.1f nm, HBonds constraints).", nonbonded_cutoff)
    system = top.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=nonbonded_cutoff * nanometer,
        constraints=HBonds,
    )

    # ------------------------------------------------------------------
    # 3. Select platform
    # ------------------------------------------------------------------
    openmm_platform, platform_properties = _select_platform(platform, precision)

    # ------------------------------------------------------------------
    # 4. Integrator
    # ------------------------------------------------------------------
    integrator = LangevinMiddleIntegrator(
        temperature * kelvin,
        1.0 / picosecond,          # friction coefficient
        timestep * picosecond,
    )

    # ------------------------------------------------------------------
    # 5. Build simulation (no barostat yet — used for minimization + NVT)
    # ------------------------------------------------------------------
    simulation = Simulation(
        top.topology,
        system,
        integrator,
        openmm_platform,
        platform_properties,
    )
    simulation.context.setPositions(gro.positions)
    simulation.context.setVelocitiesToTemperature(temperature * kelvin)

    result: dict = {
        "success": False,
        "stage": "init",
        "output_files": {},
        "minimized_energy": None,
        "final_temperature": None,
        "log": "",
    }
    log_lines: list[str] = []

    def _log(msg: str) -> None:
        logger.info(msg)
        log_lines.append(msg)

    # ------------------------------------------------------------------
    # 6. Energy minimization
    # ------------------------------------------------------------------
    if minimize:
        _log(f"[minimization] Starting energy minimization (max_iterations={min_max_iterations}) ...")
        result["stage"] = "minimization"
        simulation.minimizeEnergy(maxIterations=min_max_iterations)
        state = simulation.context.getState(getEnergy=True)
        min_energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        result["minimized_energy"] = min_energy
        _log(f"[minimization] Done. Potential energy = {min_energy:.2f} kJ/mol")

        # Save minimized structure
        min_pdb_path = output_dir / "minimized.pdb"
        from openmm.app import PDBFile
        with open(min_pdb_path, "w") as fh:
            PDBFile.writeFile(
                top.topology,
                simulation.context.getState(getPositions=True).getPositions(),
                fh,
            )
        result["output_files"]["minimized_pdb"] = str(min_pdb_path)
        _log(f"[minimization] Saved minimized structure → {min_pdb_path}")

    # ------------------------------------------------------------------
    # 7. Add position restraints (used for NVT + NPT equilibration)
    # ------------------------------------------------------------------
    _log("[restraints] Adding position restraints on protein heavy atoms ...")
    positions = simulation.context.getState(getPositions=True).getPositions()
    restraint_force = _add_position_restraints(
        system, top.topology, positions, force_constant=1000.0
    )
    # Re-create simulation after modifying the system
    integrator_nvt = LangevinMiddleIntegrator(
        temperature * kelvin,
        1.0 / picosecond,
        timestep * picosecond,
    )
    simulation_nvt = Simulation(
        top.topology,
        system,
        integrator_nvt,
        openmm_platform,
        platform_properties,
    )
    simulation_nvt.context.setPositions(positions)

    # Re-minimize with restraints to convergence — critical for membrane systems
    # where limited initial minimization may leave bad contacts near restrained atoms
    _log("[restraints] Re-minimizing with position restraints (to convergence) ...")
    simulation_nvt.minimizeEnergy(maxIterations=0)  # 0 = run until convergence
    state = simulation_nvt.context.getState(getEnergy=True)
    re_min_energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    _log(f"[restraints] Re-minimization done. Potential energy = {re_min_energy:.2f} kJ/mol")

    simulation_nvt.context.setVelocitiesToTemperature(temperature * kelvin)

    # ------------------------------------------------------------------
    # 8. NVT equilibration
    # ------------------------------------------------------------------
    if nvt_steps > 0:
        _log(f"[NVT] Starting NVT equilibration ({nvt_steps} steps = {nvt_steps * timestep:.0f} ps) ...")
        result["stage"] = "nvt"

        nvt_log_path = output_dir / "nvt.csv"
        nvt_dcd_path = output_dir / "nvt.dcd"
        nvt_chk_path = output_dir / "nvt.chk"

        simulation_nvt.reporters.append(
            StateDataReporter(
                str(nvt_log_path),
                report_interval,
                step=True,
                time=True,
                potentialEnergy=True,
                kineticEnergy=True,
                temperature=True,
                speed=True,
            )
        )
        simulation_nvt.reporters.append(DCDReporter(str(nvt_dcd_path), report_interval))

        simulation_nvt.step(nvt_steps)
        simulation_nvt.saveCheckpoint(str(nvt_chk_path))

        result["output_files"]["nvt_log"] = str(nvt_log_path)
        result["output_files"]["nvt_dcd"] = str(nvt_dcd_path)
        result["output_files"]["nvt_checkpoint"] = str(nvt_chk_path)
        _log(f"[NVT] Done. Checkpoint saved → {nvt_chk_path}")

    # ------------------------------------------------------------------
    # 9. NPT equilibration (add membrane barostat)
    # ------------------------------------------------------------------
    if npt_steps > 0:
        _log(f"[NPT] Starting NPT equilibration ({npt_steps} steps = {npt_steps * timestep:.0f} ps) ...")
        result["stage"] = "npt"

        # Add MonteCarloMembraneBarostat for semiisotropic coupling
        barostat = MonteCarloMembraneBarostat(
            pressure * bar,
            0 * bar * nanometer,          # surface tension = 0 for bilayer
            temperature * kelvin,
            MonteCarloMembraneBarostat.XYIsotropic,
            MonteCarloMembraneBarostat.ZFree,
            25,                           # barostat frequency (steps)
        )
        system.addForce(barostat)

        integrator_npt = LangevinMiddleIntegrator(
            temperature * kelvin,
            1.0 / picosecond,
            timestep * picosecond,
        )
        simulation_npt = Simulation(
            top.topology,
            system,
            integrator_npt,
            openmm_platform,
            platform_properties,
        )
        # Transfer positions + velocities from NVT
        npt_positions = simulation_nvt.context.getState(
            getPositions=True, enforcePeriodicBox=True
        ).getPositions()
        npt_velocities = simulation_nvt.context.getState(getVelocities=True).getVelocities()
        simulation_npt.context.setPositions(npt_positions)
        simulation_npt.context.setVelocities(npt_velocities)

        npt_log_path = output_dir / "npt.csv"
        npt_dcd_path = output_dir / "npt.dcd"
        npt_chk_path = output_dir / "npt.chk"

        simulation_npt.reporters.append(
            StateDataReporter(
                str(npt_log_path),
                report_interval,
                step=True,
                time=True,
                potentialEnergy=True,
                kineticEnergy=True,
                temperature=True,
                volume=True,
                speed=True,
            )
        )
        simulation_npt.reporters.append(DCDReporter(str(npt_dcd_path), report_interval))

        for step_start in range(0, npt_steps, checkpoint_interval):
            steps_this_chunk = min(checkpoint_interval, npt_steps - step_start)
            simulation_npt.step(steps_this_chunk)
            simulation_npt.saveCheckpoint(str(npt_chk_path))
            _log(f"[NPT] {step_start + steps_this_chunk}/{npt_steps} steps completed.")

        result["output_files"]["npt_log"] = str(npt_log_path)
        result["output_files"]["npt_dcd"] = str(npt_dcd_path)
        result["output_files"]["npt_checkpoint"] = str(npt_chk_path)
        _log(f"[NPT] Done. Checkpoint saved → {npt_chk_path}")
    else:
        # No NPT — keep NVT simulation as the source for production
        simulation_npt = simulation_nvt

    # ------------------------------------------------------------------
    # 10. Release position restraints for production
    # ------------------------------------------------------------------
    _log("[production] Releasing position restraints (k → 0) ...")
    # Zero the global force constant so all restrained particles are free
    simulation_npt.context.setParameter("k", 0.0)

    # ------------------------------------------------------------------
    # 11. Production MD
    # ------------------------------------------------------------------
    if production_steps > 0:
        _log(
            f"[production] Starting production MD "
            f"({production_steps} steps = {production_steps * timestep:.0f} ps) ..."
        )
        result["stage"] = "production"

        prod_log_path = output_dir / "production.csv"
        prod_dcd_path = output_dir / "production.dcd"
        prod_chk_path = output_dir / "production.chk"

        simulation_npt.reporters.clear()
        simulation_npt.reporters.append(
            StateDataReporter(
                str(prod_log_path),
                report_interval,
                step=True,
                time=True,
                potentialEnergy=True,
                kineticEnergy=True,
                temperature=True,
                volume=True,
                density=True,
                speed=True,
            )
        )
        simulation_npt.reporters.append(DCDReporter(str(prod_dcd_path), report_interval))

        for step_start in range(0, production_steps, checkpoint_interval):
            steps_this_chunk = min(checkpoint_interval, production_steps - step_start)
            simulation_npt.step(steps_this_chunk)
            simulation_npt.saveCheckpoint(str(prod_chk_path))
            _log(
                f"[production] {step_start + steps_this_chunk}/{production_steps} steps completed."
            )

        result["output_files"]["production_log"] = str(prod_log_path)
        result["output_files"]["production_dcd"] = str(prod_dcd_path)
        result["output_files"]["production_checkpoint"] = str(prod_chk_path)

        # Final state
        final_state = simulation_npt.context.getState(getEnergy=True)
        # Temperature is not directly in State; read last line of CSV
        try:
            with open(prod_log_path) as fh:
                csv_lines = [ln for ln in fh.readlines() if not ln.startswith("#")]
            if len(csv_lines) > 1:
                last = csv_lines[-1].strip().split(",")
                # Column order: step, time, potE, kinE, temp, vol, density, speed
                result["final_temperature"] = float(last[4]) if len(last) > 4 else None
        except Exception:
            pass

        _log(f"[production] Done. Trajectory saved → {prod_dcd_path}")

    result["success"] = True
    result["stage"] = "done"
    result["log"] = "\n".join(log_lines[-50:])

    _log("Simulation complete.")
    return result


# ---------------------------------------------------------------------------
# Public async interface
# ---------------------------------------------------------------------------


async def run_openmm_simulation(
    gro_path: str,
    top_path: str,
    output_dir: str = "./md_output",
    minimize: bool = True,
    min_max_iterations: int = 0,     # 0 = run until convergence (recommended)
    nvt_steps: int = 50000,          # 100 ps at 2 fs
    npt_steps: int = 500000,         # 1 ns at 2 fs
    production_steps: int = 5000000, # 10 ns at 2 fs
    temperature: float = 310.0,      # K
    pressure: float = 1.0,           # bar
    timestep: float = 0.002,         # ps
    nonbonded_cutoff: float = 1.2,   # nm
    report_interval: int = 5000,
    checkpoint_interval: int = 25000,
    platform: str = "CUDA",
    precision: str = "mixed",
    gromacs_include_dir: str | None = None,
) -> dict:
    """Run a multi-stage OpenMM MD simulation for a membrane protein system.

    Stages executed (in order):

    1. **Energy minimization** — steepest-descent until ``min_max_iterations``
       or convergence, whichever comes first.
    2. **NVT equilibration** — constant-volume, constant-temperature run with
       harmonic position restraints on all protein heavy atoms (non-hydrogen
       ATOM records, excluding water/ion residues).
    3. **NPT equilibration** — constant-pressure run using
       ``MonteCarloMembraneBarostat`` (semiisotropic XY, free Z), still with
       position restraints.
    4. **Production MD** — restraints released, full membrane barostat active.

    The simulation runs entirely on a background thread so the calling async
    event loop is never blocked.

    Parameters
    ----------
    gro_path:
        Path to the GROMACS .gro coordinate file.
    top_path:
        Path to the GROMACS .top topology file.
    output_dir:
        Directory where all output files will be written.  Created if absent.
    minimize:
        Whether to perform energy minimization before equilibration.
    min_max_iterations:
        Maximum number of minimization iterations (0 = unlimited).
    nvt_steps:
        Number of NVT equilibration steps (default 50 000 = 100 ps at 2 fs).
    npt_steps:
        Number of NPT equilibration steps (default 500 000 = 1 ns at 2 fs).
    production_steps:
        Number of production MD steps (default 5 000 000 = 10 ns at 2 fs).
    temperature:
        Simulation temperature in Kelvin (default 310 K).
    pressure:
        Target pressure in bar (default 1.0 bar).
    timestep:
        Integration timestep in picoseconds (default 0.002 ps = 2 fs).
    nonbonded_cutoff:
        PME real-space cutoff in nm (default 1.2 nm).
    report_interval:
        Steps between trajectory/energy writes (default 5 000).
    checkpoint_interval:
        Steps between checkpoint saves within a stage (default 25 000).
    platform:
        Preferred OpenMM compute platform (``"CUDA"``, ``"OpenCL"``,
        ``"CPU"``).  Falls back to OpenCL then CPU if the requested platform
        is unavailable.
    precision:
        Floating-point precision for GPU platforms (``"mixed"``, ``"single"``,
        ``"double"``).
    gromacs_include_dir:
        Optional path to a GROMACS ``include/`` directory containing force-field
        .itp files referenced by the topology.  Passed directly to
        ``GromacsTopFile(includeDir=...)``.

    Returns
    -------
    dict
        ``success`` (bool), ``stage`` (last completed stage or the failed
        stage name), ``output_files`` (dict mapping label → absolute path),
        ``minimized_energy`` (kJ/mol or None), ``final_temperature`` (K or
        None), ``log`` (last 50 lines of progress messages).

    Raises
    ------
    ImportError
        If OpenMM is not installed.  The error message contains install
        instructions.
    FileNotFoundError
        If the provided .gro or .top file does not exist.
    RuntimeError
        If no usable OpenMM platform is found.
    """
    if not _OPENMM_AVAILABLE:
        raise ImportError(
            "OpenMM is not installed in the current Python environment.\n"
            "Install it with:\n"
            "  conda install -n AmberTools25 -c conda-forge openmm\n"
            "then restart the server."
        )

    out_dir = Path(output_dir).resolve()

    try:
        return await asyncio.to_thread(
            _run_simulation,
            gro_path=gro_path,
            top_path=top_path,
            output_dir=out_dir,
            minimize=minimize,
            min_max_iterations=min_max_iterations,
            nvt_steps=nvt_steps,
            npt_steps=npt_steps,
            production_steps=production_steps,
            temperature=temperature,
            pressure=pressure,
            timestep=timestep,
            nonbonded_cutoff=nonbonded_cutoff,
            report_interval=report_interval,
            checkpoint_interval=checkpoint_interval,
            platform=platform,
            precision=precision,
            gromacs_include_dir=gromacs_include_dir,
        )
    except Exception as exc:
        logger.exception("OpenMM simulation failed: %s", exc)
        return {
            "success": False,
            "stage": "error",
            "output_files": {},
            "minimized_energy": None,
            "final_temperature": None,
            "log": str(exc),
        }
