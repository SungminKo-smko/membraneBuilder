"""AMBER to GROMACS format converter using ParmEd (AmberTools25)."""

import asyncio
import logging
import os
import subprocess
import tempfile
from pathlib import Path

logger = logging.getLogger(__name__)

CONDA_ENV = os.environ.get("MEMBRANE_CONDA_ENV", "AmberTools25")
CONDA_BASE = Path(os.environ.get("CONDA_BASE", str(Path.home() / "miniconda3")))
ENV_BIN = CONDA_BASE / "envs" / CONDA_ENV / "bin"


def _get_env() -> dict[str, str]:
    """Return env dict with conda env bin/ prepended to PATH."""
    env = os.environ.copy()
    env["PATH"] = f"{ENV_BIN}:{env.get('PATH', '')}"
    env["AMBERHOME"] = str(CONDA_BASE / "envs" / CONDA_ENV)
    return env


def _run_parmed(script_path: str) -> subprocess.CompletedProcess:
    """ParmEd를 직접 실행하고 결과를 반환."""
    parmed_exe = str(ENV_BIN / "parmed")
    cmd = [parmed_exe, "-i", script_path]
    logger.debug("Running: %s", " ".join(cmd))
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=300,
        env=_get_env(),
    )
    return result


async def convert_amber_to_gromacs(
    amber_prmtop: str,
    amber_inpcrd: str,
    output_dir: str | None = None,
    output_prefix: str = "system",
) -> dict:
    """AMBER 토폴로지+좌표를 GROMACS .top + .gro로 변환한다.

    ParmEd를 conda run으로 호출하여 변환 수행.

    Args:
        amber_prmtop: AMBER topology file (.prmtop 또는 .top)
        amber_inpcrd: AMBER coordinate file (.inpcrd 또는 .crd)
        output_dir: 출력 디렉토리 (None이면 prmtop과 같은 디렉토리)
        output_prefix: 출력 파일 접두사

    Returns:
        dict with:
            - success: bool
            - gromacs_top: str (경로)
            - gromacs_gro: str (경로)
            - log: str
    """
    prmtop_path = Path(amber_prmtop).resolve()
    inpcrd_path = Path(amber_inpcrd).resolve()

    # 1. output_dir 결정
    if output_dir is None:
        out_dir = prmtop_path.parent
    else:
        out_dir = Path(output_dir).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)

    gro_out = out_dir / f"{output_prefix}.gro"
    top_out = out_dir / f"{output_prefix}.top"

    # 2. ParmEd 스크립트 생성
    script_lines = [
        f"parm {prmtop_path}",
        f"loadRestrt {inpcrd_path}",
        f"outparm {top_out} {gro_out}",
        "quit",
    ]
    script_content = "\n".join(script_lines) + "\n"
    logger.debug("ParmEd script:\n%s", script_content)

    # 3. 스크립트를 임시 파일로 저장
    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".parmed.in",
        delete=False,
        dir=out_dir,
    ) as tmp:
        tmp.write(script_content)
        tmp_path = tmp.name

    # 4. conda run으로 ParmEd 실행 (asyncio.to_thread로 래핑)
    try:
        result = await asyncio.to_thread(_run_parmed, tmp_path)
    except FileNotFoundError as exc:
        logger.error("conda/parmed not found: %s", exc)
        return {
            "success": False,
            "gromacs_top": str(top_out),
            "gromacs_gro": str(gro_out),
            "log": str(exc),
        }
    except subprocess.TimeoutExpired as exc:
        logger.error("ParmEd timed out: %s", exc)
        return {
            "success": False,
            "gromacs_top": str(top_out),
            "gromacs_gro": str(gro_out),
            "log": "ParmEd process timed out after 300 seconds.",
        }
    finally:
        # 임시 스크립트 파일 정리
        try:
            Path(tmp_path).unlink(missing_ok=True)
        except Exception:
            pass

    log_text = (result.stdout or "") + (result.stderr or "")

    # 5. 결과 확인 및 반환
    if result.returncode != 0:
        logger.error("ParmEd failed (rc=%d):\n%s", result.returncode, log_text)
        return {
            "success": False,
            "gromacs_top": str(top_out),
            "gromacs_gro": str(gro_out),
            "log": log_text,
        }

    if not gro_out.exists() or not top_out.exists():
        missing = [str(p) for p in (gro_out, top_out) if not p.exists()]
        msg = f"ParmEd succeeded but output files missing: {missing}\n{log_text}"
        logger.error(msg)
        return {
            "success": False,
            "gromacs_top": str(top_out),
            "gromacs_gro": str(gro_out),
            "log": msg,
        }

    logger.info(
        "Conversion successful: %s, %s",
        top_out,
        gro_out,
    )
    return {
        "success": True,
        "gromacs_top": str(top_out),
        "gromacs_gro": str(gro_out),
        "log": log_text,
    }


async def validate_gromacs_files(gro_path: str, top_path: str) -> dict:
    """생성된 GROMACS 파일의 유효성을 간단히 검증.

    - .gro: 파일 존재, 원자 수 확인 (2번째 줄)
    - .top: 파일 존재, [ system ] 섹션 존재 확인

    Returns:
        dict with:
            - valid: bool
            - atom_count: int or None  (.gro에서 읽은 값)
            - has_system_section: bool
            - errors: list[str]
    """
    errors: list[str] = []
    atom_count: int | None = None
    has_system_section = False

    gro = Path(gro_path)
    top = Path(top_path)

    # --- .gro 검증 ---
    if not gro.exists():
        errors.append(f".gro file not found: {gro_path}")
    else:
        try:
            with gro.open() as fh:
                lines = fh.readlines()
            if len(lines) < 2:
                errors.append(f".gro file too short ({len(lines)} lines): {gro_path}")
            else:
                raw = lines[1].strip()
                try:
                    atom_count = int(raw)
                    if atom_count <= 0:
                        errors.append(
                            f".gro reports non-positive atom count: {atom_count}"
                        )
                except ValueError:
                    errors.append(
                        f".gro second line is not an integer (got {raw!r}): {gro_path}"
                    )
        except OSError as exc:
            errors.append(f"Cannot read .gro file: {exc}")

    # --- .top 검증 ---
    if not top.exists():
        errors.append(f".top file not found: {top_path}")
    else:
        try:
            with top.open() as fh:
                for line in fh:
                    stripped = line.strip()
                    # Match "[ system ]" with flexible spacing
                    if stripped.startswith("[") and "system" in stripped.lower():
                        has_system_section = True
                        break
            if not has_system_section:
                errors.append(
                    f"[ system ] section not found in .top file: {top_path}"
                )
        except OSError as exc:
            errors.append(f"Cannot read .top file: {exc}")

    valid = len(errors) == 0
    if valid:
        logger.info(
            "GROMACS files validated: %d atoms, [ system ] section present.",
            atom_count,
        )
    else:
        logger.warning("GROMACS validation failed: %s", errors)

    return {
        "valid": valid,
        "atom_count": atom_count,
        "has_system_section": has_system_section,
        "errors": errors,
    }
