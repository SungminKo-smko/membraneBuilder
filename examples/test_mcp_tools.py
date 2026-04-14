"""MCP 서버 3개 tool의 동등(equivalent) 테스트 스크립트.

MCP 서버를 거치지 않고 Python 모듈을 직접 import하여 테스트합니다.
테스트 대상 tool:
  1. analyze_protein_structure  → pdb_utils.suggest_membrane_params()
  2. build_membrane (기본값)    → membrane_builder.build_membrane()
  3. build_membrane (커스텀 조성) → membrane_builder.build_membrane()
"""

import asyncio
import json
import os
import sys

# 패키지 경로를 sys.path에 추가
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from membrane_builder_mcp.pdb_utils import suggest_membrane_params
from membrane_builder_mcp.membrane_builder import build_membrane

# ──────────────────────────────────────────────────────────────────────────────
# 공통 설정
# ──────────────────────────────────────────────────────────────────────────────

PDB_PATH = os.path.join(os.path.dirname(__file__), "7CFN_clean.pdb")

DIVIDER = "=" * 70


def pretty(obj) -> str:
    """dict/list를 보기 좋게 JSON 형식으로 출력."""
    return json.dumps(obj, indent=2, default=str)


# ──────────────────────────────────────────────────────────────────────────────
# 시나리오 1: analyze_protein_structure 동등 테스트
# ──────────────────────────────────────────────────────────────────────────────

def scenario_1_analyze_protein():
    print(DIVIDER)
    print("시나리오 1: analyze_protein_structure 동등 테스트")
    print(f"  pdb_path = {PDB_PATH}")
    print(DIVIDER)

    result = suggest_membrane_params(PDB_PATH)

    print("[결과]")
    print(pretty(result))
    print()
    print("  membrane_x_size       :", result["membrane_x_size"], "Å")
    print("  membrane_y_size       :", result["membrane_y_size"], "Å")
    print("  n_lipids_per_leaflet  :", result["n_lipids_per_leaflet"])
    print("  cross_section_area    :", result["cross_section_area"], "Å²")
    print("  protein_center        :", result["protein_center"])
    bbox = result["protein_bbox"]
    print("  protein_bbox.size     :", bbox["size"])
    print()


# ──────────────────────────────────────────────────────────────────────────────
# 시나리오 2: build_membrane 기본값 테스트
# ──────────────────────────────────────────────────────────────────────────────

async def scenario_2_build_default():
    print(DIVIDER)
    print("시나리오 2: build_membrane 기본값 테스트")
    print("  lipid_composition = {'POPC': 1.0}  (기본값)")
    output_dir = os.path.join(os.path.dirname(__file__), "output_default")
    print(f"  output_dir = {output_dir}")
    print(DIVIDER)

    result = await build_membrane(
        protein_pdb_path=PDB_PATH,
        output_dir=output_dir,
    )

    print("[결과]")
    print(pretty(result))
    print()

    # 생성된 .inp 파일 내용 확인
    inp_path = result.get("input_path")
    if inp_path and os.path.exists(inp_path):
        print(f"[생성된 .inp 파일: {inp_path}]")
        with open(inp_path, "r") as f:
            inp_content = f.read()
        print(inp_content)
    else:
        print(f"[경고] .inp 파일을 찾을 수 없습니다: {inp_path}")
    print()


# ──────────────────────────────────────────────────────────────────────────────
# 시나리오 3: build_membrane 커스텀 조성 테스트
# ──────────────────────────────────────────────────────────────────────────────

async def scenario_3_build_custom():
    print(DIVIDER)
    print("시나리오 3: build_membrane 커스텀 조성 테스트")
    custom_composition = {"POPC": 0.5, "POPE": 0.3, "cholesterol": 0.2}
    print(f"  lipid_composition = {custom_composition}")
    output_dir = os.path.join(os.path.dirname(__file__), "output_custom")
    print(f"  output_dir = {output_dir}")
    print(DIVIDER)

    result = await build_membrane(
        protein_pdb_path=PDB_PATH,
        output_dir=output_dir,
        lipid_composition=custom_composition,
    )

    print("[결과]")
    print(pretty(result))
    print()

    # 리피드 분포 요약
    print("[리피드 분포 요약 (upper leaflet 기준)]")
    for lipid, count in result.get("n_lipids_upper", {}).items():
        fraction = custom_composition.get(lipid, 0)
        print(f"  {lipid:12s}: {count:4d} 개  (fraction={fraction:.1f})")
    print()


# ──────────────────────────────────────────────────────────────────────────────
# 메인
# ──────────────────────────────────────────────────────────────────────────────

async def main():
    print()
    print("MCP Tool 직접 호출 테스트")
    print(f"PDB 파일: {PDB_PATH}")
    print()

    # 시나리오 1: 동기 함수
    scenario_1_analyze_protein()

    # 시나리오 2, 3: async 함수
    await scenario_2_build_default()
    await scenario_3_build_custom()

    print(DIVIDER)
    print("모든 테스트 완료.")
    print(DIVIDER)


if __name__ == "__main__":
    asyncio.run(main())
