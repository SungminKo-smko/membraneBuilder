#!/usr/bin/env python3
"""
PDB 클린업 스크립트
7CFN_raw.pdb -> 7CFN_clean.pdb

- ATOM 레코드: 모두 유지
- HETATM: HOH 및 일반 결정화 첨가물/이온 제거, 실제 리간드 유지
- TER, END 레코드 유지
- CONECT: 리간드 관련만 유지

사용법:
  python clean_pdb.py [입력파일] [출력파일]
"""

import sys
from collections import defaultdict

_DEFAULT_INPUT  = "/Users/kosungmin/workspace/membraneBuilder/examples/7CFN_raw.pdb"
_DEFAULT_OUTPUT = "/Users/kosungmin/workspace/membraneBuilder/examples/7CFN_clean.pdb"

INPUT_PDB  = sys.argv[1] if len(sys.argv) > 1 else _DEFAULT_INPUT
OUTPUT_PDB = sys.argv[2] if len(sys.argv) > 2 else _DEFAULT_OUTPUT

# 제거할 HETATM residue 목록 (물 + 결정화 첨가물 + 일반 이온)
REMOVE_RESIDUES = {
    # 물
    "HOH", "WAT", "DOD",
    # 결정화 첨가물
    "SO4", "PO4", "GOL", "EDO", "PEG", "ACT", "MPD", "DMS", "DMSO",
    "IMD", "TRS", "MES", "EPE", "HEPES", "PGE", "PG4", "1PE",
    "FMT", "ACY", "IPA", "EOH", "EGL",
    # 일반 이온
    "CL",  "NA",  "MG",  "CA",  "ZN",  "K",   "MN",  "FE",
    "CU",  "NI",  "CO",  "CD",  "HG",  "PT",  "SR",  "BA",
    "CS",  "RB",  "LI",  "IOD", "BR",  "FLU",
    # 기타 흔한 첨가물
    # SEP, TPO, PTR은 단백질 폴리펩타이드 체인의 일부(인산화 변형 잔기)이므로 기본 유지
}

def parse_pdb(path):
    with open(path, "r") as f:
        return f.readlines()

def get_residue_name(line):
    """PDB 레코드에서 residue 이름 추출 (열 18-21, 0-indexed 17-20, 최대 4글자)"""
    return line[17:21].strip()

def get_serial(line):
    """ATOM/HETATM serial number (열 1-5, 0-indexed 6-10)"""
    try:
        return int(line[6:11].strip())
    except ValueError:
        return None

def main():
    lines = parse_pdb(INPUT_PDB)

    # 원본 통계
    orig_atom   = sum(1 for l in lines if l.startswith("ATOM"))
    orig_hetatm = sum(1 for l in lines if l.startswith("HETATM"))

    # 유지할 HETATM 잔기 분석
    hetatm_residues = defaultdict(int)
    for line in lines:
        if line.startswith("HETATM"):
            resname = get_residue_name(line)
            hetatm_residues[resname] += 1

    kept_ligands   = {r for r in hetatm_residues if r not in REMOVE_RESIDUES}
    removed_residues = {r for r in hetatm_residues if r in REMOVE_RESIDUES}

    # 유지할 HETATM serial 번호 수집 (CONECT 필터링용)
    kept_serials = set()
    for line in lines:
        if line.startswith("HETATM"):
            resname = get_residue_name(line)
            if resname not in REMOVE_RESIDUES:
                serial = get_serial(line)
                if serial is not None:
                    kept_serials.add(serial)

    # ATOM serial도 kept_serials에 포함 (CONECT 유지 시 참조)
    for line in lines:
        if line.startswith("ATOM"):
            serial = get_serial(line)
            if serial is not None:
                kept_serials.add(serial)

    # 클린업 수행
    output_lines = []
    clean_atom   = 0
    clean_hetatm = 0

    for line in lines:
        record = line[:6].strip()

        if record == "ATOM":
            output_lines.append(line)
            clean_atom += 1

        elif record == "HETATM":
            resname = get_residue_name(line)
            if resname not in REMOVE_RESIDUES:
                output_lines.append(line)
                clean_hetatm += 1

        elif record in ("TER", "END"):
            output_lines.append(line)

        elif record == "CONECT":
            # CONECT 레코드: PDB 고정폭 형식으로 파싱 (열 7-11, 12-16, 17-21, 22-26, 27-31)
            # line.split()은 숫자가 붙어있을 때 파싱 실패하므로 고정폭 슬라이싱 사용
            serials_in_line = []
            for i in range(6, min(31, len(line)), 5):
                tok = line[i:i+5].strip()
                if tok:
                    try:
                        serials_in_line.append(int(tok))
                    except ValueError:
                        pass
            # 모든 serial이 kept_serials에 있으면 유지
            if serials_in_line and all(s in kept_serials for s in serials_in_line):
                output_lines.append(line)

        # REMARK, HEADER 등 기타 레코드는 건너뜀 (필요 시 아래 주석 해제)
        # else:
        #     output_lines.append(line)

    # END가 없으면 추가
    if not any(l.strip() == "END" for l in output_lines):
        output_lines.append("END\n")

    with open(OUTPUT_PDB, "w") as f:
        f.writelines(output_lines)

    # 결과 보고
    print("=" * 60)
    print(f"PDB 클린업 완료: {INPUT_PDB} → {OUTPUT_PDB}")
    print("=" * 60)
    print(f"\n[원본]")
    print(f"  ATOM   레코드 수: {orig_atom:>6}")
    print(f"  HETATM 레코드 수: {orig_hetatm:>6}")
    print(f"  합계            : {orig_atom + orig_hetatm:>6}")

    print(f"\n[클린업 후]")
    print(f"  ATOM   레코드 수: {clean_atom:>6}")
    print(f"  HETATM 레코드 수: {clean_hetatm:>6}")
    print(f"  합계            : {clean_atom + clean_hetatm:>6}")

    print(f"\n[유지된 리간드]")
    for resname in sorted(kept_ligands):
        print(f"  {resname:>4s}  ({hetatm_residues[resname]} atoms)")

    print(f"\n[제거된 항목]")
    for resname in sorted(removed_residues):
        print(f"  {resname:>4s}  ({hetatm_residues[resname]} atoms)")

    removed_total = orig_hetatm - clean_hetatm
    print(f"\n  제거된 HETATM 총계: {removed_total} atoms")
    print(f"\n출력 파일: {OUTPUT_PDB}")
    print("=" * 60)

if __name__ == "__main__":
    main()
