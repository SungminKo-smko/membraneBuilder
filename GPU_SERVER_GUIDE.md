# GPU 서버 이어서 작업 가이드

## 현재 상태 (2026-04-16)

### 완료된 작업
1. **PDB 전처리 모듈** (`src/membrane_builder_mcp/pdb_preprocessor.py`)
   - 7CFN에서 chain R (GPCR19) + 리간드 (FX0, CLR, PLM) 추출 완료
   - HOH 제거, CONECT 리맵핑, atom serial 재번호 지원
   - 테스트 5/5 통과

2. **OpenMM MD Runner** (`src/membrane_builder_mcp/openmm_runner.py`)
   - GROMACS .gro/.top 파일 → OpenMM MD simulation 파이프라인
   - 4단계: Energy Minimization → NVT (restraints) → NPT (membrane barostat) → Production
   - CUDA/OpenCL/CPU 자동 fallback
   - 코드 작성 완료, 단위 테스트 통과 (파일 로딩, 시스템 생성 검증)

3. **MCP 서버 업데이트** (`src/membrane_builder_mcp/server.py`)
   - `extract_protein_chains` 도구 추가
   - `run_md_simulation` 도구 추가

4. **기존 GROMACS 파일** (`examples/output_7cfn/gromacs/`)
   - `membrane_7CFN.gro` — 234,348 atoms (GPCR19 + POPC bilayer + TIP3P + ions)
   - `membrane_7CFN.top` — ParmEd standalone topology (AMBER ff19SB + Lipid21 + GAFF2)
   - mdp 파일: em.mdp, nvt.mdp, npt.mdp, md.mdp (GROMACS용 참조)

### 미완료 — GPU 서버에서 해야 할 일

#### 1. 환경 설정
```bash
# conda 환경에 OpenMM + CUDA 설치
conda create -n AmberTools25 python=3.12
conda install -n AmberTools25 -c conda-forge ambertools openmm cudatoolkit

# 또는 기존 환경에 추가
conda install -n AmberTools25 -c conda-forge openmm cudatoolkit

# CUDA 플랫폼 확인
python -c "
from openmm import Platform
for i in range(Platform.getNumPlatforms()):
    p = Platform.getPlatform(i)
    print(f'{p.getName()} (speed={p.getSpeed()})')
"
# 기대 출력: CUDA (speed=...) 가 목록에 있어야 함
```

#### 2. Energy Minimization 테스트 (최우선)
CPU에서 234k 원자 시스템의 minimization이 NaN 문제로 실패함.
- **원인**: CPU에서 200 iteration만으로는 에너지가 수렴하지 않음 (initial E ≈ +10^10 kJ/mol, 50iter 후 7.5M kJ/mol)
- **수정 사항**: `min_max_iterations=0` (convergence까지 실행)으로 변경 완료
- **코드 수정**: restraint 추가 후 re-minimization도 convergence 기반으로 변경

```bash
# GPU에서 빠른 convergence 테스트
python -c "
import sys, time
sys.path.insert(0, 'src')
from openmm.app import *
from openmm import *
from openmm.unit import *

gro = GromacsGroFile('examples/output_7cfn/gromacs/membrane_7CFN.gro')
top = GromacsTopFile('examples/output_7cfn/gromacs/membrane_7CFN.top',
                     periodicBoxVectors=gro.getPeriodicBoxVectors())
system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.2*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(310*kelvin, 1.0/picosecond, 0.002*picosecond)

# CUDA 플랫폼 사용
platform = Platform.getPlatformByName('CUDA')
props = {'Precision': 'mixed'}
sim = Simulation(top.topology, system, integrator, platform, props)
sim.context.setPositions(gro.positions)

state = sim.context.getState(getEnergy=True)
print(f'Initial energy: {state.getPotentialEnergy()}')

t0 = time.time()
sim.minimizeEnergy(maxIterations=0)  # convergence
t1 = time.time()
state = sim.context.getState(getEnergy=True)
print(f'Converged energy: {state.getPotentialEnergy()}')
print(f'Time: {t1-t0:.1f}s')
"
```

#### 3. 전체 MD 파이프라인 테스트
```bash
python -c "
import asyncio, sys, json
sys.path.insert(0, 'src')
from membrane_builder_mcp.openmm_runner import run_openmm_simulation

async def main():
    result = await run_openmm_simulation(
        gro_path='examples/output_7cfn/gromacs/membrane_7CFN.gro',
        top_path='examples/output_7cfn/gromacs/membrane_7CFN.top',
        output_dir='examples/output_7cfn/md_output',
        minimize=True,
        min_max_iterations=0,       # convergence
        nvt_steps=50000,            # 100ps
        npt_steps=500000,           # 1ns
        production_steps=5000000,   # 10ns
        temperature=310.0,
        pressure=1.0,
        timestep=0.002,
        nonbonded_cutoff=1.2,
        report_interval=5000,
        checkpoint_interval=25000,
        platform='CUDA',
        precision='mixed',
    )
    print(json.dumps(result, indent=2))

asyncio.run(main())
"
```

#### 4. 검증 체크리스트
- [ ] CUDA 플랫폼 사용 가능 확인
- [ ] Energy minimization convergence (에너지 음수로 수렴)
- [ ] NVT equilibration — 온도 310K ± 10K 안정화
- [ ] NPT equilibration — 압력 1 bar, 부피 안정화
- [ ] Production MD — trajectory.dcd 생성, energy/temperature CSV 정상
- [ ] 전체 파이프라인 pytest 통과

#### 5. NaN 문제 발생 시 디버깅
```python
# 1. 에너지 확인
state = sim.context.getState(getEnergy=True, getForces=True)
print(f"Energy: {state.getPotentialEnergy()}")

# 2. 최대 힘 확인 (큰 값 = 나쁜 접촉)
import numpy as np
forces = state.getForces(asNumpy=True)
max_force = np.max(np.linalg.norm(forces.value_in_unit(kilojoule_per_mole/nanometer), axis=1))
print(f"Max force: {max_force:.0f} kJ/mol/nm")

# 3. NaN 위치 찾기
positions = sim.context.getState(getPositions=True).getPositions(asNumpy=True)
nan_indices = np.where(np.isnan(positions.value_in_unit(nanometer)))[0]
print(f"NaN atom indices: {nan_indices}")
```

## 프로젝트 아키텍처

```
src/membrane_builder_mcp/
├── server.py              # MCP 서버 (도구 6개)
│   ├── build_membrane         — packmol-memgen 막 빌드
│   ├── list_available_lipids  — 지원 지질 목록
│   ├── analyze_protein        — PDB 구조 분석
│   ├── convert_to_gromacs     — AMBER → GROMACS 변환
│   ├── extract_protein_chains — PDB chain/ligand 추출 (NEW)
│   └── run_md_simulation      — OpenMM MD simulation (NEW)
├── membrane_builder.py    # packmol-memgen CLI wrapper
├── converter.py           # AMBER → GROMACS (ParmEd)
├── pdb_preprocessor.py    # PDB chain extractor (NEW)
└── openmm_runner.py       # OpenMM MD runner (NEW)
```

## 7CFN 전체 워크플로우

```
7CFN.pdb (Cryo-EM, 5 chains)
    │  extract_protein_chains(chains=["R"], keep_ligands=["FX0","CLR","PLM"])
    ▼
7CFN_GPCR19.pdb (chain R only, 2076 ATOM + 110 HETATM)
    │  build_membrane(lipids="POPC", preoriented=True, convert_to_gromacs=True)
    ▼
membrane_7CFN.gro + .top (234,348 atoms, GROMACS format)
    │  run_md_simulation(platform="CUDA")
    ▼
md_output/
    ├── minimized.pdb          — 최소화된 구조
    ├── nvt.dcd + nvt.csv      — NVT equilibration
    ├── npt.dcd + npt.csv      — NPT equilibration
    ├── production.dcd + .csv  — Production MD trajectory
    └── production.chk         — Checkpoint (restart용)
```

## 환경 변수
| 변수 | 기본값 | 설명 |
|------|--------|------|
| `MEMBRANE_CONDA_ENV` | `AmberTools25` | conda 환경 이름 |
| `CONDA_BASE` | `~/miniconda3` | conda 설치 경로 |
