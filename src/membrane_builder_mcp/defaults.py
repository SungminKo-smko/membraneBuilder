from dataclasses import dataclass
from typing import Optional

DEFAULT_LIPID_COMPOSITION: dict[str, float] = {"POPC": 1.0}
DEFAULT_BILAYER_THICKNESS: float = 40.0  # Å
DEFAULT_WATER_PADDING: float = 15.0  # Å
DEFAULT_LIPID_AREA: float = 65.0  # Å² per lipid
DEFAULT_WATER_DENSITY: float = 0.0334  # molecules/ų
DEFAULT_TOLERANCE: float = 2.0  # Å
DEFAULT_XY_PADDING: float = 15.0  # Å per side beyond protein bbox


@dataclass(frozen=True)
class LipidInfo:
    description: str
    headgroup: Optional[str]
    n_atoms: int
    area: float  # Å²


LIPID_DB: dict[str, LipidInfo] = {
    "POPC": LipidInfo(
        description="1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine",
        headgroup="PC",
        n_atoms=134,
        area=68.3,
    ),
    "POPE": LipidInfo(
        description="1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine",
        headgroup="PE",
        n_atoms=125,
        area=58.0,
    ),
    "DPPC": LipidInfo(
        description="1,2-dipalmitoyl-sn-glycero-3-phosphocholine",
        headgroup="PC",
        n_atoms=130,
        area=64.0,
    ),
    "cholesterol": LipidInfo(
        description="cholesterol",
        headgroup=None,
        n_atoms=74,
        area=40.0,
    ),
}
