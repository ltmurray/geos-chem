#!/usr/bin/env python3
"""Generate shared tagged-methane name tables for region-by-sector tracers."""

from pathlib import Path


BASE_TAGS = (
    "CH4_OIL", "CH4_GAS", "CH4_COL", "CH4_LIV", "CH4_LDF", "CH4_WST",
    "CH4_RIC", "CH4_OTA", "CH4_BBN", "CH4_WTL", "CH4_SEE", "CH4_LAK",
    "CH4_TER", "CH4_RES", "CH4_USA", "CH4_CAN", "CH4_MEX", "CH4_CAC",
    "CH4_BRA", "CH4_NSA", "CH4_SCN", "CH4_WEU", "CH4_EEU", "CH4_RUS",
    "CH4_MDE", "CH4_NAF", "CH4_WAF", "CH4_CAF", "CH4_EAF", "CH4_SIO",
    "CH4_SAC", "CH4_CAS", "CH4_IND", "CH4_SAS", "CH4_CHN", "CH4_EAS",
    "CH4_SEA", "CH4_OCE", "CH4_ROW",
)
REGIONS = (
    "USA", "CAN", "MEX", "CAC", "BRA", "NSA", "SCN", "WEU", "EEU",
    "RUS", "MDE", "NAF", "WAF", "CAF", "EAF", "SIO", "SAC", "CAS",
    "IND", "SAS", "CHN", "EAS", "SEA", "OCE", "ROW",
)
SECTORS = (
    "OIL", "GAS", "COL", "LIV", "LDF", "WST", "RIC", "OTA", "BBN",
    "WTL", "SEE", "LAK", "TER", "RES",
)
TAGS = BASE_TAGS + tuple(
    f"CH4_{region}_{sector}" for region in REGIONS for sector in SECTORS
)


def declaration(public: bool) -> str:
    visibility = ", PUBLIC" if public else ""
    lines = [
        f"  INTEGER{visibility}, PARAMETER :: N_TAGGED_CH4 = {len(TAGS)}",
        "",
        f"  CHARACTER(LEN=11){visibility}, PARAMETER :: "
        "TAGGED_CH4_NAMES(N_TAGGED_CH4) = (/ CHARACTER(LEN=11) :: &",
    ]
    for start in range(0, len(TAGS), 3):
        chunk = TAGS[start:start + 3]
        suffix = ", &" if start + 3 < len(TAGS) else " /)"
        lines.append("       " + ", ".join(f"'{name}'" for name in chunk) + suffix)
    return "\n".join(lines) + "\n"


def replace_names(path: Path, public: bool) -> None:
    text = path.read_text()
    if public:
        start = text.index("  INTEGER, PUBLIC, PARAMETER :: N_TAGGED_CH4")
        end = text.index("\n\n  !=========", start)
        replacement = declaration(public).rstrip()
    else:
        start = text.index("  INTEGER, PARAMETER :: N_TAGGED_CH4")
        end = text.index("\n\n  INTEGER, PARAMETER :: ind_CH4_tagged_base", start)
        replacement = (
            declaration(public)
            + "  INTEGER :: id_CH4_tagged(N_TAGGED_CH4)"
        )
    path.write_text(text[:start] + replacement + text[end:])


def main() -> None:
    root = Path(__file__).resolve().parents[3]
    replace_names(root / "Headers/state_diag_mod.F90", public=True)
    replace_names(root / "KPP/oxidants/oxidants_Funcs.F90", public=False)


if __name__ == "__main__":
    main()
