"""Literature-guided peak → VOC remapping for malaria breath GC-MS.

Grounded in:
- Schaber et al., J Infect Dis 2018 (ST000883 / PR000612) — pediatric Malawi;
  α-pinene, 3-carene (terpenes) elevated; thioethers largely absent in this cohort
- Berna et al., J Infect Dis 2015 — adult CHMI thioethers + cyclohexanone
- Pediatric 6-VOC panel reviews — nonanal, isoprene, tridecane, methyl-alkanes

ST000883 library names often encode RT/m/z in parentheses and use ``?-Pinene``
for α-pinene — handled here without changing global atlas physiology priors.
"""

from __future__ import annotations

import re
from typing import Optional

# Exact / prefix aliases after paren-stripping + lowercasing
_MALARIA_LIT_ALIASES: dict[str, str] = {
    # Terpenes (Schaber 2018 mosquito-attractant / infection markers)
    "alpha-pinene": "alpha_pinene",
    "α-pinene": "alpha_pinene",
    "?-pinene": "alpha_pinene",
    "a-pinene": "alpha_pinene",
    "pinene": "alpha_pinene",
    "3-carene": "delta_3_carene",
    "delta-3-carene": "delta_3_carene",
    "δ-3-carene": "delta_3_carene",
    # Limonene IUPAC-ish library string (partial after paren wipe)
    "limonene": "limonene",
    # Berna 2015 / adult CHMI
    "cyclohexanone": "cyclohexanone",
    "allyl methyl sulfide": "allyl_methyl_sulfide",
    "allyl methyl sulphide": "allyl_methyl_sulfide",
    "1-methylthio-propane": "methylthio_propane",
    "propane, 1-(methylthio)-": "methylthio_propane",
    "(z)-1-methylthio-1-propene": "methylthio_propene",
    "(e)-1-methylthio-1-propene": "methylthio_propene",
    "1-propene, 3-(methylthio)-": "allyl_methyl_sulfide",
    # Pediatric alkane / aldehyde panel adjuncts
    "tridecane": "tridecane",
    "undecane": "undecane",
    "undecane, 4-methyl": "methyl_undecane",
    "decane, 2,5,9-trimethyl-": "trimethyl_decane",
    "decane, 2,5,9-trimethyl": "trimethyl_decane",
    "p-xylene": "xylene",
    "m-xylene": "xylene",
    "o-xylene": "xylene",
    "xylene": "xylene",
}

# Substring rules (order matters — more specific first)
_MALARIA_LIT_SUBSTRINGS: list[tuple[str, str]] = [
    ("methylthio-1-propene", "methylthio_propene"),
    ("methylthio-propane", "methylthio_propane"),
    ("methylthio", "allyl_methyl_sulfide"),  # last-resort sulfur hit
    ("allyl methyl", "allyl_methyl_sulfide"),
    ("cyclohexanone", "cyclohexanone"),
    ("3-carene", "delta_3_carene"),
    ("pinene", "alpha_pinene"),
    ("1-methylethenyl", "limonene"),  # limonene IUPAC fragment
    ("tridecane", "tridecane"),
    ("trimethyl-1-nonene", "trimethyl_alkene"),
    ("4-methyl", "methyl_undecane"),  # only with undecane checked below
    ("p-xylene", "xylene"),
    ("xylene", "xylene"),
]


def _strip_lib_name(name: str) -> str:
    # Drop (m/z@rt) library suffixes; keep Greek / ?
    s = re.sub(r"\([^)]*@[^)]*\)", "", str(name))
    s = re.sub(r"\([^)]*\)", "", s)
    s = s.replace("α", "alpha").replace("δ", "delta").replace("β", "beta")
    s = re.sub(r"\s+", " ", s).strip().lower()
    return s


def map_malaria_lit_voc(name: str) -> Optional[str]:
    """Map a GC-MS library metabolite name to a malaria-literature voc_id."""
    n = _strip_lib_name(name)
    if not n:
        return None
    if n in _MALARIA_LIT_ALIASES:
        return _MALARIA_LIT_ALIASES[n]
    # undecane family
    if "undecane" in n:
        if "methyl" in n:
            return "methyl_undecane"
        return "undecane"
    if "decane" in n and "trimethyl" in n:
        return "trimethyl_decane"
    for key, vid in _MALARIA_LIT_SUBSTRINGS:
        if key in n:
            # avoid mapping random "4-methyl" alcohols
            if key == "4-methyl" and "undecane" not in n:
                continue
            if key == "methylthio" and "borane" in n:
                # ST000883 "Borane-methyl sulfide complex" is a library artifact —
                # keep as tentative sulfur proxy only when no better hit
                return "dms"
            return vid
    return None


def map_voc_malaria_expanded(name: str, base_map) -> Optional[str]:
    """Compose atlas `_map_voc` with malaria literature remapping."""
    v = base_map(name)
    if v:
        return v
    return map_malaria_lit_voc(name)


# Literature signature directions for transferable scoring on remapped peaks
# (log2fc-ish: + elevate in malaria). Used when hybrid/stack lack these voc_ids.
MALARIA_LIT_SIGNATURE: dict[str, float] = {
    "alpha_pinene": 1.2,
    "delta_3_carene": 1.1,
    "limonene": 0.6,
    "cyclohexanone": 0.8,
    "allyl_methyl_sulfide": 1.4,
    "methylthio_propane": 1.3,
    "methylthio_propene": 1.3,
    "tridecane": 0.9,
    "methyl_undecane": 0.8,
    "trimethyl_decane": 0.7,
    "undecane": 0.5,
    "nonanal": 0.7,
    "isoprene": 0.5,
    "acetone": 0.3,
    "benzene": 0.4,
    "hexanal": 0.4,
    "pentane": 0.3,
    "dms": 0.5,
    "xylene": 0.3,
}


__all__ = [
    "MALARIA_LIT_SIGNATURE",
    "map_malaria_lit_voc",
    "map_voc_malaria_expanded",
]
