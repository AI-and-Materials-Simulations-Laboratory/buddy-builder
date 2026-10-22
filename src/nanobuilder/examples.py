"""Live examples as configuration dicts: the five required ones (spec section 17,
A-E), five carbon-nanomaterial ones (F-J) and two Buddy Builder ones (K: slit
channels, L: Fibonacci-ordered dodecanethiol shell in a cubic box).

``examples/*.json`` in the tool directory hold the same configurations; a
test keeps the two in sync.
"""

from __future__ import annotations

import io
from pathlib import Path

EXAMPLES = {
    "A_bare_au_sphere": {
        "label": "A - Bare Au/FCC sphere, 3 nm",
        "config": {
            "material": {"source": "preset", "preset_id": "Au_fcc"},
            "geometry": {"shape": "sphere", "diameter_nm": 3.0, "vacuum_A": 30.0},
            "coverage": {"mode": "count", "target": 0, "region": "entire_surface"},
        },
    },
    "B_au_sphere_12_ligands": {
        "label": "B - Au sphere + 12 SCCO (Au-thiolate)",
        "config": {
            "material": {"source": "preset", "preset_id": "Au_fcc"},
            "geometry": {"shape": "sphere", "diameter_nm": 3.0, "vacuum_A": 30.0},
            "ligand": {"smiles": "SCCO", "anchor_atom_index": 0, "binding_template": "Au_thiolate",
                       "distance_A": 2.4, "conformers": 1, "seed": 42},
            "coverage": {"mode": "count", "target": 12, "region": "entire_surface",
                         "distribution": "distributed_random", "seed": 42},
        },
    },
    "C_au111_slab_density": {
        "label": "C - Au(111) slab, 0.5 ligands/nm^2 on top",
        "config": {
            "material": {"source": "preset", "preset_id": "Au_fcc"},
            "geometry": {"shape": "slab", "miller": [1, 1, 1], "slab_length_nm": 5.0,
                         "slab_width_nm": 5.0, "thickness_mode": "layers", "layers": 4,
                         "vacuum_A": 30.0},
            "ligand": {"smiles": "SCCO", "anchor_atom_index": 0, "binding_template": "Au_thiolate",
                       "distance_A": 2.4, "seed": 42},
            "coverage": {"mode": "density", "density_nm2": 0.5, "region": "top",
                         "distribution": "evenly_spaced", "seed": 42},
        },
    },
    "D_fe_bcc_sphere": {
        "label": "D - Fe/BCC sphere, 2.5 nm (bare)",
        "config": {
            "material": {"source": "preset", "preset_id": "Fe_bcc"},
            "geometry": {"shape": "sphere", "diameter_nm": 2.5, "vacuum_A": 30.0},
            "ligand": {"binding_template": "geometric"},
            "coverage": {"mode": "count", "target": 0, "region": "entire_surface"},
        },
    },
    "E_uploaded_mg_cif_slab": {
        "label": "E - Uploaded Mg HCP CIF, (001) slab + ethanol",
        "upload": "mg_hcp_fully_occupied.cif",
        "config": {
            "material": {"source": "uploaded_file", "file_format": "cif",
                         "filename": "mg_hcp_fully_occupied.cif"},
            "geometry": {"shape": "slab", "miller": [0, 0, 1], "slab_length_nm": 3.0,
                         "slab_width_nm": 3.0, "thickness_mode": "layers", "layers": 3,
                         "vacuum_A": 30.0},
            "ligand": {"smiles": "OCC", "anchor_atom_index": 0, "binding_template": "geometric",
                       "distance_A": 2.5, "surface_species": ["Mg"], "seed": 42},
            "coverage": {"mode": "count", "target": 6, "region": "top",
                         "distribution": "evenly_spaced", "seed": 42},
        },
    },
    "F_graphene_oxide_flake": {
        "label": "F - Graphene oxide flake, 30 % O (fresh 2:1, AMBER GO)",
        "config": {
            "material": {"source": "carbon"},
            "geometry": {"vacuum_A": 20.0, "max_atoms": 100000},
            "carbon": {"kind": "graphene", "length_nm": 5.0, "width_nm": 5.0, "periodic": False,
                       "oxygen_percent": 30.0, "oxidation_profile": "fresh",
                       "hydroxyl_epoxy_ratio": "2:1", "seed": 7},
            "coverage": {"mode": "count", "target": 0, "region": "both"},
        },
    },
    "G_aged_go_amide": {
        "label": "G - Aged GO with holes + ethylenediamine amides",
        "config": {
            "material": {"source": "carbon"},
            "geometry": {"vacuum_A": 20.0, "max_atoms": 100000},
            "carbon": {"kind": "graphene", "length_nm": 6.0, "width_nm": 6.0,
                       "vacancy_percent": 8.0, "stone_wales_percent": 2.0, "oxygen_percent": 20.0,
                       "oxidation_profile": "aged_1200K", "oxidation_pattern": "edges_first",
                       "seed": 11},
            "ligand": {"smiles": "NCCN", "anchor_atom_index": 0, "binding_template": "C_carboxyl",
                       "conformer_mode": "extended", "seed": 42, "library_id": "ethylenediamine"},
            "coverage": {"mode": "count", "target": 6, "region": "edges", "seed": 42},
        },
    },
    "H_cnt_aryl": {
        "label": "H - Periodic SWCNT (10,10) + 4-carboxyphenyl",
        "config": {
            "material": {"source": "carbon"},
            "geometry": {"vacuum_A": 20.0, "max_atoms": 100000},
            "carbon": {"kind": "nanotube", "cnt_n": 10, "cnt_m": 10, "cnt_length_nm": 5.0,
                       "periodic": True, "seed": 7},
            "ligand": {"smiles": "c1ccc(cc1)C(=O)O", "anchor_atom_index": 0,
                       "binding_template": "C_basal", "seed": 42, "library_id": "carboxyphenyl"},
            "coverage": {"mode": "count", "target": 10, "region": "side",
                         "distribution": "evenly_spaced", "seed": 42},
        },
    },
    "I_c60_amine": {
        "label": "I - C60 fullerene + 6 butylamines",
        "config": {
            "material": {"source": "carbon"},
            "geometry": {"vacuum_A": 20.0, "max_atoms": 100000},
            "carbon": {"kind": "carbon_sphere", "sphere_type": "fullerene", "fullerene_order": 1},
            "ligand": {"smiles": "NCCCC", "anchor_atom_index": 0, "binding_template": "C_basal",
                       "conformer_mode": "extended", "seed": 42},
            "coverage": {"mode": "count", "target": 6, "region": "entire_surface",
                         "distribution": "evenly_spaced", "seed": 42},
        },
    },
    "J_periodic_go_epoxide": {
        "label": "J - Periodic GO sheet, epoxides opened by ethanolamine",
        "config": {
            "material": {"source": "carbon"},
            "geometry": {"vacuum_A": 20.0, "max_atoms": 100000},
            "carbon": {"kind": "graphene", "length_nm": 5.0, "width_nm": 5.0, "periodic": True,
                       "oxygen_percent": 15.0, "oxidation_profile": "fresh",
                       "hydroxyl_epoxy_ratio": "1:1", "seed": 3},
            "ligand": {"smiles": "NCCO", "anchor_atom_index": 0, "binding_template": "C_epoxide",
                       "conformer_mode": "extended", "seed": 42, "library_id": "ethanolamine"},
            "coverage": {"mode": "count", "target": 6, "region": "top", "seed": 42},
        },
    },
    "K_au_slit_channels": {
        "label": "K - Au(111) slit channels x2, 16 SCCO on the walls, 3.5 x 3 x 8 nm box",
        "config": {
            "material": {"source": "preset", "preset_id": "Au_fcc"},
            "geometry": {"shape": "channel", "channel_type": "slit", "channel_count": 2,
                         "channel_width_nm": 1.5, "layers": 4, "miller": [1, 1, 1],
                         "box_shape": "rectangular", "box_lengths_nm": [3.5, 3.0, 8.0]},
            "ligand": {"smiles": "SCCO", "anchor_atom_index": 0, "binding_template": "Au_thiolate",
                       "distance_A": 2.4, "seed": 42},
            "coverage": {"mode": "count", "target": 16, "region": "inner",
                         "distribution": "fibonacci", "seed": 42},
        },
    },
    "L_au_sphere_fibonacci_c12": {
        "label": "L - Au 4 nm sphere + 60 dodecanethiols, Fibonacci order, 12.5 nm cube",
        "config": {
            "material": {"source": "preset", "preset_id": "Au_fcc"},
            "geometry": {"shape": "sphere", "diameter_nm": 4.0, "box_shape": "cube",
                         "box_lengths_nm": [12.5, 12.5, 12.5]},
            "ligand": {"smiles": "SCCCCCCCCCCCC", "anchor_atom_index": 0,
                       "binding_template": "Au_thiolate", "distance_A": 2.4, "seed": 42,
                       "conformer_mode": "extended", "library_id": "dodecanethiol"},
            "coverage": {"mode": "count", "target": 60, "region": "entire_surface",
                         "distribution": "fibonacci", "seed": 42},
        },
    },
}

# JSON copies in examples/ (kept in sync by tests/test_integration.py)
EXAMPLE_FILES = {
    "A_bare_au_sphere": "au_sphere_bare.json",
    "B_au_sphere_12_ligands": "au_sphere_12ligands.json",
    "C_au111_slab_density": "au111_slab_density.json",
    "D_fe_bcc_sphere": "fe_bcc_sphere.json",
    "E_uploaded_mg_cif_slab": "uploaded_crystal_slab.json",
    "F_graphene_oxide_flake": "graphene_oxide_flake.json",
    "G_aged_go_amide": "aged_go_holes_amide.json",
    "H_cnt_aryl": "cnt_aryl.json",
    "I_c60_amine": "c60_amine.json",
    "J_periodic_go_epoxide": "periodic_go_epoxide_opening.json",
    "K_au_slit_channels": "au_slit_channels.json",
    "L_au_sphere_fibonacci_c12": "au_sphere_fibonacci_c12.json",
}


def fixture_dir():
    """``examples/fixtures`` of the tool directory, if this copy can find it."""
    here = Path(__file__).resolve().parent
    for parent in [here, *here.parents][:5]:
        candidate = parent / "examples" / "fixtures"
        if candidate.is_dir():
            return candidate
    return None


def mg_hcp_cif_bytes() -> bytes:
    """The valid noncubic fixture; regenerated with ASE if the file is not found."""
    folder = fixture_dir()
    if folder and (folder / "mg_hcp_fully_occupied.cif").is_file():
        return (folder / "mg_hcp_fully_occupied.cif").read_bytes()
    from ase.build import bulk
    from ase.io import write

    buf = io.BytesIO()
    write(buf, bulk("Mg", "hcp", a=3.21, c=5.21), format="cif")
    return buf.getvalue()


def example_upload(key: str):
    name = EXAMPLES[key].get("upload")
    if name == "mg_hcp_fully_occupied.cif":
        return name, mg_hcp_cif_bytes()
    return None
