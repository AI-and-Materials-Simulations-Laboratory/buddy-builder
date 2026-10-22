"""Attachment templates and the Au-thiolate availability rule (spec stage C1).

This module is the single definition of which chemistry is offered for which
substrate. Nothing here assigns force-field parameters.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .schema import BuilderError, CrystalModel


@dataclass(frozen=True)
class AttachmentTemplate:
    key: str
    label: str
    default_distance_A: float
    covalent_surface_bond: bool
    bond_type: str = ""
    required_anchor_element: str = ""
    surface_species: tuple = field(default_factory=tuple)
    description: str = ""
    substrate: str = "any"  # "metal" (Au preset), "carbon", or "any"
    site_kind: str = "on_top"  # carbon: basal, edge, carboxyl, epoxide
    anchor_elements: tuple = field(default_factory=tuple)  # allowed anchor elements (carbon)
    removes_anchor_h: bool = False  # the anchor loses one H to form the bond


# Typical bond lengths (A) from a surface carbon to the ligand anchor, by site kind:
# sp3 basal carbon, sp2 edge (aryl) carbon, acyl carbon of a carboxyl (amide/ester),
# and the sp3 carbon of an opened epoxide.
COVALENT_LENGTHS = {
    "basal": {"C": 1.53, "N": 1.47, "O": 1.43, "S": 1.82, "P": 1.87, "Si": 1.89, "B": 1.58},
    "edge": {"C": 1.50, "N": 1.40, "O": 1.36, "S": 1.77, "P": 1.82, "Si": 1.86, "B": 1.56},
    "carboxyl": {"N": 1.34, "O": 1.34, "S": 1.77},
    "epoxide": {"N": 1.47, "O": 1.43, "S": 1.82},
}
CARBON_TEMPLATES = ("C_basal", "C_edge", "C_carboxyl", "C_epoxide")


TEMPLATES = {
    "Au_thiolate": AttachmentTemplate(
        key="Au_thiolate",
        label="Au-thiolate template",
        default_distance_A=2.4,
        covalent_surface_bond=True,
        bond_type="Au-S",
        required_anchor_element="S",
        surface_species=("Au",),
        description=("Thiol sulfur loses its hydrogen and is placed on top of a surface Au atom "
                     "with an explicit Au-S topology bond. Starting geometry only; the Au-S "
                     "interaction still needs a force field."),
    ),
    "geometric": AttachmentTemplate(
        key="geometric",
        label="Physisorption (no bond)",
        default_distance_A=2.5,
        covalent_surface_bond=False,
        description=("The intact ligand sits at the chosen anchor-to-surface distance, held by "
                     "non-bonded interactions only: no surface bond is added and no atoms are "
                     "removed."),
    ),
    "chemisorption": AttachmentTemplate(
        key="chemisorption",
        label="Chemisorption (anchor-surface bond)",
        default_distance_A=2.3,
        covalent_surface_bond=True,
        description=("The anchor bonds the surface atom below it: an explicit anchor-surface "
                     "bond is written to the topology, at the sum of covalent radii by default. "
                     "A thiol, alcohol or acid anchor (S-H, O-H) gives up its H (thiolate, "
                     "alkoxide, carboxylate); N and P anchors bind through their lone pair and "
                     "keep their H. Starting geometry only; the bond needs force-field "
                     "parameters."),
    ),
    "C_basal": AttachmentTemplate(
        key="C_basal",
        label="Covalent: basal-plane addition (C-X, sp2 to sp3)",
        default_distance_A=1.53,
        covalent_surface_bond=True,
        surface_species=("C",),
        substrate="carbon",
        site_kind="basal",
        anchor_elements=("C", "N", "O", "S", "P", "Si", "B"),
        removes_anchor_h=True,
        description=("Radical or nitrene-type addition to an unmodified basal carbon (e.g. aryl "
                     "groups from diazonium salts, alkyl or amine radicals). The anchor loses one H; "
                     "the surface carbon becomes sp3 (moved 0.40 A outward) and bonds the anchor "
                     "along the local normal."),
    ),
    "C_edge": AttachmentTemplate(
        key="C_edge",
        label="Covalent: edge substitution (replaces an edge H)",
        default_distance_A=1.50,
        covalent_surface_bond=True,
        surface_species=("C",),
        substrate="carbon",
        site_kind="edge",
        anchor_elements=("C", "N", "O", "S", "P", "Si", "B"),
        removes_anchor_h=True,
        description=("The ligand replaces a hydrogen on a sheet, hole or nanotube-rim edge and "
                     "bonds the sp2 edge carbon in the plane. Needs a finite flake/tube or holes "
                     "with hydrogen-terminated edges; choose the Edges (or Tips) region."),
    ),
    "C_carboxyl": AttachmentTemplate(
        key="C_carboxyl",
        label="Covalent: amide/ester coupling to GO carboxyl",
        default_distance_A=1.34,
        covalent_surface_bond=True,
        surface_species=("C",),
        substrate="carbon",
        site_kind="carboxyl",
        anchor_elements=("N", "O", "S"),
        removes_anchor_h=True,
        description=("EDC/NHS-type coupling: a carboxylic acid of the oxidized substrate loses its "
                     "OH and the ligand's amine (amide), alcohol (ester) or thiol (thioester) "
                     "loses one H. Needs carboxyl groups (aged or custom oxidation profile)."),
    ),
    "C_epoxide": AttachmentTemplate(
        key="C_epoxide",
        label="Covalent: epoxide ring-opening (amine, alcohol, thiol)",
        default_distance_A=1.47,
        covalent_surface_bond=True,
        surface_species=("C",),
        substrate="carbon",
        site_kind="epoxide",
        anchor_elements=("N", "O", "S"),
        removes_anchor_h=True,
        description=("Nucleophilic opening of a GO epoxide: the ligand bonds one epoxide carbon "
                     "and the oxygen stays on the other carbon as a hydroxyl (an H is added). On "
                     "flat sheets the ligand adds from the opposite face (anti, SN2-like); on "
                     "curved walls both groups point outward, tilted apart."),
    ),
}


def template_distance(template_key: str, anchor_element: str = "", surface_element: str = "") -> float:
    """Default anchor-to-surface distance for a template, anchor and surface element."""
    t = TEMPLATES[template_key]
    if template_key == "chemisorption":
        from ase.data import atomic_numbers, covalent_radii
        if anchor_element in atomic_numbers and surface_element in atomic_numbers:
            return round(float(covalent_radii[atomic_numbers[anchor_element]]
                               + covalent_radii[atomic_numbers[surface_element]]), 2)
        return t.default_distance_A
    if t.substrate == "carbon":
        table = COVALENT_LENGTHS[t.site_kind]
        if anchor_element in table:
            return table[anchor_element]
        from ase.data import atomic_numbers, covalent_radii
        if anchor_element in atomic_numbers:
            return round(float(covalent_radii[6] + covalent_radii[atomic_numbers[anchor_element]]), 2)
    return t.default_distance_A


def thiolate_allowed(crystal: CrystalModel) -> bool:
    """Au-thiolate availability rule (spec section 8, C1).

    Offered only for the built-in Au/FCC preset (an edited lattice constant
    keeps it, with a warning elsewhere). Every other preset and every upload,
    including an Au CIF, gets geometric adsorption only.
    """
    return crystal.source_type == "preset" and crystal.preset_id == "Au_fcc"


def template_keys(crystal: CrystalModel) -> list:
    """Attachment templates offered for this substrate, in GUI order."""
    if crystal.source_type == "carbon":
        return list(CARBON_TEMPLATES) + ["geometric"]
    if thiolate_allowed(crystal):
        return ["Au_thiolate", "chemisorption", "geometric"]
    return ["geometric", "chemisorption"]


def available_templates(crystal: CrystalModel) -> list:
    return [TEMPLATES[k] for k in template_keys(crystal)]


def thiolate_unavailable_reason(crystal: CrystalModel) -> str:
    if thiolate_allowed(crystal):
        return ""
    if crystal.source_type == "carbon":
        return ""
    if crystal.source_type == "materials_project":
        return ("Au-thiolate is only reviewed for the built-in Au/FCC preset; Materials Project "
                "structures use physisorption or generic chemisorption")
    if crystal.source_type == "uploaded_file":
        return ("Au-thiolate is disabled for uploaded crystals (no reviewed structure or "
                "termination); use physisorption or generic chemisorption")
    return (f"Au-thiolate is only reviewed for the Au/FCC preset, not {crystal.structure_label}; "
            "use physisorption or generic chemisorption")


def check_template(template_key: str, crystal: CrystalModel) -> AttachmentTemplate:
    if template_key not in TEMPLATES:
        raise BuilderError(f"Unknown attachment template {template_key!r}")
    if template_key == "Au_thiolate" and not thiolate_allowed(crystal):
        if crystal.source_type == "carbon":
            raise BuilderError("Au-thiolate needs the Au/FCC preset; for carbon substrates choose a "
                               "covalent carbon template or physisorption")
        raise BuilderError(thiolate_unavailable_reason(crystal))
    if template_key == "chemisorption" and crystal.source_type == "carbon":
        raise BuilderError("For carbon substrates choose one of the covalent carbon templates "
                           "(they form the C-X bond) or physisorption")
    if TEMPLATES[template_key].substrate == "carbon" and crystal.source_type != "carbon":
        raise BuilderError(f"{TEMPLATES[template_key].label} needs a carbon nanomaterial substrate "
                           "(graphene, graphene oxide, nanotube or fullerene)")
    return TEMPLATES[template_key]


def resolve_surface_species(template: AttachmentTemplate, requested, crystal: CrystalModel) -> list:
    """Which substrate elements may host an anchor."""
    candidates = list(crystal.anchor_species_candidates)
    if template.surface_species:
        allowed = [s for s in template.surface_species if s in candidates]
        if requested and set(requested) - set(allowed):
            raise BuilderError(f"{template.label} only binds to {', '.join(allowed)}")
        return allowed
    if requested is None:  # not specified: every element of the crystal may anchor
        return candidates
    if not requested:  # explicitly empty, e.g. every GUI box unticked
        raise BuilderError("No surface anchor species selected: tick at least one of "
                           f"{', '.join(candidates)} under Surface anchor species")
    unknown = set(requested) - set(candidates)
    if unknown:
        raise BuilderError(f"Surface species {', '.join(sorted(unknown))} not in this crystal "
                           f"({', '.join(candidates)})")
    return sorted(requested)
