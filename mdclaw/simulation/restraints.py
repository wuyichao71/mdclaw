"""Shared positional-restraint atom selection for simulation nodes."""

import json
import math
import re
from pathlib import Path
from typing import Any, Optional

from mdclaw.chemistry_constants import WATER_NAMES, is_standard_bare_ion_resname


RESTRAINT_SELECTIONS = ("solute_heavy", "CA", "backbone", "heavy")
_RESTRAINT_TYPES = ("distance", "angle", "dihedral")
_RESTRAINT_GROUP_COUNT = {"distance": 2, "angle": 3, "dihedral": 4}
_RESTRAINT_TARGET_FIELD = {
    "distance": "target_distance_nm",
    "angle": "target_angle_deg",
    "dihedral": "target_angle_deg",
}
_RESTRAINT_FORCE_CONSTANT_FIELD = {
    "distance": "force_constant_kj_mol_nm2",
    "angle": "force_constant_kj_mol_rad2",
    "dihedral": "force_constant_kj_mol_rad2",
}
# A dihedral difference must be folded into (-pi, pi] before it is squared:
# a plain harmonic on the raw difference is discontinuous at the +/-180 deg
# wrap and produces a large spurious force there.
_TWO_PI_LITERAL = "6.283185307179586"
_RESTRAINT_ENERGY = {
    "distance": "0.5*k*(distance(g1,g2)-r0)^2",
    "angle": "0.5*k*(angle(g1,g2,g3)-r0)^2",
    "dihedral": (
        "0.5*k*dphi^2;"
        f"dphi=dtheta-{_TWO_PI_LITERAL}*floor(0.5+dtheta/{_TWO_PI_LITERAL});"
        "dtheta=dihedral(g1,g2,g3,g4)-r0"
    ),
}

_DISTANCE_RESTRAINT_NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
_BACKBONE_NAMES = {"N", "CA", "C", "O"}
_SOLUTE_COMPONENT_TYPES = {"protein", "nucleic", "glycan", "ligand", "ion"}
# Used only when prep provenance is unavailable. The canonical path selects
# prep-derived chains by index and therefore does not classify residue names.
_COMMON_LIPID_RESNAMES = {
    "PA", "PC", "PE", "PGR", "PG", "PS", "PSER", "OL",
    "POPC", "POPE", "POPG", "POPS", "DOPC", "DOPE", "DOPG", "DOPS",
    "DPPC", "DPPE", "DPPG", "DMPC", "DSPC", "DLPC", "CHL", "CHL1",
}


class DistanceRestraintError(RuntimeError):
    """Structured validation error for declarative distance restraints."""

    def __init__(self, *, code: str, message: str):
        super().__init__(message)
        self.code = code


def _restraint_fields(restraint_type: str) -> set[str]:
    """Required field names for one restraint type (``type`` stays optional)."""
    return {
        "name",
        _RESTRAINT_FORCE_CONSTANT_FIELD[restraint_type],
        _RESTRAINT_TARGET_FIELD[restraint_type],
    } | {
        f"selection_group{index}"
        for index in range(1, _RESTRAINT_GROUP_COUNT[restraint_type] + 1)
    }


def normalize_distance_restraints(
    distance_restraints: Optional[list[dict]],
) -> Optional[list[dict]]:
    """Validate and normalize the topology-independent restraint schema."""
    if distance_restraints is None:
        return None
    if not isinstance(distance_restraints, list) or not distance_restraints:
        raise DistanceRestraintError(
            code="distance_restraints_invalid",
            message="distance_restraints must be a non-empty list of objects.",
        )

    normalized: list[dict] = []
    names: set[str] = set()
    for index, item in enumerate(distance_restraints):
        if not isinstance(item, dict):
            raise DistanceRestraintError(
                code="distance_restraints_invalid",
                message=f"distance_restraints[{index}] must be an object.",
            )
        restraint_type = item.get("type", "distance")
        if restraint_type not in _RESTRAINT_TYPES:
            raise DistanceRestraintError(
                code="restraint_type_invalid",
                message=(
                    f"distance_restraints[{index}].type must be one of "
                    + ", ".join(_RESTRAINT_TYPES)
                    + f"; got {restraint_type!r}."
                ),
            )
        required = _restraint_fields(restraint_type)
        # Only ``type`` is optional. Nothing derived may enter the normalized
        # payload: node conditions are cross-checked against it verbatim, so an
        # extra key the caller did not declare fails the node.
        allowed = required | {"type"}
        missing = required - set(item)
        unknown = set(item) - allowed
        if missing or unknown:
            details = []
            if missing:
                details.append("missing " + ", ".join(sorted(missing)))
            if unknown:
                details.append("unknown " + ", ".join(sorted(unknown)))
            raise DistanceRestraintError(
                code="distance_restraints_invalid",
                message=f"distance_restraints[{index}] has " + "; ".join(details) + ".",
            )

        name = item["name"]
        if not isinstance(name, str) or not _DISTANCE_RESTRAINT_NAME_RE.fullmatch(name):
            raise DistanceRestraintError(
                code="distance_restraints_invalid",
                message=(
                    f"distance_restraints[{index}].name must match "
                    "[A-Za-z][A-Za-z0-9_]*."
                ),
            )
        if name in names:
            raise DistanceRestraintError(
                code="distance_restraints_invalid",
                message=f"distance restraint name {name!r} is duplicated.",
            )
        names.add(name)

        selections = {}
        for group_index in range(1, _RESTRAINT_GROUP_COUNT[restraint_type] + 1):
            key = f"selection_group{group_index}"
            value = item[key]
            if not isinstance(value, str) or not value.strip():
                raise DistanceRestraintError(
                    code="distance_restraints_invalid",
                    message=f"distance_restraints[{index}].{key} must be a non-empty string.",
                )
            selections[key] = value.strip()

        force_constant_field = _RESTRAINT_FORCE_CONSTANT_FIELD[restraint_type]
        target_field = _RESTRAINT_TARGET_FIELD[restraint_type]
        force_constant = item[force_constant_field]
        target = item[target_field]
        if (
            isinstance(force_constant, bool)
            or not isinstance(force_constant, (int, float))
            or not math.isfinite(float(force_constant))
            or float(force_constant) <= 0.0
        ):
            raise DistanceRestraintError(
                code="distance_restraints_invalid",
                message=(
                    f"distance_restraints[{index}].{force_constant_field} "
                    "must be finite and greater than 0."
                ),
            )
        if (
            isinstance(target, bool)
            or not isinstance(target, (int, float))
            or not math.isfinite(float(target))
        ):
            raise DistanceRestraintError(
                code="distance_restraints_invalid",
                message=(
                    f"distance_restraints[{index}].{target_field} must be "
                    "a finite number."
                ),
            )
        target = float(target)
        if restraint_type == "distance" and target < 0.0:
            raise DistanceRestraintError(
                code="distance_restraints_invalid",
                message=(
                    f"distance_restraints[{index}].target_distance_nm must be "
                    "finite and greater than or equal to 0."
                ),
            )
        if restraint_type == "angle" and not 0.0 <= target <= 180.0:
            raise DistanceRestraintError(
                code="restraint_target_out_of_range",
                message=(
                    f"distance_restraints[{index}].target_angle_deg must be "
                    "between 0 and 180 for an angle restraint."
                ),
            )
        if restraint_type == "dihedral" and not -180.0 <= target <= 180.0:
            raise DistanceRestraintError(
                code="restraint_target_out_of_range",
                message=(
                    f"distance_restraints[{index}].target_angle_deg must be "
                    "between -180 and 180 for a dihedral restraint."
                ),
            )

        entry = {"name": name}
        if restraint_type != "distance":
            # Kept out of the payload for plain distance restraints so their
            # normalized form and signature stay byte-identical to schema v1.
            entry["type"] = restraint_type
        entry.update(selections)
        entry[force_constant_field] = float(force_constant)
        entry[target_field] = target
        normalized.append(entry)
    return normalized


def distance_restraint_signature(
    distance_restraints: Optional[list[dict]],
) -> Optional[dict]:
    """Return the reproducibility signature for declarative distance bias."""
    normalized = normalize_distance_restraints(distance_restraints)
    if normalized is None:
        return None
    types = sorted({item.get("type", "distance") for item in normalized})
    if types == ["distance"]:
        # Schema v1 signature preserved byte-for-byte so previously recorded
        # node conditions keep cross-checking.
        return {
            "kind": "openmm_centroid_distance_restraints",
            "mass_weighting": "physical_element",
            "restraints": normalized,
        }
    return {
        "kind": "openmm_centroid_restraints",
        "mass_weighting": "physical_element",
        "types": types,
        "angle_cv_unit": "radian",
        "restraints": normalized,
    }


def load_distance_restraints(
    *,
    system,
    topology,
    distance_restraints: list[dict],
    is_periodic: bool,
) -> dict:
    """Build one native OpenMM harmonic COM-distance bias force."""
    import mdtraj as md
    import numpy as np
    from openmm import CustomCentroidBondForce
    from openmm.unit import dalton

    normalized = normalize_distance_restraints(distance_restraints)
    mdtraj_topology = md.Topology.from_openmm(topology)
    if mdtraj_topology.n_atoms != system.getNumParticles():
        raise DistanceRestraintError(
            code="distance_restraint_topology_mismatch",
            message=(
                f"topology.pdb has {mdtraj_topology.n_atoms} atoms but system.xml "
                f"has {system.getNumParticles()} particles."
            ),
        )

    topology_atoms = list(topology.atoms())
    groups: list[tuple[list[list[int]], list[list[float]]]] = []
    for item in normalized:
        restraint_type = item.get("type", "distance")
        selected: list[list[int]] = []
        selected_weights: list[list[float]] = []
        for group_index in range(1, _RESTRAINT_GROUP_COUNT[restraint_type] + 1):
            key = f"selection_group{group_index}"
            try:
                indices = [int(value) for value in mdtraj_topology.select(item[key])]
            except Exception as exc:
                raise DistanceRestraintError(
                    code="distance_restraint_selection_invalid",
                    message=(
                        f"distance restraint {item['name']!r} has an invalid "
                        f"{key}: {exc}"
                    ),
                ) from exc
            if not indices:
                raise DistanceRestraintError(
                    code="restraint_selection_empty",
                    message=(
                        f"distance restraint {item['name']!r} {key}={item[key]!r} "
                        "matched zero atoms."
                    ),
                )
            selected_residues = {
                topology_atoms[atom_index].residue.index:
                topology_atoms[atom_index].residue
                for atom_index in indices
            }
            water_count = sum(
                residue.name.strip().upper() in WATER_NAMES
                for residue in selected_residues.values()
            )
            ion_count = sum(
                len(list(residue.atoms())) == 1
                and is_standard_bare_ion_resname(residue.name.strip())
                for residue in selected_residues.values()
            )
            if water_count or ion_count:
                raise DistanceRestraintError(
                    code="distance_restraints_invalid",
                    message=(
                        f"distance restraint {item['name']!r} {key} matched "
                        f"{water_count} water residue(s) and {ion_count} "
                        "bare-ion residue(s). On a solvated topology, use "
                        "resid rather than resSeq so wrapped/reused PDB residue "
                        "numbers cannot select solvent."
                    ),
                )
            weights = [
                (
                    topology_atoms[atom_index].element.mass.value_in_unit(dalton)
                    if topology_atoms[atom_index].element is not None
                    else 0.0
                )
                for atom_index in indices
            ]
            if sum(weights) <= 0.0:
                raise DistanceRestraintError(
                    code="distance_restraints_invalid",
                    message=(
                        f"distance restraint {item['name']!r} {key} has zero "
                        "total particle mass."
                    ),
                )
            selected.append(indices)
            selected_weights.append(weights)
        for first in range(len(selected)):
            for second in range(first + 1, len(selected)):
                overlap = sorted(set(selected[first]) & set(selected[second]))
                if overlap:
                    raise DistanceRestraintError(
                        code="distance_restraint_groups_overlap",
                        message=(
                            f"distance restraint {item['name']!r} groups "
                            f"{first + 1} and {second + 1} overlap at "
                            f"{len(overlap)} atoms; use disjoint groups."
                        ),
                    )
        groups.append((selected, selected_weights))

    # One CustomCentroidBondForce carries a single energy expression and a
    # single groups-per-bond count, so each restraint type needs its own force.
    forces = []
    for restraint_type in _RESTRAINT_TYPES:
        indices = [
            position
            for position, item in enumerate(normalized)
            if item.get("type", "distance") == restraint_type
        ]
        if not indices:
            continue
        force = CustomCentroidBondForce(
            _RESTRAINT_GROUP_COUNT[restraint_type],
            _RESTRAINT_ENERGY[restraint_type],
        )
        force.addPerBondParameter("k")
        force.addPerBondParameter("r0")
        force.setUsesPeriodicBoundaryConditions(bool(is_periodic))
        for position in indices:
            item = normalized[position]
            selected, selected_weights = groups[position]
            # Use physical elemental masses rather than the System particle
            # masses: HMR changes both hydrogen and bonded-heavy-atom masses,
            # but the scientific COM coordinate must not change under HMR.
            group_ids = [
                force.addGroup(atom_indices, weights)
                for atom_indices, weights in zip(selected, selected_weights)
            ]
            force.addBond(
                group_ids,
                [
                    item[_RESTRAINT_FORCE_CONSTANT_FIELD[restraint_type]],
                    item["target_distance_nm"]
                    if restraint_type == "distance"
                    # The force works in radians; the schema takes degrees.
                    else math.radians(item["target_angle_deg"]),
                ],
            )
        forces.append(force)

    def _minimum_image(displacement, box_np):
        if is_periodic and box_np is not None:
            fractional = displacement @ np.linalg.inv(box_np)
            displacement = displacement - np.rint(fractional) @ box_np
        return displacement

    def _evaluator(positions_np, box_np):
        values = {}
        for item, group in zip(normalized, groups):
            restraint_type = item.get("type", "distance")
            selected, selected_weights = group
            centers = [
                np.average(positions_np[atom_indices], axis=0, weights=weights)
                for atom_indices, weights in zip(selected, selected_weights)
            ]
            if restraint_type == "distance":
                displacement = _minimum_image(centers[1] - centers[0], box_np)
                values[item["name"]] = float(np.linalg.norm(displacement))
                continue
            # Angle and dihedral CVs are reported in radians so that the
            # logged value pairs directly with force_constant_kj_mol_rad2.
            if restraint_type == "angle":
                first = _minimum_image(centers[0] - centers[1], box_np)
                second = _minimum_image(centers[2] - centers[1], box_np)
                cosine = float(
                    np.dot(first, second)
                    / (np.linalg.norm(first) * np.linalg.norm(second))
                )
                values[item["name"]] = float(
                    np.arccos(min(1.0, max(-1.0, cosine)))
                )
                continue
            # Sign convention must match OpenMM's dihedral(): the first bond
            # vector points from atom 2 back to atom 1, otherwise the reported
            # CV comes out negated relative to the bias actually applied.
            b0 = _minimum_image(centers[0] - centers[1], box_np)
            b1 = _minimum_image(centers[2] - centers[1], box_np)
            b2 = _minimum_image(centers[3] - centers[2], box_np)
            axis = b1 / np.linalg.norm(b1)
            v = b0 - np.dot(b0, axis) * axis
            w = b2 - np.dot(b2, axis) * axis
            values[item["name"]] = float(
                np.arctan2(np.dot(np.cross(axis, v), w), np.dot(v, w))
            )
        return values

    signature = distance_restraint_signature(normalized)
    return {
        "forces": forces,
        "evaluator": _evaluator,
        "cv_names": [item["name"] for item in normalized],
        "kind": signature["kind"],
        "signature": signature,
        "restraints": normalized,
    }


def _is_heavy_atom(atom) -> bool:
    return atom.element is not None and atom.element.symbol != "H"


def _is_legacy_solute_atom(atom) -> bool:
    residue = atom.residue
    resname = residue.name.strip()
    if resname.upper() in WATER_NAMES:
        return False
    residue_atoms = list(residue.atoms())
    return not (
        len(residue_atoms) == 1 and is_standard_bare_ion_resname(resname)
    )


def _component_label(component: dict[str, Any]) -> str:
    component_type = (
        component.get("source_chain_type")
        or component.get("prepared_fragment_role")
        or "unknown"
    )
    if component_type == "nucleic":
        return str(component.get("source_nucleic_subtype") or "nucleic").lower()
    if component_type == "ion":
        return "structural_ion"
    return str(component_type)


def _load_component_map(path: Optional[str]) -> tuple[list[dict], list[str]]:
    """Solute components from prep, in prep's own atom order."""
    if not path:
        return [], []
    try:
        payload = json.loads(Path(path).read_text())
        components = [
            component
            for component in payload.get("components", [])
            if (
                component.get("source_chain_type")
                or component.get("prepared_fragment_role")
            ) in _SOLUTE_COMPONENT_TYPES
        ]
        components.sort(key=lambda c: int(c.get("atom_index_start") or 0))
        return components, []
    except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
        return [], [f"Could not read prep chain_identity_map: {exc}"]


def select_restraint_atoms(
    topology,
    selection: str,
    *,
    chain_identity_map_file: Optional[str] = None,
) -> dict[str, Any]:
    """Return atom indices and provenance for a restraint selection."""
    if selection not in RESTRAINT_SELECTIONS:
        raise ValueError(f"Unknown restraint selection: {selection}")

    if selection == "solute_heavy":
        components, warnings = _load_component_map(chain_identity_map_file)
        if components:
            # Address components by their prep atom-index range, not by chain
            # index. Topology generation does not preserve prep's chain
            # decomposition: when Pablo identifies every residue it emits each
            # ACE/NME cap as a chain of its own, and when it falls back to
            # PDBFile it does not. Chain index N in prep is then a different
            # molecule in the built topology -- and the failure is silent,
            # because the wrong chains still yield plausible-looking counts.
            # Solute atoms keep prep's order and lead the topology (solvent and
            # its virtual sites are appended), so the ranges do carry over.
            atoms = list(topology.atoms())
            indices: list[int] = []
            counts: dict[str, int] = {}
            for component in components:
                start = component.get("atom_index_start")
                end = component.get("atom_index_end_exclusive")
                if start is None or end is None:
                    warnings.append(
                        "prep chain_identity_map component "
                        f"{component.get('component_id')} carries no atom range"
                    )
                    continue
                start, end = int(start), int(end)
                if end > len(atoms):
                    warnings.append(
                        "prep chain_identity_map component "
                        f"{component.get('component_id')} ends at atom {end}, "
                        f"past the {len(atoms)}-atom topology"
                    )
                    continue
                label = _component_label(component)
                for atom in atoms[start:end]:
                    if atom.residue.name.strip().upper() in WATER_NAMES:
                        warnings.append(
                            "prep chain_identity_map component "
                            f"{component.get('component_id')} covers solvent at "
                            f"atom {atom.index}; solute atom order did not carry "
                            "over to the built topology"
                        )
                        break
                    if _is_heavy_atom(atom):
                        indices.append(atom.index)
                        counts[label] = counts.get(label, 0) + 1
            return {
                "success": True,
                "atom_indices": indices,
                "counts_by_component": counts,
                "selection_source": "prep_chain_identity_map",
                "warnings": warnings,
                "errors": [],
            }

        warnings.append(
            "prep chain_identity_map is unavailable; structural and solvent "
            "ions cannot be distinguished, so ions are excluded"
        )
        indices = []
        for atom in topology.atoms():
            resname = atom.residue.name.strip().upper()
            if resname in WATER_NAMES or resname in _COMMON_LIPID_RESNAMES:
                continue
            if not _is_legacy_solute_atom(atom) or not _is_heavy_atom(atom):
                continue
            indices.append(atom.index)
        return {
            "success": True,
            "atom_indices": indices,
            "counts_by_component": {"unclassified_solute": len(indices)},
            "selection_source": "topology_fallback",
            "warnings": warnings,
            "errors": [],
        }

    indices = []
    for atom in topology.atoms():
        if not _is_legacy_solute_atom(atom):
            continue
        if selection == "heavy":
            if not _is_heavy_atom(atom):
                continue
        elif selection == "CA":
            if atom.name != "CA":
                continue
        elif atom.name not in _BACKBONE_NAMES:
            continue
        indices.append(atom.index)
    return {
        "success": True,
        "atom_indices": indices,
        "counts_by_component": {"legacy_solute": len(indices)},
        "selection_source": "legacy_selection",
        "warnings": [],
        "errors": [],
    }


def select_lipid_headgroup_anchors(topology) -> dict[str, Any]:
    """Phosphorus atoms of the lipid headgroups, for a flat-bilayer restraint.

    Minimisation and the first thermalisation exist to close what assembly left
    open: the gap where lipids were carved away from the solute, the seams
    between stacked water slabs, the thin few angstroms at each end of the
    cell. A bilayer with nothing holding it can answer that by bending or
    thinning into those spaces instead, and the picture that would show it is
    the one whose defects took a day to find.

    Restraining the headgroup phosphorus in z alone is what CHARMM-GUI's
    membrane protocol does, and it is the restraint that matches the intent:
    the bilayer keeps its thickness and stays flat, while lipids remain free to
    move in the membrane plane and pack back around the solute — which is the
    relaxation being asked for. A full positional restraint on lipid heavy
    atoms would stop that too.

    Sterols carry no phosphorus and are not anchors; they follow the
    phospholipids they sit between.
    """
    from mdclaw.solvation.constants import lipid21_template_contract

    contract = lipid21_template_contract()
    head_names = {name.upper() for name in contract.head_names}
    whole_names = {name.upper() for name in contract.full_names}
    indices: list[int] = []
    for atom in topology.atoms():
        if atom.element is None or atom.element.symbol != "P":
            continue
        resname = atom.residue.name.strip().upper()
        if resname in head_names or resname in whole_names:
            indices.append(atom.index)
    return {
        "atom_indices": indices,
        "count": len(indices),
        "selection": "lipid_headgroup_phosphorus",
    }
