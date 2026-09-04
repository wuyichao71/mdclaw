# Production MD: Harmonic Distance / Angle / Dihedral Restraints

Use `--distance-restraints` for harmonic restraints on a distance, an angle, or
a dihedral between mdtraj selections (each group is a mass-weighted centroid).
This route uses native OpenMM forces; do not write a `--custom-force-script`
for a potential this page can express - measured on a 10 k-atom system the
custom-force route is ~4.6x slower for the same harmonic distance bias, and it
cannot run at all on GPUs newer than the container's PyTorch build supports.

`type` selects the restraint and decides which other fields are required. It is
optional and defaults to `distance`, so schema-v1 payloads stay valid.

| `type` | Groups | Force constant | Target |
|---|---|---|---|
| `distance` (default) | `selection_group1`, `selection_group2` | `force_constant_kj_mol_nm2` | `target_distance_nm` (>= 0) |
| `angle` | `selection_group1` .. `selection_group3` | `force_constant_kj_mol_rad2` | `target_angle_deg` (0 .. 180) |
| `dihedral` | `selection_group1` .. `selection_group4` | `force_constant_kj_mol_rad2` | `target_angle_deg` (-180 .. 180) |

Every restraint also needs `name` - a unique identifier that becomes the CV
column name; use letters, digits, and underscores. Pass no other keys: node
conditions are cross-checked against the normalized payload verbatim, so an
undeclared key fails the node.

Groups within one restraint must be pairwise disjoint. Angle and dihedral
groups are ordered: `g1-g2-g3` spans the angle at `g2`, and `g1-g2-g3-g4`
follows the usual dihedral convention about the `g2-g3` bond.

**Units.** Angle targets are given in degrees but `force_constant_kj_mol_rad2`
is per radian squared, and the CV column is logged in **radians** so that
`0.5 * k * (cv - target)^2` reproduces the bias directly. Distances stay in nm
throughout.

**Dihedral periodicity is handled for you.** The angle difference is folded
into (-pi, pi] before it is squared, so the bias stays continuous across the
+/-180 deg wrap and remains a plain harmonic that WHAM/MBAR can reweight with
the same `0.5*k*dx^2` expression used for distances.

Every selection must be non-empty. A one-atom selection per group gives the
plain atom-atom / atom-triple / atom-quadruple coordinate. Multi-atom groups
use physical elemental masses, so the CV is unchanged by HMR. Periodic boundary
handling follows the topology automatically.

**Do not use `resSeq` selections on a solvated topology.** PDB residue numbers
wrap at 9999 and solvated chains reuse them, so a `resSeq` range can silently
select water with the same numbers. Use the topology-wide, 0-based `resid`
index together with `protein`. See `scripts/cv_selection.py` when author residue
numbers must be mapped onto a built topology. Distance-restraint selections
that contain water or bare ions are rejected.

## Create and run one window

Declare the complete restraint list in the prod node's conditions and pass the
same JSON list to `run_production`:

```bash
RESTRAINTS='[{"name":"tm3_tm6","selection_group1":"protein and resid 100 to 120 and not element H","selection_group2":"protein and resid 240 to 260 and not element H","force_constant_kj_mol_nm2":1000.0,"target_distance_nm":1.2}]'

mdclaw create_node --job-dir <job_dir> --node-type prod \
  --parent-node-ids <eq_node_id> --label "tm3_tm6_r0_1.2" \
  --conditions "{\"simulation_time_ns\":100,\"distance_restraints\":$RESTRAINTS}"

mdclaw --job-dir <job_dir> --node-id <prod_node_id> run_production \
  --simulation-time-ns 100 --distance-restraints "$RESTRAINTS"
```

A dihedral restraint looks like this - note that `--conditions` must declare
exactly what the tool will receive:

```bash
RESTRAINTS='[{"name":"phi_ala2","type":"dihedral",
  "selection_group1":"index 4","selection_group2":"index 6",
  "selection_group3":"index 8","selection_group4":"index 14",
  "force_constant_kj_mol_rad2":145.0,"target_angle_deg":-75.0}]'
```

Create one sibling prod node per target distance when running umbrella windows.
Completed windows are immutable; add a later window as another child of the
same completed eq node. Use `submit_array_job` when submitting many independent
window nodes.

### Seeding (required when running umbrella sampling)

**When you are running umbrella sampling, seed every window with a staircase
slow pull first; do not branch all windows from one shared `eq` state.**
Branching every window off the same `eq` node applies the full
restraint instantly, so the CV reaches its target within ~1 ps while the slow
degrees of freedom (backbone dihedrals, secondary structure) cannot follow;
every window then carries the same starting conformation and the PMF is
systematically distorted. Instead chain short biased runs with
`--continue-from`, advancing `target_distance_nm` one window spacing per step
(e.g. 0.04 nm / 500 ps), and seed each window from the pull step at its own
`r0`. Each step's end state is already a valid restart state for that window,
so no artifact is moved outside the DAG.

## Outputs and continuation

The prod node records `metadata.distance_restraints`,
`metadata.distance_restraint_signature`, `artifacts.collective_variables`, and
`artifacts.collective_variables_meta`. The CSV contains the total bias energy
and one exact OpenMM distance column per restraint name.

Use `--continue-from <biased_prod_id>` to extend a window. Omitted
`--distance-restraints` inherits the parent's declaration. Biased continuation
requires the parent's XML `state` artifact; binary checkpoint restart is not
supported. Do not combine `--distance-restraints` with
`--custom-force-script`.

## Codes

| Code | Fix |
|---|---|
| `distance_restraints_invalid` | Supply a non-empty list with all five fields and valid finite values. |
| `distance_restraint_selection_invalid` | Fix the mdtraj selection syntax. |
| `restraint_selection_empty` | Use selections that each match at least one atom. |
| `distance_restraint_groups_overlap` | Make the two groups disjoint. |
| `distance_restraint_topology_mismatch` | Use the matching topology/System artifact pair. |
| `production_bias_conflict` | Choose the distance route or custom script, not both. |
| `production_bias_checkpoint_unsupported` | Restart from the portable XML state. |
