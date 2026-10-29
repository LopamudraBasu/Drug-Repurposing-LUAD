# Molecular Docking Methods and Results

## Objective

Molecular docking was used as a computational follow-up to the LUAD drug-repurposing screen. Palbociclib was evaluated against CDK6 using a known inhibitor-bound structure and against CDK4 using a separate inhibitor-bound structure.

## Software

- AutoDock Vina
- Meeko
- Gemmi
- RDKit
- Open Babel

Record the exact installed software versions in the final reproducibility record.

## CDK6 redocking validation

The human CDK6 structure PDB 5L2I contains palbociclib in the experimentally observed binding site. The crystal ligand coordinates were used to prepare the ligand, and the receptor was prepared for docking with Meeko.

The Vina search used:
- Exhaustiveness: 16
- Maximum poses: 9
- Random seed: 2026
- Search box centered on the native palbociclib ligand

The best-ranked pose had a Vina score of -10.640 kcal/mol and a reported crystal-pose RMSD of 0.646 Å. Two of nine poses had RMSD <= 2 Å.

These results are encouraging for pose recovery, but the atom mapping, protonation, receptor preparation, and RMSD calculation should be documented and independently checked before making a formal validation claim.

## CDK4 exploratory cross-docking

The human CDK4 structure PDB 9CSK contains the inhibitor atirmociclib. CDK4 chain D was extracted and prepared with Meeko. The docking box was centered on the native atirmociclib binding site.

Palbociclib was docked into the prepared CDK4 receptor using:
- Exhaustiveness: 16
- Maximum poses: 9
- Random seed: 2026

The best-ranked pose had a Vina score of -9.164 kcal/mol. The initial heavy-atom contact screen found no ligand–receptor heavy-atom pairs below 2.5 Å.

The CDK4 experiment is exploratory cross-docking. It is not a redocking validation because the experimental CDK4 structure used here does not contain palbociclib.

## Outputs

- `results/tables/CDK6_palbo_redocking_validation.csv`
- `results/tables/CDK4_palbo_crossdocking_scores.csv`
- `results/tables/CDK4_palbo_pose1_residue_contacts.csv`
- `results/tables/CDK4_palbo_pose1_closest_atom_pairs.csv`
- `results/tables/Palbociclib_docking_summary.csv`
- `data/processed/docking/CDK6_5L2I_receptor.pdbqt`
- `data/processed/docking/palbociclib_native.pdbqt`
- `data/processed/docking/palbociclib_redocked.pdbqt`
- `data/processed/docking/CDK4_9CSK_receptor.pdbqt`
- `data/processed/docking/CDK4_palbo_crossdocked.pdbqt`

## Limitations

1. Vina scores are approximate scoring-function outputs, not experimental binding free energies.
2. The CDK4 result is exploratory and cannot be directly compared with the CDK6 result as though the receptor structures were equivalent.
3. The CDK4 receptor preparation used the isolated kinase chain and omitted cyclin D1.
4. A lack of severe steric overlaps does not establish a biologically correct binding pose.
5. Docking does not demonstrate drug sensitivity, clinical efficacy, or safety in lung adenocarcinoma.
6. Follow-up work should inspect hydrogen bonds, ligand strain, protonation states, binding-site geometry, and the effects of receptor preparation.

## Conclusion

The CDK6 experiment provides preliminary evidence that the docking workflow can recover the known palbociclib pose. The CDK4 cross-docking result provides a testable structural hypothesis that requires additional validation. Neither result alone establishes therapeutic efficacy in LUAD.
