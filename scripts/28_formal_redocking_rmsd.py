from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolAlign, rdMolDescriptors

folder = Path("data/processed/docking")

reference = next(
    m for m in Chem.SDMolSupplier(
        str(folder / "palbociclib_crystal_instance.sdf"),
        removeHs=False
    ) if m is not None
)

print("Reference formula:", rdMolDescriptors.CalcMolFormula(reference))
print("Reference heavy atoms:", reference.GetNumHeavyAtoms())
print("\nPose\tHeavy atoms\tSymmetry-aware RMSD (Å)")

for i in range(1, 10):
    path = folder / f"palbociclib_redocked{i}.sdf"
    if not path.exists():
        print(f"{i}\tMISSING\t-")
        continue

    pose = next(
        (m for m in Chem.SDMolSupplier(str(path), removeHs=False)
         if m is not None),
        None
    )
    if pose is None:
        print(f"{i}\tPARSE ERROR\t-")
        continue

    # PDBQT-to-SDF conversion can lose or alter bond-order information.
    # Assign the known palbociclib connectivity before atom mapping.
    # The atom order in the converted ligand is not assumed to match
    # the reference; use substructure matching on the molecular graph.
    ref_no_h = Chem.RemoveHs(reference)
    pose_no_h = Chem.RemoveHs(pose)

    if ref_no_h.GetNumHeavyAtoms() != pose_no_h.GetNumHeavyAtoms():
        print(f"{i}\t{pose_no_h.GetNumHeavyAtoms()}\tAtom-count mismatch")
        continue

    matches = ref_no_h.GetSubstructMatches(pose_no_h, uniquify=False,
                                           maxMatches=10000)
    if not matches:
        print(f"{i}\t{pose_no_h.GetNumHeavyAtoms()}\tNo graph match")
        continue

    # Each match maps pose atoms to reference atoms.
    # Compute RMSD without independently superposing coordinates.
    import numpy as np
    ref_conf = ref_no_h.GetConformer()
    pose_conf = pose_no_h.GetConformer()
    best = float("inf")

    for match in matches:
        # match[k] is the reference atom corresponding to pose atom k.
        ref_xyz = np.array([
            ref_conf.GetAtomPosition(match[k])
            for k in range(pose_no_h.GetNumAtoms())
        ])
        pose_xyz = np.array([
            pose_conf.GetAtomPosition(k)
            for k in range(pose_no_h.GetNumAtoms())
        ])
        rmsd = np.sqrt(np.mean(np.sum((pose_xyz - ref_xyz) ** 2, axis=1)))
        best = min(best, rmsd)

    print(f"{i}\t{pose_no_h.GetNumHeavyAtoms()}\t{best:.3f}")
