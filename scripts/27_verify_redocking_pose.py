from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

folder = Path("data/processed/docking")

reference = next(
    m for m in Chem.SDMolSupplier(
        str(folder / "palbociclib_crystal_instance.sdf"),
        removeHs=False
    ) if m is not None
)

print("Reference heavy atoms:", reference.GetNumHeavyAtoms())
print("Reference formula:", rdMolDescriptors.CalcMolFormula(reference))

docked_file = folder / "palbociclib_redocked.pdbqt"
docked_sdf = folder / "palbociclib_redocked.sdf"

supplier = Chem.SDMolSupplier(str(docked_sdf), removeHs=False)
print("\nDocked poses:")

for i, pose in enumerate(supplier, start=1):
    if pose is None:
        print(f"Pose {i}: could not parse")
        continue

    print(
        f"Pose {i}: atoms={pose.GetNumAtoms()}, "
        f"heavy_atoms={pose.GetNumHeavyAtoms()}, "
        f"formula={rdMolDescriptors.CalcMolFormula(pose)}"
    )
