from Bio.PDB import MMCIFParser
from pathlib import Path
import numpy as np

pdb_path = Path("data/raw/pdb/5L2I.cif")
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("CDK6_PALBO", str(pdb_path))

model = structure[0]

print("Structure:", pdb_path)
print("\nChains and residue counts:")
for chain in model:
    residues = list(chain.get_residues())
    standard = sum(1 for r in residues if r.id[0] == " ")
    ligands = [
        (r.resname, r.id[1], r.id[0])
        for r in residues if r.id[0] != " "
    ]
    print(
        f"Chain {chain.id!r}: {standard} standard residues; "
        f"nonstandard residues: {ligands}"
    )

print("\nPalbociclib candidates:")
found = False
for chain in model:
    for residue in chain:
        if residue.resname == "LQQ":
            atoms = [a for a in residue.get_atoms() if a.element != "H"]
            coords = np.array([a.coord for a in atoms])
            center = coords.mean(axis=0)
            print(f"Chain: {chain.id}")
            print(f"Residue ID: {residue.id}")
            print(f"Heavy atoms: {len(atoms)}")
            print(
                "Ligand centroid (Å): "
                + ", ".join(f"{x:.3f}" for x in center)
            )
            found = True

if not found:
    print("LQQ was not found; inspect the nonstandard residue list above.")
