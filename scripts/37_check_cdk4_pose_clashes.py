from pathlib import Path
import math
import csv

base = Path("data/processed/docking")
pose_file = base / "CDK4_palbo_crossdocked.pdbqt"
receptor_file = base / "CDK4_9CSK_receptor.pdbqt"

def read_atoms(path, first_model=False):
    lines = Path(path).read_text().splitlines()
    atoms = []
    inside = not first_model

    for line in lines:
        if first_model and line.startswith("MODEL"):
            if atoms:
                break
            inside = True
            continue
        if first_model and line.startswith("ENDMDL") and inside:
            break

        if inside and line.startswith(("ATOM  ", "HETATM")):
            fields = line.split()
            atom_type = fields[-1]
            # Ignore hydrogen atoms in this heavy-atom distance screen.
            if atom_type.upper().startswith("H"):
                continue
            atoms.append({
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "resid": line[22:26].strip(),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "type": atom_type,
            })
    return atoms

ligand = read_atoms(pose_file, first_model=True)
protein = read_atoms(receptor_file)

pairs = []
for la in ligand:
    for pa in protein:
        d = math.dist(
            (la["x"], la["y"], la["z"]),
            (pa["x"], pa["y"], pa["z"])
        )
        pairs.append((d, la, pa))

pairs.sort(key=lambda x: x[0])

out = Path("results/tables/CDK4_palbo_pose1_closest_atom_pairs.csv")
out.parent.mkdir(parents=True, exist_ok=True)

with out.open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "distance_A", "ligand_atom", "ligand_type",
        "protein_chain", "protein_residue", "protein_resname",
        "protein_atom", "protein_type"
    ])
    for d, la, pa in pairs:
        if d > 4.0:
            break
        writer.writerow([
            round(d, 3), la["name"], la["type"],
            pa["chain"], pa["resid"], pa["resname"],
            pa["name"], pa["type"]
        ])

print("Ligand heavy atoms:", len(ligand))
print("Receptor heavy atoms:", len(protein))
print("Atom pairs below 2.0 Å:", sum(d < 2.0 for d, _, _ in pairs))
print("Atom pairs below 2.5 Å:", sum(d < 2.5 for d, _, _ in pairs))
print("Saved:", out)
print("\n10 closest atom pairs:")
for d, la, pa in pairs[:10]:
    print(
        f"{d:5.2f} Å | ligand {la['name']:>4s} ({la['type']}) "
        f"-- {pa['resname']} {pa['resid']} {pa['name']} ({pa['type']})"
    )
