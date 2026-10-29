import csv
import math
import re
from pathlib import Path

base = Path("data/processed/docking")
log_file = base / "CDK4_palbo_crossdocking.log"
pose_file = base / "CDK4_palbo_crossdocked.pdbqt"
receptor_file = base / "CDK4_9CSK_receptor.pdbqt"
outdir = Path("results/tables")
outdir.mkdir(parents=True, exist_ok=True)

for p in (log_file, pose_file, receptor_file):
    if not p.exists():
        raise FileNotFoundError(f"Missing required file: {p}")

# Parse the Vina score table from the log.
scores = []
for line in log_file.read_text().splitlines():
    match = re.match(
        r"^\s*(\d+)\s+(-?\d+\.\d+)\s+"
        r"(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$",
        line
    )
    if match:
        rank, affinity, rmsd_lb, rmsd_ub = match.groups()
        scores.append((int(rank), float(affinity),
                       float(rmsd_lb), float(rmsd_ub)))

if not scores:
    raise RuntimeError("No docking score rows found in the Vina log.")

score_path = outdir / "CDK4_palbo_crossdocking_scores.csv"
with score_path.open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "pose_rank", "vina_affinity_kcal_mol",
        "rmsd_lower_bound_A", "rmsd_upper_bound_A"
    ])
    writer.writerows(scores)

# Read fixed-column PDBQT atom records.
def parse_atom(line):
    return {
        "name": line[12:16].strip(),
        "resname": line[17:20].strip(),
        "chain": line[21:22].strip(),
        "resseq": line[22:26].strip(),
        "x": float(line[30:38]),
        "y": float(line[38:46]),
        "z": float(line[46:54]),
        "type": line.split()[-1],
    }

# Extract ligand atoms from the first ranked pose.
ligand_atoms = []
in_first_model = False

for line in pose_file.read_text().splitlines():
    if line.startswith("MODEL"):
        if ligand_atoms:
            break
        in_first_model = True
        continue
    if line.startswith("ENDMDL"):
        if in_first_model:
            break
    if in_first_model and line.startswith(("ATOM  ", "HETATM")):
        ligand_atoms.append(parse_atom(line))

if not ligand_atoms:
    raise RuntimeError("Could not parse ligand atoms from pose 1.")

receptor_atoms = []
for line in receptor_file.read_text().splitlines():
    if line.startswith(("ATOM  ", "HETATM")):
        receptor_atoms.append(parse_atom(line))

if not receptor_atoms:
    raise RuntimeError("Could not parse receptor atoms.")

# Count ligand–protein heavy-atom contacts and find the closest
# ligand–protein distance for each receptor residue.
residue_stats = {}

for atom in receptor_atoms:
    # Ignore hydrogen atoms if present.
    if atom["type"].upper().startswith("H"):
        continue

    key = (atom["chain"], atom["resseq"], atom["resname"])
    stats = residue_stats.setdefault(
        key, {"min_distance": float("inf"), "contacts_4A": 0}
    )

    for lig in ligand_atoms:
        if lig["type"].upper().startswith("H"):
            continue

        dx = atom["x"] - lig["x"]
        dy = atom["y"] - lig["y"]
        dz = atom["z"] - lig["z"]
        d = math.sqrt(dx*dx + dy*dy + dz*dz)

        if d < stats["min_distance"]:
            stats["min_distance"] = d
        if d <= 4.0:
            stats["contacts_4A"] += 1

# Save residues within 5 Å of the docked ligand.
nearby = [
    (key, stats)
    for key, stats in residue_stats.items()
    if stats["min_distance"] <= 5.0
]
nearby.sort(key=lambda item: item[1]["min_distance"])

contact_path = outdir / "CDK4_palbo_pose1_residue_contacts.csv"
with contact_path.open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "chain", "residue_number", "residue_name",
        "minimum_distance_A", "heavy_atom_contacts_within_4A"
    ])
    for (chain, resseq, resname), stats in nearby:
        writer.writerow([
            chain, resseq, resname,
            round(stats["min_distance"], 3),
            stats["contacts_4A"]
        ])

print(f"Saved docking scores: {score_path}")
print(f"Saved pose-1 contacts: {contact_path}")
print(f"Docking poses recorded: {len(scores)}")
print(f"Top pose affinity: {scores[0][1]:.3f} kcal/mol")
print(f"Ligand atoms parsed in pose 1: {len(ligand_atoms)}")
print("\nNearest CDK4 residues to pose 1:")
print("Chain ResNo ResName MinDist(A) Contacts<=4A")
for (chain, resseq, resname), stats in nearby[:15]:
    print(
        f"{chain or '-':>5s} {resseq:>5s} {resname:>7s} "
        f"{stats['min_distance']:>10.3f} {stats['contacts_4A']:>12d}"
    )
