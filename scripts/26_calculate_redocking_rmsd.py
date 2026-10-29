from pathlib import Path
import numpy as np

folder = Path("data/processed/docking")
reference_file = folder / "palbociclib_native.pdbqt"
docked_file = folder / "palbociclib_redocked.pdbqt"


def read_pdbqt_models(path):
    models = {}
    current_model = 1
    models[current_model] = []

    for line in path.read_text().splitlines():
        if line.startswith("MODEL"):
            current_model = int(line.split()[1])
            models.setdefault(current_model, [])
        elif line.startswith(("ATOM  ", "HETATM")):
            atom_type = line[77:79].strip() if len(line) >= 79 else line.split()[-1]
            if atom_type.upper() in {"H", "HD"}:
                continue

            try:
                xyz = [
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ]
            except ValueError:
                continue

            models[current_model].append(xyz)

    return {
        model_id: np.asarray(coords, dtype=float)
        for model_id, coords in models.items()
        if coords
    }


reference_models = read_pdbqt_models(reference_file)
docked_models = read_pdbqt_models(docked_file)

reference = reference_models[1]

print("Reference heavy atoms:", len(reference))
print("Docked poses:", len(docked_models))
print("\nPose\tHeavy atoms\tDirect RMSD (Å)")

for model_id, coords in sorted(docked_models.items()):
    if len(coords) != len(reference):
        print(f"{model_id}\t{len(coords)}\tAtom count mismatch; not calculated")
        continue

    # The PDBQT preparation preserves atom order across the poses.
    # This is a direct coordinate comparison, without ligand alignment.
    rmsd = np.sqrt(np.mean(np.sum((coords - reference) ** 2, axis=1)))
    print(f"{model_id}\t{len(coords)}\t{rmsd:.3f}")

print("\nNote: direct RMSD assumes identical atom ordering in the")
print("reference and docked PDBQT files. Verify atom mapping before")
print("using these values as a formal validation result.")
