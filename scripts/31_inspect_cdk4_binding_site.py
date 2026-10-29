from pathlib import Path
import gemmi

path = Path("data/raw/pdb/9CSK.cif")
structure = gemmi.make_structure_from_block(
    gemmi.cif.read(str(path)).sole_block()
)
model = structure[0]

# Inspect the first CDK4 chain and its nearby native inhibitor.
chain = model["B"]
ligand_chain = model["B"]

ligands = [
    residue for residue in ligand_chain
    if residue.name == "A1AZ4"
]
if not ligands:
    # In mmCIF, the ligand may be represented under a separate chain.
    ligands = [
        residue
        for other_chain in model
        for residue in other_chain
        if residue.name == "A1AZ4"
    ]

if not ligands:
    raise RuntimeError("Could not find the native atirmociclib ligand.")

ligand = ligands[0]
ligand_atoms = [
    atom for atom in ligand
    if atom.element.name not in ("H", "D")
]

print(f"CDK4 chain: {chain.name}")
print(f"CDK4 residues present: {sum(1 for r in chain if r.entity_type == gemmi.EntityType.Polymer)}")
print(f"Native ligand: {ligand.name}, chain {ligand.seqid}")
print(f"Native ligand heavy atoms: {len(ligand_atoms)}")

# Calculate minimum distance from each CDK4 residue to the ligand.
nearby = []
for residue in chain:
    if residue.entity_type != gemmi.EntityType.Polymer:
        continue
    min_distance = min(
        atom.pos.dist(lig_atom.pos)
        for atom in residue
        if atom.element.name not in ("H", "D")
        for lig_atom in ligand_atoms
    )
    if min_distance <= 5.0:
        nearby.append((residue.seqid.num, residue.name, min_distance))

print("\nCDK4 residues within 5 Å of atirmociclib:")
for number, name, distance in sorted(nearby):
    print(f"{number:4d} {name:3s} {distance:5.2f} Å")

xs = [a.pos.x for a in ligand_atoms]
ys = [a.pos.y for a in ligand_atoms]
zs = [a.pos.z for a in ligand_atoms]
center = (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))
print("\nNative ligand centroid (Å):")
print(f"x={center[0]:.3f}, y={center[1]:.3f}, z={center[2]:.3f}")
