import gemmi

path = "data/raw/pdb/9CSK.cif"
structure = gemmi.make_structure_from_block(
    gemmi.cif.read(str(path)).sole_block()
)
chain = structure[0]["B"]

targets = {14, 35, 235, 240}
center = gemmi.Position(8.258, 0.987, 36.964)

for residue in chain:
    if residue.seqid.num not in targets:
        continue

    print(f"\nResidue {residue.seqid.num} {residue.name}")
    heavy_atoms = [
        atom for atom in residue
        if atom.element.name not in ("H", "D")
    ]

    for atom in heavy_atoms:
        print(
            f"  {atom.name:>4s} "
            f"{atom.pos.x:8.3f} {atom.pos.y:8.3f} {atom.pos.z:8.3f}"
        )

    min_dist = min(atom.pos.dist(center) for atom in heavy_atoms)
    print(f"  Closest heavy atom to ligand centroid: {min_dist:.2f} Å")
