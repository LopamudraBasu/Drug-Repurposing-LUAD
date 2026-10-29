import gemmi

path = "data/raw/pdb/9CSK.cif"
structure = gemmi.make_structure_from_block(
    gemmi.cif.read(str(path)).sole_block()
)
model = structure[0]

for chain_id in ("B", "D"):
    chain = model[chain_id]
    ligand = next(
        r for r in chain
        if r.name == "A1AZ4"
    ) if any(r.name == "A1AZ4" for r in chain) else None

    # The ligand may be in a separate chain; locate the copy
    # whose coordinates are closest to this CDK4 chain.
    all_ligands = [
        r for c in model for r in c
        if r.name == "A1AZ4"
    ]
    if ligand is None:
        protein_atoms = [
            a for r in chain for a in r
            if a.element.name not in ("H", "D")
        ]
        ligand = min(
            all_ligands,
            key=lambda lig: min(
                a.pos.dist(b.pos)
                for a in protein_atoms
                for b in lig
                if b.element.name not in ("H", "D")
            )
        )

    ligand_atoms = [
        a for a in ligand
        if a.element.name not in ("H", "D")
    ]
    center = gemmi.Position(
        sum(a.pos.x for a in ligand_atoms) / len(ligand_atoms),
        sum(a.pos.y for a in ligand_atoms) / len(ligand_atoms),
        sum(a.pos.z for a in ligand_atoms) / len(ligand_atoms),
    )

    print(f"\n=== CDK4 chain {chain_id} ===")
    print(f"Residues: {sum(r.entity_type == gemmi.EntityType.Polymer for r in chain)}")
    print(f"Associated ligand centroid: {center.x:.3f}, {center.y:.3f}, {center.z:.3f}")

    for number in (12, 14, 20, 35, 93, 95, 96, 158, 235, 240):
        residue = next(
            (r for r in chain if r.seqid.num == number), None
        )
        if residue is None:
            print(f"Residue {number}: absent")
            continue

        atoms = [
            a for a in residue
            if a.element.name not in ("H", "D")
        ]
        names = {a.name for a in atoms}
        distance = min(a.pos.dist(center) for a in atoms) if atoms else float("nan")
        print(
            f"{number:3d} {residue.name:3s}: "
            f"{len(names):2d} heavy atoms; "
            f"nearest to ligand centroid {distance:.2f} Å"
        )
