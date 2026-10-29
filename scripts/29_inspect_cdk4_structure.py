from pathlib import Path
import gemmi

path = Path("data/raw/pdb/9CSK.cif")
if not path.exists():
    raise FileNotFoundError(f"Missing structure: {path}")

doc = gemmi.cif.read(str(path))
block = doc.sole_block()

print("=== Entity descriptions ===")
entity_ids = block.find_values("_entity.id")
descriptions = block.find_values("_entity.pdbx_description")
for entity_id, description in zip(entity_ids, descriptions):
    print(f"Entity {entity_id}: {description}")

print("\n=== Polymer chains ===")
structure = gemmi.make_structure_from_block(block)
model = structure[0]

for chain in model:
    residues = [r for r in chain if r.entity_type == gemmi.EntityType.Polymer]
    if residues:
        standard = [
            r for r in residues
            if r.name in gemmi.find_tabulated_residue(r.name).one_letter_code
        ]
        print(
            f"Chain {chain.name}: polymer residues={len(residues)}, "
            f"first={residues[0].name if residues else 'NA'}, "
            f"last={residues[-1].name if residues else 'NA'}"
        )

print("\n=== Non-polymer residues / ligands ===")
for chain in model:
    for residue in chain:
        if residue.entity_type == gemmi.EntityType.NonPolymer:
            atoms = [a for a in residue if a.element.name not in ("H", "D")]
            if atoms:
                xs = [a.pos.x for a in atoms]
                ys = [a.pos.y for a in atoms]
                zs = [a.pos.z for a in atoms]
                center = (
                    sum(xs) / len(xs),
                    sum(ys) / len(ys),
                    sum(zs) / len(zs),
                )
                print(
                    f"Chain {chain.name}, residue {residue.name} "
                    f"{residue.seqid}: heavy_atoms={len(atoms)}, "
                    f"center=({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})"
                )
