from pathlib import Path
import gemmi

source = Path("data/raw/pdb/9CSK.cif")
outdir = Path("data/processed/docking")
outdir.mkdir(parents=True, exist_ok=True)

structure = gemmi.make_structure_from_block(
    gemmi.cif.read(str(source)).sole_block()
)
chain = structure[0]["D"]

new_structure = gemmi.Structure()
new_structure.name = "CDK4_9CSK_chainD"
new_model = gemmi.Model("1")
new_chain = gemmi.Chain("D")

for residue in chain:
    if residue.entity_type == gemmi.EntityType.Polymer:
        new_chain.add_residue(residue.clone())

new_model.add_chain(new_chain)
new_structure.add_model(new_model)
new_structure.setup_entities()

output = outdir / "CDK4_9CSK_chainD.pdb"
new_structure.write_pdb(str(output))

print(f"Saved: {output}")
print(f"Polymer residues: {len(new_chain)}")
print(f"File size: {output.stat().st_size:,} bytes")
