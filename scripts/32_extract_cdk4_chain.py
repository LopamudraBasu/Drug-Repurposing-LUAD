from pathlib import Path
import gemmi

source = Path("data/raw/pdb/9CSK.cif")
outdir = Path("data/processed/docking")
outdir.mkdir(parents=True, exist_ok=True)

structure = gemmi.make_structure_from_block(
    gemmi.cif.read(str(source)).sole_block()
)
model = structure[0]

# Keep CDK4 chain B and exclude non-polymer residues.
chain = model["B"]
new_structure = gemmi.Structure()
new_structure.name = "CDK4_chain_B_9CSK"
new_model = gemmi.Model("1")
new_chain = gemmi.Chain("B")

for residue in chain:
    if residue.entity_type == gemmi.EntityType.Polymer:
        new_chain.add_residue(residue.clone())

new_model.add_chain(new_chain)
new_structure.add_model(new_model)
new_structure.setup_entities()

output = outdir / "CDK4_9CSK_chainB.pdb"
new_structure.write_pdb(str(output))

print(f"Saved: {output}")
print(f"Polymer residues extracted: {len(new_chain)}")
print(f"Output size: {output.stat().st_size:,} bytes")
