import gemmi

path = "data/raw/pdb/9CSK.cif"
block = gemmi.cif.read(path).sole_block()

print("=== Chain-to-entity mapping ===")
labels = block.find_values("_struct_asym.id")
entities = block.find_values("_struct_asym.entity_id")

entity_desc = dict(zip(
    block.find_values("_entity.id"),
    block.find_values("_entity.pdbx_description")
))

for chain_id, entity_id in zip(labels, entities):
    print(
        f"Chain {chain_id}: entity {entity_id} — "
        f"{entity_desc.get(entity_id, 'description unavailable')}"
    )

print("\n=== Polymer sequence lengths ===")
poly_ids = block.find_values("_entity_poly.entity_id")
sequences = block.find_values("_entity_poly.pdbx_seq_one_letter_code_can")

for entity_id, sequence in zip(poly_ids, sequences):
    clean_sequence = "".join(sequence.split()).replace(";", "")
    print(
        f"Entity {entity_id}: "
        f"{entity_desc.get(entity_id, 'unknown')} — "
        f"{len(clean_sequence)} residues"
    )
