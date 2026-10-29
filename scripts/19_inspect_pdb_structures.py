from pathlib import Path
import collections
import re

pdb_dir = Path("data/raw/pdb")
pdb_ids = ["9CSK", "2W99", "5L2I"]

for pdb_id in pdb_ids:
    path = pdb_dir / f"{pdb_id}.cif"
    text = path.read_text(errors="replace")

    print("\n" + "=" * 75)
    print(f"PDB: {pdb_id}")
    print("=" * 75)

    # Print structure title
    title_match = re.search(
        r"_struct\.title\s+(?:'([^']*)'|\"([^\"]*)\"|([^\n]+))",
        text
    )
    if title_match:
        print("Title:", next(g for g in title_match.groups() if g is not None).strip())

    # Read the non-polymer component table
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.strip() == "_chem_comp.id":
            print("\nNon-polymer chemical components:")
            j = i + 1
            while j < len(lines) and not lines[j].startswith("loop_"):
                row = lines[j].strip()
                if row and not row.startswith("_") and not row.startswith("#"):
                    print(" ", row[:180])
                j += 1
            break

    # Summarize polymer chains from atom-site records
    chain_counts = collections.Counter()
    entity_counts = collections.Counter()

    for line in lines:
        if not line.startswith("ATOM "):
            continue
        fields = line.split()
        # Atom-site records are variable-width; report only safely
        # identifiable chain/entity details in a compact summary.
        if len(fields) > 18:
            chain_counts[fields[18]] += 1

    if chain_counts:
        print("\nAtom-site chain token counts (preliminary; verify before use):")
        for chain, count in chain_counts.most_common():
            print(f"  {chain}: {count:,} atom records")

    print("\nFile:", path, f"({path.stat().st_size:,} bytes)")
