from pathlib import Path
import urllib.request

outdir = Path("data/raw/pdb")
outdir.mkdir(parents=True, exist_ok=True)

pdb_ids = ["9CSK", "2W99", "5L2I"]

for pdb_id in pdb_ids:
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    destination = outdir / f"{pdb_id}.cif"

    print(f"\nDownloading {pdb_id}...")
    try:
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "LUAD-drug-repurposing-project"}
        )
        with urllib.request.urlopen(request, timeout=60) as response:
            content = response.read()

        if len(content) < 1000:
            raise ValueError("Downloaded file is unexpectedly small.")

        destination.write_bytes(content)
        print(f"Saved: {destination} ({len(content):,} bytes)")

    except Exception as exc:
        print(f"FAILED: {pdb_id}: {exc}")
