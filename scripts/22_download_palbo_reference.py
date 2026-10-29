from pathlib import Path
import urllib.request

outdir = Path("data/processed/docking")
outdir.mkdir(parents=True, exist_ok=True)

urls = {
    "palbociclib_RCSB_ideal.sdf":
        "https://files.rcsb.org/ligands/download/LQQ_ideal.sdf",
    "palbociclib_crystal_instance.sdf":
        "https://models.rcsb.org/v1/5l2i/ligand?auth_asym_id=A&auth_seq_id=900&encoding=sdf",
}

for filename, url in urls.items():
    destination = outdir / filename
    print(f"\nDownloading {filename}")

    try:
        request = urllib.request.Request(
            url, headers={"User-Agent": "LUAD-drug-repurposing-project"}
        )
        with urllib.request.urlopen(request, timeout=60) as response:
            content = response.read()

        if len(content) < 100:
            raise ValueError("Response unexpectedly small.")

        destination.write_bytes(content)
        print(f"Saved: {destination} ({len(content):,} bytes)")

    except Exception as exc:
        print(f"FAILED: {exc}")
