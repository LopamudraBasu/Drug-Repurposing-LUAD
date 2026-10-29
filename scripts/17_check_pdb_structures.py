import json
import urllib.request

pdb_ids = ["9CSK", "2W99", "5L2I", "7JUR", "7JUX"]

for pdb_id in pdb_ids:
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"

    print(f"\n{'=' * 65}\nPDB: {pdb_id}")

    try:
        req = urllib.request.Request(
            url, headers={"User-Agent": "LUAD-drug-repurposing-project"}
        )
        with urllib.request.urlopen(req, timeout=30) as response:
            data = json.load(response)

        info = data.get("rcsb_entry_info", {})
        print("Title:", data.get("struct", {}).get("title", "Unknown"))
        print("Experimental method:", info.get("experimental_method"))
        print("Resolution (Å):", info.get("resolution_combined"))
        print("URL:", f"https://www.rcsb.org/structure/{pdb_id}")

    except Exception as e:
        print("Could not retrieve metadata:", e)
