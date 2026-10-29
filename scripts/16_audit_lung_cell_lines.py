import csv
import re
from collections import defaultdict

files = {
    "LUAD_up_vs_drug_down":
        "results/tables/LINCS_reversal_LUAD_up_vs_drug_down.csv",
    "LUAD_down_vs_drug_up":
        "results/tables/LINCS_reversal_LUAD_down_vs_drug_up.csv",
}

candidates = {
    "palbociclib", "trametinib", "erlotinib", "mitoxantrone",
    "selumetinib", "afatinib", "gefitinib"
}

# Candidate lung cancer cell lines; the output will show which
# actually occur in the source records.
lung_lines = {
    "A549", "HCC515", "H1975", "HCC827", "PC9", "HCC4006",
    "H1650", "H1299", "H358", "H23", "H441", "HCC15",
    "HCC44", "HCC78", "H460", "H322", "H727", "Calu-3",
    "Calu-6", "NCIH1975", "NCIH2228", "NCIH23"
}

for label, filename in files.items():
    print("\n" + "=" * 85)
    print(label)

    grouped = defaultdict(list)

    with open(filename, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    for row in rows:
        term = row.get("term", "")
        # Expected format: plate cell-line duration-compound-dose
        parts = term.split()
        if len(parts) < 3:
            continue

        cell_line = parts[1]
        if cell_line not in lung_lines:
            continue

        # Remove the plate, cell-line and duration prefix.
        compound_part = re.sub(
            r"^[^ ]+ [^ ]+ [0-9]+H-", "", term
        )
        compound = re.sub(r"-[0-9.]+$", "", compound_part).lower()

        if compound in candidates:
            grouped[compound].append({
                "cell_line": cell_line,
                "term": term,
                "adjusted_p_value": row.get("adjusted_p_value", ""),
                "overlapping_genes": row.get("overlapping_genes", "")
            })

    if not grouped:
        print("No candidate records found for the specified lung cell lines.")
        continue

    for drug in sorted(grouped):
        records = grouped[drug]
        records.sort(
            key=lambda r: float(r["adjusted_p_value"])
            if r["adjusted_p_value"] else 1
        )
        print(f"\n{drug}: {len(records)} records")
        print("Cell lines:", ", ".join(sorted(
            {r["cell_line"] for r in records}
        )))

        for r in records[:3]:
            print(
                f"  {r['term']} | adjusted p="
                f"{r['adjusted_p_value']}"
            )
            genes = r["overlapping_genes"].split(";")
            print("  Example overlapping genes:", ", ".join(genes[:12]))
