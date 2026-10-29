
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/portfolio/LINCS_compounds_in_both_screens.csv"
OUTPUT = ROOT / "figures/LINCS_drug_candidate_ranking.png"

df = pd.read_csv(INPUT)

required = {"compound", "combined_evidence_score"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

df = df.dropna(subset=["compound", "combined_evidence_score"])
df["combined_evidence_score"] = pd.to_numeric(
    df["combined_evidence_score"], errors="coerce"
)
df = df.dropna(subset=["combined_evidence_score"])
df = df.sort_values("combined_evidence_score", ascending=False).head(10)
df = df.sort_values("combined_evidence_score", ascending=True)

fig, ax = plt.subplots(figsize=(9, 5.5))
ax.barh(df["compound"], df["combined_evidence_score"])
ax.set_xlabel("Combined evidence score (heuristic)")
ax.set_ylabel("")
ax.set_title("Top 10 LUAD Drug-Repurposing Candidates")
ax.grid(axis="x", alpha=0.25)
ax.set_axisbelow(True)

fig.text(
    0.01,
    0.01,
    "Exploratory LINCS screen; scores are not clinical efficacy estimates.",
    fontsize=8,
)
fig.tight_layout(rect=(0, 0.06, 1, 1))

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUTPUT, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {OUTPUT}")
print(
    df[["compound", "combined_evidence_score"]]
    .sort_values("combined_evidence_score", ascending=False)
    .to_string(index=False)
)
