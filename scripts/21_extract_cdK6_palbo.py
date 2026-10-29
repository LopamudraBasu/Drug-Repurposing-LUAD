from Bio.PDB import MMCIFParser, PDBIO, Select
from pathlib import Path

source = Path("data/raw/pdb/5L2I.cif")
outdir = Path("data/processed/docking")
structure = MMCIFParser(QUIET=True).get_structure("CDK6", str(source))

class ProteinSelect(Select):
    def accept_chain(self, chain):
        return chain.id == "A"

    def accept_residue(self, residue):
        return residue.id[0] == " "  # standard protein residues only

class LigandSelect(Select):
    def accept_chain(self, chain):
        return chain.id == "A"

    def accept_residue(self, residue):
        return residue.resname == "LQQ"

io = PDBIO()
io.set_structure(structure)
io.save(str(outdir / "CDK6_5L2I_protein.pdb"), ProteinSelect())
io.save(str(outdir / "palbociclib_native_5L2I.pdb"), LigandSelect())

for name in ["CDK6_5L2I_protein.pdb", "palbociclib_native_5L2I.pdb"]:
    path = outdir / name
    print(f"{path}: {path.stat().st_size:,} bytes")
    print("  First lines:")
    with path.open() as handle:
        for _, line in zip(range(3), handle):
            print(" ", line.rstrip())
