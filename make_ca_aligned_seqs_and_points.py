import os
import pickle as pkl
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

PID_LIST_PKL = "./data/test_pid_list.pkl"
PDB_DIR = "./data/PDB/PDB_folder"

OUT_SEQS = "./processed_file/pdb_seqs.pkl"
OUT_POINTS = "./processed_file/pdb_points.pkl"

PREFERRED_CHAINS = ("A", "H", "L", "h", "l")

def read_pkl(fp):
    with open(fp, "rb") as f:
        return pkl.load(f)

def save_pkl(fp, obj):
    os.makedirs(os.path.dirname(fp), exist_ok=True)
    with open(fp, "wb") as f:
        pkl.dump(obj, f)

def extract_chain_ca_aligned(chain):
    """
    Returns (sequence, points) where:
      - sequence contains only standard amino acids WITH CA atoms
      - points is list of (x,y,z) CA coordinates aligned to sequence positions
    """
    seq = []
    pts = []

    for residue in chain:
        if residue.id[0] != " ":
            continue
        if "CA" not in residue:
            continue
        try:
            aa = seq1(residue.resname)  # standard AA
        except KeyError:
            continue

        c = residue["CA"].get_coord()
        seq.append(aa)
        pts.append((float(c[0]), float(c[1]), float(c[2])))

    return "".join(seq), pts

def choose_best_chain(structure):
    """
    Choose a chain deterministically:
      1) preferred chain IDs if present and non-empty
      2) otherwise the chain with the longest CA-aligned sequence
    Returns (chain_id, seq, points)
    """
    # collect all chains (first model only is fine for typical PDBs)
    model = next(iter(structure))

    # try preferred chains first
    chain_map = {ch.id: ch for ch in model}
    for cid in PREFERRED_CHAINS:
        if cid in chain_map:
            seq, pts = extract_chain_ca_aligned(chain_map[cid])
            if seq and pts:
                return cid, seq, pts

    # fallback: longest CA-aligned chain
    best = (None, "", [])
    for ch in model:
        seq, pts = extract_chain_ca_aligned(ch)
        if len(seq) > len(best[1]):
            best = (ch.id, seq, pts)

    return best

def main():
    pid_list = read_pkl(PID_LIST_PKL)

    pdb_seqs = {}
    pdb_points = {}

    for pid in pid_list:
        pdb_path = os.path.join(PDB_DIR, f"{pid}.pdb")
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"Missing PDB: {pdb_path}")

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pid, pdb_path)

        cid, seq, pts = choose_best_chain(structure)
        if not seq or not pts:
            raise ValueError(f"Could not extract CA-aligned residues from any chain in {pdb_path}")

        # By construction these must match:
        assert len(seq) == len(pts), f"Internal error: seq/pts mismatch for {pid}"

        pdb_seqs[pid] = seq
        pdb_points[pid] = pts
        print(pid, f"chain={cid}", f"L={len(seq)}")

    save_pkl(OUT_SEQS, pdb_seqs)
    save_pkl(OUT_POINTS, pdb_points)
    print("Saved:", OUT_SEQS)
    print("Saved:", OUT_POINTS)

if __name__ == "__main__":
    main()
