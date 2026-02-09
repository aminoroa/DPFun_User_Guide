import os
import pickle as pkl
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

PID_LIST_PKL = "./data/test_pid_list.pkl"
PDB_DIR = "./data/PDB/PDB_folder"
OUT_POINTS = "./processed_file/pdb_points.pkl"

def read_pkl(fp):
    with open(fp, "rb") as f:
        return pkl.load(f)

def save_pkl(fp, obj):
    os.makedirs(os.path.dirname(fp), exist_ok=True)
    with open(fp, "wb") as f:
        pkl.dump(obj, f)

def extract_points_best_chain(pdb_file, preferred=("A","H","L","h","l")):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    chain_points = {}
    for model in structure:
        for chain in model:
            pts = []
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                try:
                    _ = seq1(residue.resname)  # standard AA
                except KeyError:
                    continue
                if "CA" not in residue:
                    continue
                c = residue["CA"].get_coord()
                pts.append((float(c[0]), float(c[1]), float(c[2])))
            if pts:
                chain_points[chain.id] = pts

    for cid in preferred:
        if cid in chain_points:
            return cid, chain_points[cid]

    if chain_points:
        cid, pts = max(chain_points.items(), key=lambda x: len(x[1]))
        return cid, pts

    return None, []

pid_list = read_pkl(PID_LIST_PKL)
points_info = {}

for pid in pid_list:
    pdb_path = os.path.join(PDB_DIR, f"{pid}.pdb")
    if not os.path.exists(pdb_path):
        raise FileNotFoundError(f"Missing PDB: {pdb_path}")

    chain_id, pts = extract_points_best_chain(pdb_path)
    if not pts:
        raise ValueError(f"Could not extract CA points from ANY chain in {pdb_path}")

    points_info[pid] = pts
    print(pid, f"chain={chain_id}", "points=", len(pts))

save_pkl(OUT_POINTS, points_info)
print("Saved:", OUT_POINTS, "proteins:", len(points_info))
