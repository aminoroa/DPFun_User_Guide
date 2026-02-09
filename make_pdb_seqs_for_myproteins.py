import os
import pickle as pkl
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

PDB_DIR = "./data/PDB/PDB_folder"
PID_LIST_PKL = "./data/test_pid_list.pkl"
OUT_SEQS = "./processed_file/pdb_seqs.pkl"

def read_pkl(fp):
    with open(fp, "rb") as f:
        return pkl.load(f)

def save_pkl(fp, obj):
    os.makedirs(os.path.dirname(fp), exist_ok=True)
    with open(fp, "wb") as f:
        pkl.dump(obj, f)

def extract_best_chain_sequence(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    best_chain = None
    best_seq = ""

    for model in structure:
        for chain in model:
            seq = []
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                try:
                    seq.append(seq1(residue.resname))
                except KeyError:
                    # non-standard residue; skip
                    pass
            seq = "".join(seq)
            if len(seq) > len(best_seq):
                best_seq = seq
                best_chain = chain.id

    return best_chain, best_seq

pid_list = read_pkl(PID_LIST_PKL)

pdb_seqs = {}
for pid in pid_list:
    pdb_path = os.path.join(PDB_DIR, f"{pid}.pdb")
    if not os.path.exists(pdb_path):
        raise FileNotFoundError(f"Missing PDB: {pdb_path}")
    chain_id, seq = extract_best_chain_sequence(pdb_path)
    if not seq:
        raise ValueError(f"Could not extract sequence from ANY chain in {pdb_path}")
    pdb_seqs[pid] = seq
    print(pid, f"chain={chain_id}", len(seq))
    print(pid, len(seq))

save_pkl(OUT_SEQS, pdb_seqs)
print("Saved:", OUT_SEQS, "proteins:", len(pdb_seqs))
