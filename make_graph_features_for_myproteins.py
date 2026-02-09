import os
import math
import pickle as pkl
import torch
import dgl
from tqdm.auto import tqdm

PID_LIST_PKL = "./data/test_pid_list.pkl"

PDB_POINTS_PKL = "./processed_file/pdb_points.pkl"
PROTEIN_MAP_PKL = "./processed_file/protein_map.pkl"
ESM_DIR = "./processed_file/esm_emds"

OUT_GRAPH = "./processed_file/graph_features/mf_test_whole_pdb_part0.pkl"  # name can be mf/cc/bp; structure is same

THRESHOLD = 12.0  # Angstrom, same as tutorial
DTYPE_NODE = torch.float32  # keep consistent

def read_pkl(fp):
    with open(fp, "rb") as f:
        return pkl.load(f)

def save_pkl(fp, obj):
    os.makedirs(os.path.dirname(fp), exist_ok=True)
    with open(fp, "wb") as f:
        pkl.dump(obj, f)

def get_dis(p1, p2):
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    dz = p1[2] - p2[2]
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def main():
    pid_list = read_pkl(PID_LIST_PKL)
    pdb_points = read_pkl(PDB_POINTS_PKL)
    protein_map = read_pkl(PROTEIN_MAP_PKL)

    graphs = []

    for pid in tqdm(pid_list):
        pts = pdb_points[pid]
        n = len(pts)

        # load esm features for this protein
        esm_part_idx = protein_map[pid]
        esm_part_file = os.path.join(ESM_DIR, f"esm_part_{esm_part_idx}.pkl")
        esm_features = read_pkl(esm_part_file)
        node_feat = esm_features[pid]  # numpy array (L, 1280)

        if node_feat.shape[0] != n:
            raise ValueError(f"Length mismatch for {pid}: points={n}, esm={node_feat.shape[0]}")

        # build edges (O(n^2) but fine for 7 proteins)
        u_list, v_list, dis_list = [], [], []
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                d = get_dis(pts[i], pts[j])
                if d <= THRESHOLD:
                    u_list.append(i)
                    v_list.append(j)
                    dis_list.append(d)

        u = torch.tensor(u_list, dtype=torch.int64)
        v = torch.tensor(v_list, dtype=torch.int64)
        dis = torch.tensor(dis_list, dtype=torch.float32)

        g = dgl.graph((u, v), num_nodes=n)
        g.edata["dis"] = dis
        g.ndata["x"] = torch.from_numpy(node_feat).to(DTYPE_NODE)

        graphs.append(g)

        # minimal print
        print(f"{pid}: nodes={n}, edges={g.num_edges()}")

    save_pkl(OUT_GRAPH, graphs)
    print("Saved graphs to:", OUT_GRAPH)
    print("Graph count:", len(graphs))

if __name__ == "__main__":
    main()
