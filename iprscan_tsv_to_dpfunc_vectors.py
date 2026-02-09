import os
import pickle as pkl
import numpy as np

PID_LIST = pkl.load(open("./data/test_pid_list.pkl","rb"))
IPRSCAN_TSV = "./processed_file/interproscan_out/myproteins.interproscan.tsv"

# EDIT THIS to the discovered 22369-list path:
IPR_LIST_PKL = "./PATH/TO/interpro_list_22369.pkl"

ipr_list = pkl.load(open(IPR_LIST_PKL, "rb"))
ipr_idx = {ipr:i for i, ipr in enumerate(ipr_list)}
N = len(ipr_list)

# parse TSV -> pid -> set(IPR)
hits = {pid:set() for pid in PID_LIST}
with open(IPRSCAN_TSV, "r") as f:
    for line in f:
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 12:
            continue
        pid = cols[0]
        ipr = cols[11]
        if pid in hits and ipr and ipr != "-" and ipr in ipr_idx:
            hits[pid].add(ipr)

out_dir = "./processed_file/interpro_22369"
os.makedirs(out_dir, exist_ok=True)

X = np.zeros((len(PID_LIST), N), dtype=np.float32)
for i, pid in enumerate(PID_LIST):
    vec = np.zeros(N, dtype=np.float32)
    for ipr in hits[pid]:
        vec[ipr_idx[ipr]] += 1.0
    pkl.dump(vec, open(f"{out_dir}/{pid}.pkl","wb"))
    X[i,:] = vec

out_mat = "./processed_file/mf_mytest_interpro_22369.pkl"
pkl.dump(X, open(out_mat, "wb"))

print("Wrote per-protein vectors to:", out_dir)
print("Wrote matrix to:", out_mat, "shape:", X.shape)