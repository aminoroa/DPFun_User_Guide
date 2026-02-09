import pickle as pkl
e0 = pkl.load(open("./processed_file/esm_emds/esm_part_0.pkl","rb"))
points = pkl.load(open("./processed_file/pdb_points.pkl","rb"))

mismatches = []
for pid in sorted(e0.keys()):
    L_esm = e0[pid].shape[0]
    L_pts = len(points[pid])
    ok = (L_esm == L_pts)
    print(f"{pid}: ESM={L_esm} points={L_pts} {'OK' if ok else 'MISMATCH'}")
    if not ok:
        mismatches.append(pid)

print("Mismatches:", mismatches)