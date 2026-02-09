import pickle as pkl, os

pdb_seqs = pkl.load(open("./processed_file/pdb_seqs.pkl","rb"))
out_fa = "./processed_file/myproteins.fasta"

with open(out_fa, "w") as f:
    for pid, seq in pdb_seqs.items():
        f.write(f">{pid}\n")
        # wrap lines to 60 chars (safe for tools)
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

print("Wrote:", out_fa, "proteins:", len(pdb_seqs))