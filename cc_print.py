import pickle as pkl
from goatools.obo_parser import GODag

df = pkl.load(open("./results/DPFunc_model_cc_final.pkl", "rb"))
godag = GODag("./data/go.obo")

TOP_K = 10

for _, row in df.iterrows():
    pid = row["protein_id"]
    preds = row["predictions"]

    top = sorted(preds.items(), key=lambda x: x[1], reverse=True)[:TOP_K]

    print("=" * 90)
    print(f"Protein: {pid}")
    print("-" * 90)
    for go, score in top:
        name = godag[go].name if go in godag else "NA"
        print(f"{go:>12s}  score={score:.4f}  {name}")