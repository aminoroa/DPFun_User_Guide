import pickle

pid_list = [
    "3H3B",
    "3KDM",
    "4MN8",
    "8HND",
    "8IKW",
    "8IQS",
    "8JYR"
]

with open("./data/test_pid_list.pkl", "wb") as f:
    pickle.dump(pid_list, f)

print("Saved test_pid_list.pkl with", len(pid_list), "proteins")
