import csv
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['figure.constrained_layout.use'] = True

# with open("./results/mcm_ids.csv", "r") as file:
#     data = csv.reader(file.readlines())
#     header = data.__next__()
#     mcm_history: list[list[int]] = [[] for i in range(len(header))]
#     for line in data:
#         for i, mcm in zip(range(len(header)), line):
#             # mcm_history[int(mcm)].append(i)
#             mcm_history[int(mcm)].append(-i)

# fig, axes = plt.subplots(len(header), 1)
# for mcm, data in zip(range(len(header)), mcm_history):
#     if not mcm == len(header) - 1:
#         axes[mcm].set(xticks=[])
#     axes[mcm].plot(range(len(data)), data, "-", label=f"{mcm}")
# plt.show()

chunk_size = 1200
with open("./results/mutation_history.csv", "r") as file:
    data = csv.reader(file.readlines())
    header = data.__next__()

    mutation_history: list[dict[str, int]] = []
    mutation_accepted_history: list[dict[str, int]] = []
    counter = 0
    for line in data: 
        if not counter:
            mutation_history.append({}.fromkeys(["Swap", "Split", "Merge"], 0))
            mutation_accepted_history.append({}.fromkeys(["Swap", "Split", "Merge"], 0))

        mutation_history[-1][line[0]] += 1
        if line[1] == "true":
            mutation_accepted_history[-1][line[0]] += 1
            
        counter += 1
        if counter == chunk_size: 
            counter = 0

plot_data: dict[str, list[int]] = {"Merge": [], "Swap": [], "Split": [], }
plot_accepted_data: dict[str, list[int]] = {"Merge": [], "Swap": [], "Split": [], }
btm = []

for chunk in mutation_history:
    for key, value in chunk.items():
        plot_data[key].append(value)
    btm.append(0)
for chunk in mutation_accepted_history:
    for key, value in chunk.items():
        plot_accepted_data[key].append(value)

bottom = np.array(btm)
# bottom_accepted = np.array(bottom)

fig, axes = plt.subplots(2, 1)
for merge_type, counts in plot_data.items():
    axes[0].bar(range(len(counts)), counts, 1.0, label=merge_type, bottom=bottom)  
    bottom += np.array(counts)
axes[0].legend()
bottom[:] = 0

for merge_type, counts in plot_accepted_data.items():
    axes[1].bar(range(len(counts)), counts, 1.0, label=merge_type, bottom=bottom)  
    bottom += np.array(counts)
axes[1].legend()
plt.show()