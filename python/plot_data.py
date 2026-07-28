import csv
import matplotlib.pyplot as plt
plt.rcParams['figure.constrained_layout.use'] = True

with open("./results/mcm_ids.csv", "r") as file:
    data = csv.reader(file.readlines())
    header = data.__next__()
    mcm_history: list[list[int]] = [[] for i in range(len(header))]
    for line in data:
        for i, mcm in zip(range(len(header)), line):
            # mcm_history[int(mcm)].append(i)
            mcm_history[int(mcm)].append(-i)

fig, axes = plt.subplots(len(header), 1)
for mcm, data in zip(range(len(header)), mcm_history):
    if not mcm == len(header) - 1:
        axes[mcm].set(xticks=[])
    axes[mcm].plot(range(len(data)), data, "-", label=f"{mcm}")
plt.show()