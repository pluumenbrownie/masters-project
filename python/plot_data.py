import csv
import matplotlib.pyplot as plt
import numpy as np
import typer
from time import localtime

plt.rcParams["figure.constrained_layout.use"] = True

app = typer.Typer()


@app.command()
def par_temp_pools(save: bool = False) -> None:
    with open("./results/mcm_ids.csv", "r") as file:
        data = csv.reader(file.readlines())
        header = data.__next__()
        mcm_history: list[list[int]] = [[] for i in range(len(header))]
        for line in data:
            for i, mcm_str in zip(range(len(header)), line):
                # mcm_history[int(mcm)].append(i)
                mcm_history[int(mcm_str)].append(-i)

    fig, axes = plt.subplots(len(header), 1)
    for mcm, pool_data in zip(range(len(header)), mcm_history):
        if not mcm == len(header) - 1:
            axes[mcm].set(xticks=[])
        axes[mcm].plot(range(len(pool_data)), pool_data, "-", label=f"{mcm}")

    if save:
        plt.savefig(f"./results/par_temp_pools_{localtime()}.svg")
    plt.show()


@app.command()
def mutation_history(save: bool = False) -> None:
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
                mutation_accepted_history.append(
                    {}.fromkeys(["Swap", "Split", "Merge"], 0)
                )

            mutation_history[-1][line[0]] += 1
            if line[1] == "true":
                mutation_accepted_history[-1][line[0]] += 1

            counter += 1
            if counter == chunk_size:
                counter = 0

    plot_data: dict[str, list[int]] = {
        "Merge": [],
        "Swap": [],
        "Split": [],
    }
    plot_accepted_data: dict[str, list[int]] = {
        "Merge": [],
        "Swap": [],
        "Split": [],
    }
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

    if save:
        plt.savefig(f"./results/mutation_history_{localtime()}.svg")
    plt.show()


@app.command()
def mutation_pool_history(save: bool = False) -> None:
    chunk_size = 200

    data_per_temperature: dict[float, list[list[str]]] = {}
    with open("./results/mutation_history.csv", "r") as file:
        data = csv.reader(file.readlines())
        header = data.__next__()
        for line in data:
            data_per_temperature.setdefault(float(line[2]), []).append(line[:2])

    plot_data_temperature: dict[float, dict[str, list[int]]] = {}
    plot_accepted_data_temperature: dict[float, dict[str, list[int]]] = {}
    for temp, data_for_temp in sorted(
        data_per_temperature.items(), key=lambda i: float(i[0]), reverse=True
    ):
        mutation_history: list[dict[str, int]] = []
        mutation_accepted_history: list[dict[str, int]] = []
        counter = 0
        for line in data_for_temp:
            if not counter:
                mutation_history.append({}.fromkeys(["Swap", "Split", "Merge"], 0))
                mutation_accepted_history.append(
                    {}.fromkeys(["Swap", "Split", "Merge"], 0)
                )

            mutation_history[-1][line[0]] += 1
            if line[1] == "true":
                mutation_accepted_history[-1][line[0]] += 1

            counter += 1
            if counter == chunk_size:
                counter = 0

        plot_data: dict[str, list[int]] = {
            "Merge": [],
            "Swap": [],
            "Split": [],
        }
        plot_accepted_data: dict[str, list[int]] = {
            "Merge": [],
            "Swap": [],
            "Split": [],
        }
        btm = []

        for chunk in mutation_history:
            for key, value in chunk.items():
                plot_data[key].append(value)
            btm.append(0)
        for chunk in mutation_accepted_history:
            for key, value in chunk.items():
                plot_accepted_data[key].append(value)

        plot_data_temperature[temp] = plot_data
        plot_accepted_data_temperature[temp] = plot_accepted_data

    bottom = np.array(btm)

    plot_amount = len(plot_accepted_data_temperature)
    fig, axes = plt.subplots(plot_amount, 2)
    for n, ((temp, plot_data), plot_accepted_data) in enumerate(
        zip(plot_data_temperature.items(), plot_accepted_data_temperature.values())
    ):
        for merge_type, counts in plot_data.items():
            axes[n, 0].bar(
                range(len(counts)), counts, 1.0, label=merge_type, bottom=bottom
            )
            bottom += np.array(counts)
        axes[n, 0].set(title=f"{temp = }")
        if not n == plot_amount - 1:
            axes[n, 0].set(xticks=[])
        else:
            axes[n, 0].legend()

        bottom[:] = 0

        for merge_type, counts in plot_accepted_data.items():
            axes[n, 1].bar(
                range(len(counts)), counts, 1.0, label=merge_type, bottom=bottom
            )
            bottom += np.array(counts)
        if not n == plot_amount - 1:
            axes[n, 1].set(xticks=[])
        bottom[:] = 0

    if save:
        plt.savefig(f"./results/mutation_pool_history_{localtime()}.svg")
    plt.show()


if __name__ == "__main__":
    app()
