import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator
import matplotlib as mpl
import numpy as np
import typer
from time import asctime

plt.rcParams["figure.constrained_layout.use"] = True

app = typer.Typer()


@app.command()
def par_temp_pools(save: bool = False) -> None:
    header, mcm_history = read_par_temp_ids()
    log_e = read_par_temp_log_e(header)

    fig, axes = plt.subplots(len(header), 1, figsize=(10, 10))
    for mcm, pool_data, mcm_log_e in zip(range(len(header)), mcm_history, log_e):
        axes[mcm].set_yticks(
            [0, -len(header) + 1],
            [round(float(header[0]), ndigits=3), round(float(header[-1]), ndigits=3)],
        )
        axes[mcm].yaxis.set_minor_locator(MultipleLocator(1))
        axes[mcm].grid(axis="y", which="both")
        axes[mcm].set(ylim=(-len(header) + 0.7, 0.3), xlim=(0, len(pool_data)))
        axes[mcm].plot(range(len(pool_data)), pool_data, "b-", label="pool")
        axes[mcm].set_ylabel("T")

        ax2 = axes[mcm].twinx()
        ax2.plot(range(len(mcm_log_e)), mcm_log_e, "r-", label="log E")
        ax2.set_ylabel("log(E)")

        if not mcm == len(header) - 1:
            axes[mcm].set(xticks=[])
        else:
            axes[mcm].plot([0, 0], [10, 10], "r-", label="log E")
            axes[mcm].legend()

    if save:
        plt.savefig(f"./results/par_temp_pools_{asctime()}.svg")
    plt.show()


def read_par_temp_log_e(header) -> list[list[float]]:
    with open("./results/par_temp_log_e.csv", "r") as file:
        data = csv.reader(file.readlines())
        _ = next(data)
        log_e: list[list[float]] = [[] for i in header]
        for line in data:
            log_e[int(line[0])].append(float(line[1]))
    return log_e


def read_par_temp_ids() -> tuple[list[str], list[list[int]]]:
    with open("./results/mcm_ids.csv", "r") as file:
        data = csv.reader(file.readlines())
        header = next(data)
        mcm_history: list[list[int]] = [[] for i in range(len(header))]
        for line in data:
            for i, mcm_str in zip(range(len(header)), line):
                mcm_history[int(mcm_str)].append(-i)
    return header, mcm_history


@app.command()
def par_temp_log_e(save: bool = False) -> None:
    header, _ = read_par_temp_ids()
    log_e = read_par_temp_log_e(header)

    for mcm, mcm_log_e in enumerate(log_e):
        plt.plot(range(len(mcm_log_e)), mcm_log_e, "-", label=f"{mcm}")

    if save:
        plt.savefig(f"./results/par_temp_log_e_{asctime()}.svg")
    plt.legend()
    plt.show()


@app.command()
def mutation_history(save: bool = False) -> None:
    chunk_size = 100
    with open("./results/mutation_history.csv", "r") as file:
        data = csv.reader(file.readlines())
        header = next(data)

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

    fig, axes = plt.subplots(2, 1, figsize=(10, 10))
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
        plt.savefig(f"./results/mutation_history_{asctime()}.svg")
    plt.show()


@app.command()
def par_temp_mutations(save: bool = False) -> None:
    chunk_size = 200

    data_per_temperature: dict[float, list[list[str]]] = {}
    with open("./results/mutation_history.csv", "r") as file:
        data = csv.reader(file.readlines())
        header = next(data)
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
    fig, axes = plt.subplots(plot_amount, 2, figsize=(10, 10))
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
        plt.savefig(f"./results/mutation_pool_history_{(asctime())}.svg")
    plt.show()


@app.command()
def log_e(save: bool = False) -> None:
    with open("./results/log_e_history.csv") as file:
        data = csv.reader(file.readlines())
        header = next(data)

    log_e_data = [float(i[0]) for i in data]
    plt.plot(log_e_data)

    if save:
        plt.savefig(f"./results/log_e_history_{(asctime())}.svg")
    plt.show()


@app.command()
def log_e_greedy(save: bool = False) -> None:
    with open("./results/log_e_greedy.csv") as file:
        data = csv.reader(file.readlines())
        header = next(data)

    data_list = [i for i in data]

    algorithm = None
    log_e_data = []
    steps = []
    for nr, line in enumerate(data_list):
        if not algorithm == line[0]:
            if algorithm:
                plt.plot(steps, log_e_data, label=algorithm)
                log_e_data.clear()
                steps.clear()
            algorithm = line[0]
        log_e_data.append(float(line[2]))
        steps.append(nr)
    plt.plot(steps, log_e_data, label=algorithm)
    plt.grid(True)
    plt.legend()

    if save:
        plt.savefig(f"./results/log_e_greedy_history_{(asctime())}.svg")
    plt.show()


@app.command()
def log_e_annealing(save: bool = False) -> None:
    with open("./results/log_e_annealing.csv") as file:
        data = csv.reader(file.readlines())
        header = next(data)

    log_e_data = []
    temperature_data = []
    for line in data:
        log_e_data.append(float(line[0]))
        temperature_data.append(float(line[1]))

    fig, ax = plt.subplots(1, 1)
    ax2 = ax.twinx()
    ax.plot(range(len(log_e_data)), log_e_data, label="log(E)")
    ax2.semilogy(range(len(log_e_data)), temperature_data, "r-", label="T")
    fig.legend()

    if save:
        plt.savefig(f"./results/log_e_annealing_history_{(asctime())}.svg")
    plt.show()


@app.command()
def log_e_hill_climbing(save: bool = False) -> None:
    with open("./results/log_e_hill_climbing.csv") as file:
        data = csv.reader(file.readlines())
        header = next(data)

    log_e_data = []
    history_data = []
    for line in data:
        log_e_data.append(float(line[0]))
        history_data.append(float(line[1]))

    fig, ax = plt.subplots(1, 1)
    ax2 = ax.twinx()
    ax.plot(range(len(log_e_data)), log_e_data, label="log(E)")
    ax2.plot(range(len(log_e_data)), history_data, "r-", label="History log(E)")
    fig.legend()

    if save:
        plt.savefig(f"./results/log_e_hill_climbing_{(asctime())}.svg")
    plt.show()


@app.command()
def par_temp(save: bool = False) -> None:
    par_temp_log_e(save)
    par_temp_pools(save)
    par_temp_mutations(save)


@app.command()
def annealing(save: bool = False) -> None:
    log_e_annealing(save)
    mutation_history(save)


@app.command()
def hill_climbing(save: bool = False) -> None:
    log_e_hill_climbing(save)
    mutation_history(save)


if __name__ == "__main__":
    app()
