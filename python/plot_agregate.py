import csv
from typing import Any
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator
import matplotlib as mpl
import numpy as np
import typer
from time import asctime
import json
import os
import re
from scipy.stats import linregress  # type: ignore
from itertools import cycle

plt.rcParams["figure.constrained_layout.use"] = True

app = typer.Typer()


def find_floats(text: str) -> list[float]:
    pattern = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
    return [float(match) for match in re.findall(pattern, text)]


@app.command()
def prove_serial_time_dependence(save: bool = False):
    prove_time_dependence(save, "./results/simulated_annealing")


def prove_time_dependence(save: bool, folder: str) -> None:
    time_seconds = []
    iccs_covered = []

    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        time_seconds.append(float(content["time_seconds"]))
        iccs_covered.append(int(content["iccs_checked"]))

    x_range = [min(iccs_covered) * 0.995, max(iccs_covered) * 1.005]
    fit = linregress(iccs_covered, time_seconds)
    fitline = lambda x: fit.slope * x + fit.intercept

    plt.plot(iccs_covered, time_seconds, "o", label="datapoints")
    plt.plot(x_range, list(map(fitline, x_range)), "--", label="linear fit")

    plt.xlabel("ICCs calculated")
    plt.ylabel("time (s)")
    plt.title("Runtime over ICCs calculated")
    plt.legend()
    plt.xlim(*x_range)
    plt.show()


@app.command()
def plot_sa(save: bool = True) -> None:
    folder = "./results/simulated_annealing"
    temp_plotted = False
    fig, ax = plt.subplots(1, 1)
    ax2 = ax.twinx()

    for filename in filter(lambda name: name[-4:] == ".csv", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = csv.reader(file.readlines())
        log_e = []
        temperature = []

        next(content)
        for line in content:
            log_e.append(float(line[0]))
            temperature.append(float(line[1]))

        ax.plot(range(len(log_e)), log_e, "-", color="#1f77b4", alpha=0.3)
        if not temp_plotted:
            ax2.semilogy(range(len(temperature)), temperature, "r-", label="T")
            temp_plotted = True

    ax.set_xlabel("step")
    ax.set_ylabel("log(E)")
    ax2.set_ylabel("log(T)")
    fig.suptitle("Log(E) of Simulated Annealing")
    plt.legend()
    ax.set_xlim(-50, len(log_e))
    plt.show()


@app.command()
def prove_parallel_time_dependence(save: bool = False) -> None:
    prove_time_dependence(save, "./results/parallel_tempering_timing")


@app.command()
def sa_shapes(save: bool = True) -> None:
    folder = "./results/simulated_annealing_shapes"
    data: dict[str, dict[str, list[float | int]]] = {}

    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        name = filename[7:-14]
        data.setdefault(name, {}).setdefault("log_e", []).append(
            float(content["log_e"])
        )
        data.setdefault(name, {}).setdefault("iccs", []).append(
            int(content["iccs_checked"])
        )

    sorted_data = [x for x in data.items()]
    sorted_data = sorted(sorted_data, key=lambda x: x[0])
    colors = cycle(plt.cm.tab10.colors)  # type: ignore

    def shapes(name):
        if "exp" in name:
            return "^"
        if "lin" in name:
            return "s"
        if "cus" in name:
            return "x"
        return "o"

    best: None | tuple[str, float] = None
    for (name, points), color in zip(sorted_data, colors):
        log_e = np.array(points["log_e"])
        iccs = np.array(points["iccs"])
        plt.plot(iccs, log_e, shapes(name), label=str(name), color=color, alpha=0.75)

        best_log_e = log_e.max()
        if not best or best_log_e > best[1]:
            best = (name, best_log_e)

        log_e_percentiles = np.percentile(log_e, [25, 50, 75])
        iccs_percentiles = np.percentile(iccs, [25, 50, 75])
        plt.errorbar(
            iccs_percentiles[1],
            log_e_percentiles[1],
            xerr=[
                [iccs_percentiles[1] - iccs_percentiles[0]],
                [iccs_percentiles[2] - iccs_percentiles[1]],
            ],
            yerr=[
                [log_e_percentiles[1] - log_e_percentiles[0]],
                [log_e_percentiles[2] - log_e_percentiles[1]],
            ],
            fmt=shapes(name),
            color=color,
            ms=9,
            mec="k",
            capsize=4,
            lw=2,
        )

    plt.legend()
    plt.show()

    print(f"best: {best[0]} with {best[1]}")  # type: ignore

    # plt.boxplot([d[1] for d in log_e_data], tick_labels=[str(d[0]) for d in log_e_data])
    # plt.show()


@app.command()
def sa_exp_decay(save: bool = True) -> None:
    folder = "./results/simulated_annealing_longer"
    data: dict[float, dict[str, list[float | int]]] = {}

    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        decay = find_floats(filename)[1]
        data.setdefault(decay, {}).setdefault("log_e", []).append(
            float(content["log_e"])
        )
        data.setdefault(decay, {}).setdefault("iccs", []).append(
            int(content["iccs_checked"])
        )

    sorted_data = [x for x in data.items()]
    sorted_data = sorted(sorted_data, key=lambda x: x[0])
    colors = plt.cm.tab10.colors  # type: ignore

    best: None | tuple[str, float] = None
    for (name, points), color in zip(sorted_data, colors):
        log_e = np.array(points["log_e"])
        iccs = np.array(points["iccs"])
        plt.plot(iccs, log_e, "o", label=str(name), color=color, alpha=0.75)

        best_log_e = log_e.max()
        if not best or best_log_e > best[1]:
            best = (str(name), best_log_e)

        log_e_percentiles = np.percentile(log_e, [25, 50, 75])
        iccs_percentiles = np.percentile(iccs, [25, 50, 75])
        plt.errorbar(
            iccs_percentiles[1],
            log_e_percentiles[1],
            xerr=[
                [iccs_percentiles[1] - iccs_percentiles[0]],
                [iccs_percentiles[2] - iccs_percentiles[1]],
            ],
            yerr=[
                [log_e_percentiles[1] - log_e_percentiles[0]],
                [log_e_percentiles[2] - log_e_percentiles[1]],
            ],
            fmt="o",
            color=color,
            ms=9,
            mec="k",
            capsize=4,
            lw=2,
        )

    plt.legend()
    plt.show()

    print(f"best: {best[0]} with {best[1]}")  # type: ignore

    # plt.boxplot([d[1] for d in log_e_data], tick_labels=[str(d[0]) for d in log_e_data])
    # plt.show()


@app.command()
def lahc_hist_size(save: bool = True) -> None:
    folder = "./results/lahc_hist_size"
    data: dict[float, dict[str, list[float | int]]] = {}

    best: None | tuple[str, float] = None
    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        history_size = int(find_floats(filename)[1])

        log_e = content["log_e"]
        if not best or log_e > best[1]:
            best = (filename, log_e)

        data.setdefault(history_size, {}).setdefault("log_e", []).append(float(log_e))
        data.setdefault(history_size, {}).setdefault("iccs", []).append(
            int(content["iccs_checked"])
        )

    sorted_data = [x for x in data.items()]
    sorted_data = sorted(sorted_data, key=lambda x: x[0])
    colors = plt.cm.tab10.colors  # type: ignore

    for (name, points), color in zip(sorted_data, colors):
        log_e = np.array(points["log_e"])
        iccs = np.array(points["iccs"])
        plt.plot(iccs, log_e, "o", label=str(name), color=color, alpha=0.75)

        log_e_percentiles = np.percentile(log_e, [25, 50, 75])
        iccs_percentiles = np.percentile(iccs, [25, 50, 75])
        plt.errorbar(
            iccs_percentiles[1],
            log_e_percentiles[1],
            xerr=[
                [iccs_percentiles[1] - iccs_percentiles[0]],
                [iccs_percentiles[2] - iccs_percentiles[1]],
            ],
            yerr=[
                [log_e_percentiles[1] - log_e_percentiles[0]],
                [log_e_percentiles[2] - log_e_percentiles[1]],
            ],
            fmt="o",
            color=color,
            ms=9,
            mec="k",
            capsize=4,
            lw=2,
        )

    plt.legend()
    plt.show()

    for (name, points), color in zip(sorted_data, colors):
        log_e_ranked = np.array(points["log_e"])

        log_e_ranked.sort()
        probability = np.arange(0.0, 1.0, step=1 / len(log_e_ranked))
        # probability = np.arange(1 / len(log_e_ranked), 1.01, step=1 / len(log_e_ranked))
        plt.plot(-log_e_ranked, 1 - probability, "o-", color=color, label=str(name))
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.show()

    print(f"best: {best[0]} with {best[1]}")  # type: ignore

    # plt.boxplot([d[1] for d in log_e_data], tick_labels=[str(d[0]) for d in log_e_data])
    # plt.show()


class MapToShape:
    def __init__(self, keys: list[str], values: list[str], names: list[str]) -> None:
        self.keys = keys
        self.values = values
        self.names = names

    def __call__(self, name: str) -> str:
        for key, value in zip(self.keys, self.values):
            if key in name:
                return value
        raise ValueError

    def plot_legend(self):
        for value, name in zip(self.values, self.names):
            plt.plot([], [], value, color="#00000000", ms=9, mec="k", lw=2, label=name)

    def plot_legend_custom(self, ax):
        for value, name in zip(self.values, self.names):
            ax.plot([], [], value, color="#00000000", ms=9, mec="k", lw=2, label=name)


class MapToColor:
    def __init__(
        self,
        keys: list[str],
        values: list[Any],
        names: list[str],
    ) -> None:
        self.keys = keys
        self.values = values
        self.names = names

    def __call__(self, name: str) -> Any:
        for key, value in zip(self.keys, self.values):
            if key in name:
                return value
        raise ValueError

    def plot_legend(self):
        for value, name in zip(self.values, self.names):
            plt.plot([], [], "s", color=value, label=name)

    def plot_legend_custom(self, ax):
        for value, name in zip(self.values, self.names):
            ax.plot([], [], "s", color=value, label=name)


@app.command()
def greedy(save: bool = True) -> None:
    folder = "./results/greedy"
    data: dict[str, dict[str, list[float | int]]] = {}

    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        name = filename[15:-5]
        data.setdefault(name, {}).setdefault("log_e", []).append(
            float(content["log_e"])
        )
        data.setdefault(name, {}).setdefault("iccs", []).append(
            int(content["iccs_checked"])
        )

    sorted_data = [x for x in data.items()]
    sorted_data = sorted(sorted_data, key=lambda x: x[0])
    colors = MapToColor(
        ["n", "m", "c"], plt.cm.tab10.colors[0:3], ["None", "Merge", "Constructive"]  # type: ignore
    )

    shapes = MapToShape(["t", "l", "nr"], ["^", "s", "o"], ["Tabu", "Local", "None"])

    for name, points in sorted_data:
        log_e = np.array(points["log_e"])
        iccs = np.array(points["iccs"])
        plt.plot(iccs, log_e, shapes(name), color=colors(name), ms=9, mec="k", lw=2)

    colors.plot_legend()
    shapes.plot_legend()

    plt.legend()
    plt.show()

    fig, (merge, constructive, none) = plt.subplots(3, 1)

    best: None | tuple[str, float] = None
    for name, points in sorted_data:
        log_e = np.array(points["log_e"])
        iccs = np.array(points["iccs"])

        best_log_e = log_e.max()
        if not best or best_log_e > best[1]:
            best = (name, best_log_e)

        if "t" in name:
            if "m" in name:
                merge.plot(iccs, log_e, shapes(name), ms=9, mec="k", lw=2, label=name)
            elif "c" in name:
                constructive.plot(
                    iccs, log_e, shapes(name), ms=9, mec="k", lw=2, label=name
                )
            else:
                none.plot(iccs, log_e, shapes(name), ms=9, mec="k", lw=2, label=name)
    merge.legend()
    constructive.legend()
    none.legend()
    plt.show()

    print(f"best: {best[0]} with {best[1]}")  # type: ignore


@app.command()
def adaptive(save: bool = True) -> None:
    folder = "./results/adaptive_annealing"
    data: dict[str, dict[str, list[float | int]]] = {}

    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        name = filename[7:-14]
        data.setdefault(name, {}).setdefault("log_e", []).append(
            float(content["log_e"])
        )
        data.setdefault(name, {}).setdefault("iccs", []).append(
            int(content["iccs_checked"])
        )

    sorted_data = [x for x in data.items()]
    sorted_data = sorted(sorted_data, key=lambda x: x[0])

    colors = MapToColor(
        ["0.004", "0.003", "0.002", "0.0015", "0.001", "0.0007", "0.0005"],
        plt.cm.tab10.colors[0:7],  # type: ignore
        ["0.004", "0.003", "0.002", "0.0015", "0.001", "0.0007", "0.0005"],
    )

    shapes = MapToShape(
        ["0.95", "0.9", "0.85", "0.8", "0.75"],
        ["o", "s", "^", "*", "P"],
        ["0.95", "0.90", "0.85", "0.80", "0.75"],
    )

    best: None | tuple[str, float] = None
    for name, points in sorted_data:
        log_e = np.array(points["log_e"])
        iccs = np.array(points["iccs"])
        plt.plot(iccs, log_e, shapes(name), color=colors(name), alpha=0.75)

        best_log_e = log_e.max()
        if not best or best_log_e > best[1]:
            best = (name, best_log_e)

        log_e_percentiles = np.percentile(log_e, [25, 50, 75])
        iccs_percentiles = np.percentile(iccs, [25, 50, 75])
        plt.errorbar(
            iccs_percentiles[1],
            log_e_percentiles[1],
            xerr=[
                [iccs_percentiles[1] - iccs_percentiles[0]],
                [iccs_percentiles[2] - iccs_percentiles[1]],
            ],
            yerr=[
                [log_e_percentiles[1] - log_e_percentiles[0]],
                [log_e_percentiles[2] - log_e_percentiles[1]],
            ],
            fmt=shapes(name),
            color=colors(name),
            ms=9,
            mec="k",
            capsize=4,
            lw=2,
        )
    colors.plot_legend()
    shapes.plot_legend()

    plt.legend()
    plt.show()

    print(f"best: {best[0]} with {best[1]}")  # type: ignore


if __name__ == "__main__":
    app()
