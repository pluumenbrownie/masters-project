import csv
from typing import Any
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import typer
from time import asctime
import json
import os
import re
from scipy.stats import linregress  # type: ignore
from itertools import cycle, chain

plt.rcParams["figure.constrained_layout.use"] = True
plt.rcParams["figure.figsize"] = (10, 8)

app = typer.Typer()


def find_floats(text: str) -> list[float]:
    pattern = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
    return [float(match) for match in re.findall(pattern, text)]


def sort_datapoints(sorted_data):
    return sorted(sorted_data, key=lambda x: find_floats(x[0])[0])  # type: ignore


@app.command()
def prove_annealing_time_dependence(save: bool = False):
    prove_time_dependence(save, "./results/simulated_annealing", "Annealing")


def prove_time_dependence(save: bool, folder: str, name: str) -> None:
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
    plt.title(f"{name} time dependence")
    plt.legend()
    plt.xlim(*x_range)
    if save:
        plt.savefig(f"./results/plots/{name}_time_dependence_{asctime()}.svg")
    plt.show()


@app.command()
def plot_sa(save: bool = False) -> None:
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
    if save:
        plt.savefig(f"./results/plots/simulated_annealing_{asctime()}.svg")
    plt.show()


@app.command()
def plot_sa_adaptive(save: bool = False) -> None:
    folder = "./results/adaptive_annealing"
    temp_plotted = False
    fig, ax = plt.subplots(1, 1)
    ax2 = ax.twinx()

    for filename in filter(
        lambda name: name[-4:] == ".csv" and "0.0015_0.75_bU6" in name,
        os.listdir(folder),
    ):
        with open(folder + "/" + filename, "r") as file:
            content = csv.reader(file.readlines())
        log_e = []
        temperature = []

        next(content)
        for line in content:
            log_e.append(float(line[0]))
            temperature.append(float(line[1]))

        ax.plot(range(len(log_e)), log_e, "-", color="#1f77b4")
        ax2.semilogy(range(len(temperature)), temperature, "r-", label="T")

    ax.set_xlabel("step")
    ax.set_ylabel("log(E)")
    ax2.set_ylabel("log(T)")
    fig.suptitle("Log(E) of Adaptive Annealing")
    plt.legend()
    ax.set_xlim(-50, len(log_e))
    if save:
        plt.savefig(f"./results/plots/adaptive_annealing_{asctime()}.svg")
    plt.show()


@app.command()
def prove_parallel_time_dependence(save: bool = False) -> None:
    prove_time_dependence(
        save, "./results/parallel_tempering_timing", "Parallel tempering"
    )


@app.command()
def sa_shapes(save: bool = False) -> None:
    folder = "./results/simulated_annealing_shapes"
    data: dict[str, dict[str, list[float | int]]] = {}

    best: None | tuple[str, float] = None
    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        name = filename[7:-14]
        best_log_e = float(content["log_e"])
        data.setdefault(name, {}).setdefault("log_e", []).append(best_log_e)
        data.setdefault(name, {}).setdefault("iccs", []).append(
            int(content["iccs_checked"])
        )
        if not best or best_log_e > best[1]:
            best = (filename[7:-5], best_log_e)

    print(f"best: {best[0]} with {best[1]}")  # type: ignore

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

    for (name, points), color in zip(sorted_data, colors):
        log_e = np.array(points["log_e"])
        iccs = np.array(points["iccs"])
        plt.plot(iccs, log_e, shapes(name), label=str(name), color=color, alpha=0.75)

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

    plt.legend(ncol=2)
    plt.xlabel("ICCs calculated")
    plt.ylabel("log(E)")
    plt.title("Annealing shapes comparison")
    if save:
        plt.savefig(f"./results/plots/annealing_shapes_{asctime()}.svg")
    plt.show()


@app.command()
def sa_exp_decay(save: bool = False) -> None:
    folder = "./results/simulated_annealing_longer"
    data: dict[float, dict[str, list[float | int]]] = {}

    best: None | tuple[str, float] = None
    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        decay = find_floats(filename)[1]
        best_log_e = float(content["log_e"])
        data.setdefault(decay, {}).setdefault("log_e", []).append(best_log_e)
        data.setdefault(decay, {}).setdefault("iccs", []).append(
            int(content["iccs_checked"])
        )
        if not best or best_log_e > best[1]:
            best = (filename[7:-5], best_log_e)

    print(f"best: {best[0]} with {best[1]}")  # type: ignore

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
    plt.xlabel("ICCs calculated")
    plt.ylabel("log(E)")
    plt.title("Exponential decay speed")
    if save:
        plt.savefig(f"./results/plots/annealing_exp_decay_{asctime()}.svg")
    plt.show()


@app.command()
def lahc_hist_size(save: bool = False) -> None:
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
    print(f"best: {best[0]} with {best[1]}")  # type: ignore

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
    plt.xlabel("ICCs calculated")
    plt.ylabel("log(E)")
    plt.title("Late Acceptance Hill Climbing")
    if save:
        plt.savefig(f"./results/plots/lahc_hist_size_{asctime()}.svg")
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
    if save:
        plt.savefig(f"./results/plots/lahc_loglog_{asctime()}.svg")
    plt.show()


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
def greedy(save: bool = False) -> None:
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
    sorted_data = sort_datapoints(sorted_data)
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
    plt.xlabel("ICCs calculated")
    plt.ylabel("log(E)")
    plt.title("Greedy Strategy Comparison")
    if save:
        plt.savefig(f"./results/plots/greedy_{asctime()}.svg")
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
    merge.legend(ncols=3)
    constructive.legend(ncols=3)
    none.legend(ncols=3)
    if save:
        plt.savefig(f"./results/plots/greedy_split_{asctime()}.svg")
    plt.show()

    print(f"best: {best[0]} with {best[1]}")  # type: ignore


@app.command()
def tabu(save: bool = False) -> None:
    folder = "./results/tabu_divergence"
    data: dict[str, dict[str, list[float | int]]] = {}

    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        name = filename[7:-5]
        data.setdefault(name, {}).setdefault("log_e", []).append(
            float(content["log_e"])
        )
        data.setdefault(name, {}).setdefault("iccs", []).append(
            int(content["iccs_checked"])
        )

    sorted_data = [x for x in data.items()]
    sorted_data = sort_datapoints(sorted_data)

    best: None | tuple[str, float] = None
    for name, points in sorted_data:

        log_e = np.array(points["log_e"])
        iccs = np.array(points["iccs"])

        best_log_e = log_e.max()
        if not best or best_log_e > best[1]:
            best = (name, best_log_e)

        plt.plot(iccs, log_e, "o", ms=9, mec="k", lw=2, label=name)
    plt.legend(ncols=2)
    plt.xlabel("ICCs calculated")
    plt.ylabel("log(E)")
    plt.title("Tabu list size")
    if save:
        plt.savefig(f"./results/plots/tabu_{asctime()}.svg")
    plt.show()

    colors = chain(plt.cm.tab20c.colors, plt.cm.tab20b.colors)  # type: ignore

    for filename in list_csvs(folder):
        with open(folder + "/" + filename, "r") as file:
            content = csv.reader(file.readlines())
        next(content)

        name = filename[6:-4]

        log_e_list = []
        for solver, size, log_e in content:
            if not solver == "Constructive":
                log_e_list.append(float(log_e))
        plt.plot(
            range(len(log_e_list)), log_e_list, "-", color=next(colors), label=name
        )
    plt.legend(ncols=4)
    plt.xlabel("ICCs calculated")
    plt.ylabel("log(E)")
    plt.title("Tabu list size trajectories")
    plt.ylim((-2430000, -2380000))
    if save:
        plt.savefig(f"./results/plots/tabu_trajectory_{asctime()}.svg")
    plt.show()

    print(f"best: {best[0]} with {best[1]}")  # type: ignore


@app.command()
def tabu_2(save: bool = False) -> None:
    folder = "./results/tabu_divergence_lower"
    data: dict[str, dict[str, list[float | int]]] = {}

    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        name = filename[7:-5]
        data.setdefault(name, {}).setdefault("log_e", []).append(
            float(content["log_e"])
        )
        data.setdefault(name, {}).setdefault("iccs", []).append(
            int(content["iccs_checked"])
        )

    sorted_data = [x for x in data.items()]
    sorted_data = sort_datapoints(sorted_data)

    best: None | tuple[str, float] = None
    for name, points in sorted_data:

        log_e = np.array(points["log_e"])
        iccs = np.array(points["iccs"])

        best_log_e = log_e.max()
        if not best or best_log_e > best[1]:
            best = (name, best_log_e)

        plt.plot(iccs, log_e, "o", ms=9, mec="k", lw=2, label=name)
    plt.legend(ncols=2)
    plt.xlabel("ICCs calculated")
    plt.ylabel("log(E)")
    plt.title("Exponential decay speed")
    if save:
        plt.savefig(f"./results/plots/tabu2_{asctime()}.svg")
    plt.show()

    colors = chain(plt.cm.tab20c.colors, plt.cm.tab20b.colors)  # type: ignore

    for filename in list_csvs(folder):
        with open(folder + "/" + filename, "r") as file:
            content = csv.reader(file.readlines())
        next(content)

        name = filename[6:-4]

        log_e_list = []
        for solver, size, log_e in content:
            if not solver == "Constructive":
                log_e_list.append(float(log_e))
        plt.plot(
            range(len(log_e_list)), log_e_list, "-", color=next(colors), label=name
        )
    plt.legend(ncols=4)
    plt.xlabel("ICCs calculated")
    plt.ylabel("log(E)")
    plt.title("Tabu list size trajectories")
    plt.ylim((-2430000, -2380000))
    if save:
        plt.savefig(f"./results/plots/tabu2_trajectory_{asctime()}.svg")
    plt.show()

    print(f"best: {best[0]} with {best[1]}")  # type: ignore


def list_csvs(folder):
    return sorted(
        list(filter(lambda name: name[-4:] == ".csv", os.listdir(folder))),
        key=lambda f: int(f[8:-4]),
    )  # type: ignore


@app.command()
def tabu_only(save: bool = False) -> None:
    folder = "./results/pure_tabu"
    data: dict[str, dict[str, list[float | int]]] = {}

    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        name = filename[7:-5]
        data.setdefault(name, {}).setdefault("log_e", []).append(
            float(content["log_e"])
        )
        data.setdefault(name, {}).setdefault("iccs", []).append(
            int(content["iccs_checked"])
        )

    sorted_data = [x for x in data.items()]
    sorted_data = sort_datapoints(sorted_data)
    cmap = mpl.colormaps["gist_rainbow"]
    colors = cmap(np.linspace(0, 1, len(sorted_data)))

    best: None | tuple[str, float] = None
    for (name, points), color in zip(sorted_data, colors):

        log_e = np.array(points["log_e"])
        iccs = np.array(points["iccs"])

        best_log_e = log_e.max()
        if not best or best_log_e > best[1]:
            best = (name, best_log_e)

        plt.plot(iccs, log_e, "o", color=color, ms=9, mec="k", lw=2)

    plt.xlabel("ICCs calculated")
    plt.ylabel("log(E)")
    plt.title("Tabu only strategy")
    if save:
        plt.savefig(f"./results/plots/tabu_only_{asctime()}.svg")
    plt.show()
    print(f"best: {best[0]} with {best[1]}")  # type: ignore

    for filename, color in zip(list_csvs(folder), colors):
        with open(folder + "/" + filename, "r") as file:
            content = csv.reader(file.readlines())
        next(content)

        name = filename[6:-4]

        log_e_list = []
        for solver, size, log_e in content:
            log_e_list.append(float(log_e))

        if name == best[0]:  # type: ignore
            plt.plot(
                range(len(log_e_list)),
                log_e_list,
                "-",
                color="#000000",
                label=f"best: {name}",
            )
        else:
            plt.plot(range(len(log_e_list)), log_e_list, "-", color=color)

    plt.legend()
    plt.xlabel("ICCs calculated")
    plt.ylabel("log(E)")
    plt.title("Tabu only trajectories")
    plt.ylim((-2450000, -2380000))
    if save:
        plt.savefig(f"./results/plots/tabu_only_trajectory_{asctime()}.svg")
    plt.show()


@app.command()
def merge_beam(save: bool = False) -> None:
    folder = "./results/merge_beam"
    data: dict[str, dict[str, list[float | int]]] = {}

    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        name = filename[7:-5]
        data.setdefault(name, {}).setdefault("log_e", []).append(
            float(content["log_e"])
        )
        data.setdefault(name, {}).setdefault("iccs", []).append(
            int(content["iccs_checked"])
        )

    sorted_data = [x for x in data.items()]
    sorted_data = sort_datapoints(sorted_data)
    cmap = mpl.colormaps["gist_rainbow"]
    colors = cmap(np.linspace(0, 1, len(sorted_data)))

    best: None | tuple[str, float] = None
    for (name, points), color in zip(sorted_data, colors):

        log_e = np.array(points["log_e"])
        iccs = np.array(points["iccs"])

        best_log_e = log_e.max()
        if not best or best_log_e > best[1]:
            best = (name, best_log_e)

        plt.plot(iccs, log_e, "o", color=color, ms=9, mec="k", lw=2, label=name)

    plt.xlabel("ICCs calculated")
    plt.ylabel("log(E)")
    plt.title("Merge beam size")
    if save:
        plt.savefig(f"./results/plots/merge_beam_{asctime()}.svg")
    plt.show()
    print(f"best: {best[0]} with {best[1]}")  # type: ignore


@app.command()
def prove_greedy_time_dependence(save: bool = False):
    prove_time_dependence(save, "./results/construct_beam", "Greedy")


@app.command()
def construct_beam(save: bool = False) -> None:
    folder = "./results/construct_beam"
    data: dict[str, dict[str, list[float | int]]] = {}

    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        name = filename[7:-5]
        data.setdefault(name, {}).setdefault("log_e", []).append(
            float(content["log_e"])
        )
        data.setdefault(name, {}).setdefault("iccs", []).append(
            int(content["iccs_checked"])
        )

    sorted_data = [x for x in data.items()]
    sorted_data = sort_datapoints(sorted_data)
    cmap = mpl.colormaps["gist_rainbow"]
    colors = cmap(np.linspace(0, 1, len(sorted_data)))

    best: None | tuple[str, float] = None
    for (name, points), color in zip(sorted_data, colors):

        log_e = np.array(points["log_e"])
        iccs = np.array(points["iccs"])

        best_log_e = log_e.max()
        if not best or best_log_e > best[1]:
            best = (name, best_log_e)

        plt.plot(iccs, log_e, "o", color=color, ms=9, mec="k", lw=2)

    plt.xlabel("ICCs calculated")
    plt.ylabel("log(E)")
    plt.title("Constructive beam size")
    if save:
        plt.savefig(f"./results/plots/construct_beam_{asctime()}.svg")
    plt.show()
    print(f"best: {best[0]} with {best[1]}")  # type: ignore


@app.command()
def adaptive(save: bool = False) -> None:
    folder = "./results/adaptive_annealing"
    data: dict[str, dict[str, list[float | int]]] = {}

    best: None | tuple[str, float] = None
    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        name = filename[7:-14]
        best_log_e = float(content["log_e"])
        data.setdefault(name, {}).setdefault("log_e", []).append(best_log_e)
        data.setdefault(name, {}).setdefault("iccs", []).append(
            int(content["iccs_checked"])
        )
        if not best or best_log_e > best[1]:
            best = (filename[7:-5], best_log_e)

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

    for name, points in sorted_data:
        log_e = np.array(points["log_e"])
        iccs = np.array(points["iccs"])
        plt.plot(iccs, log_e, shapes(name), color=colors(name), alpha=0.75)

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
    plt.xlabel("ICCs calculated")
    plt.ylabel("log(E)")
    plt.title("Adaptive annealing")
    if save:
        plt.savefig(f"./results/plots/adaptive_annealing_{asctime()}.svg")
    plt.show()

    print(f"best: {best[0]} with {best[1]}")  # type: ignore


@app.command()
def parallel_tempering(save: bool = False) -> None:
    folder = "./results/parallel_tempering"
    data: dict[str, dict[str, list[float | int]]] = {}

    best: None | tuple[str, float] = None
    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        name = filename[7:-14]
        best_log_e = float(content["log_e"])
        data.setdefault(name, {}).setdefault("log_e", []).append(best_log_e)
        data.setdefault(name, {}).setdefault("iccs", []).append(
            int(content["iccs_checked"])
        )
        if not best or best_log_e > best[1]:
            best = (filename[7:-5], best_log_e)

    print(f"best: {best[0]} with {best[1]}")  # type: ignore

    sorted_data = [x for x in data.items()]
    sorted_data = sorted(sorted_data, key=lambda x: x[0])

    colors = MapToColor(
        ["32", "64", "128", "256", "512", "1024"],
        plt.cm.tab10.colors[0:6],  # type: ignore
        ["32", "64", "128", "256", "512", "1024"],
    )

    shapes = MapToShape(
        ["e", "g", "i", "l"],
        ["o", "s", "^", "*"],
        ["exponential", "geometric", "inverse linear", "linear"],
    )

    for name, points in sorted_data:
        log_e = np.array(points["log_e"])
        iccs = np.array(points["iccs"])
        plt.plot(iccs, log_e, shapes(name), color=colors(name), alpha=0.75)

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
    plt.xlabel("ICCs calculated")
    plt.ylabel("log(E)")
    plt.title("Parallel tempering swap amount")
    if save:
        plt.savefig(f"./results/plots/parallel_tempering_{asctime()}.svg")
    plt.show()


@app.command()
def all(save: bool = False):
    prove_greedy_time_dependence(save)
    prove_annealing_time_dependence(save)
    prove_parallel_time_dependence(save)

    greedy(save)
    tabu(save)
    tabu_only(save)
    merge_beam(save)
    construct_beam(save)

    plot_sa(save)
    sa_exp_decay(save)
    sa_shapes(save)

    lahc_hist_size(save)
    adaptive(save)
    parallel_tempering(save)


if __name__ == "__main__":
    app()
