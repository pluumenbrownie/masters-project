import csv
from functools import cache
from typing import Any
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Rectangle
import numpy as np
import typer
from time import asctime
import json
import os
import re
from scipy.stats import linregress  # type: ignore
from itertools import cycle, chain

from plot_data import par_temp_pools, par_temp_log_e, log_e_hill_climbing, log_e_annealing

plt.rcParams["figure.constrained_layout.use"] = True
plt.rcParams["figure.figsize"] = (10, 8)

EXPORT_FORMAT = "svg"

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
    print(fit)
    fitline = lambda x: fit.slope * x + fit.intercept

    plt.plot(iccs_covered, time_seconds, "o", label="datapoints")
    plt.plot(x_range, list(map(fitline, x_range)), "--", label="linear fit")

    plt.xlabel("ICCs calculated")
    plt.ylabel("time (s)")
    plt.title(f"{name} time dependence")
    plt.legend()
    plt.xlim(*x_range)
    if save:
        plt.savefig(
            f"./results/plots/{name}_time_dependence_{asctime()}.{EXPORT_FORMAT}"
        )
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
        plt.savefig(f"./results/plots/simulated_annealing_{asctime()}.{EXPORT_FORMAT}")
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
    add_runtime_axis("simulated_annealing")
    if save:
        plt.savefig(f"./results/plots/annealing_shapes_{asctime()}.{EXPORT_FORMAT}")
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
    add_runtime_axis("simulated_annealing")
    if save:
        plt.savefig(f"./results/plots/annealing_exp_decay_{asctime()}.{EXPORT_FORMAT}")
    plt.show()


@app.command()
def lahc_hist_size(save: bool = False) -> None:
    folder = "./results/lahc_hist_size"
    data: dict[str, dict[str, list[float | int]]] = {}

    best: None | tuple[str, float] = None
    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        history_size = int(find_floats(filename)[1])

        log_e = content["log_e"]
        if not best or log_e > best[1]:
            best = (filename, log_e)

        if "_a_" in filename:
            data.setdefault(f"{history_size}_a", {}).setdefault("log_e", []).append(
                float(log_e)
            )
            data.setdefault(f"{history_size}_a", {}).setdefault("iccs", []).append(
                int(content["iccs_checked"])
            )
        else:
            data.setdefault(str(history_size), {}).setdefault("log_e", []).append(
                float(log_e)
            )
            data.setdefault(str(history_size), {}).setdefault("iccs", []).append(
                int(content["iccs_checked"])
            )
    assert best
    print(f"best: {best[0]} with {best[1]}")  

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
    add_runtime_axis("simulated_annealing")
    if save:
        plt.savefig(f"./results/plots/lahc_hist_size_{asctime()}.{EXPORT_FORMAT}")
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
        plt.savefig(f"./results/plots/lahc_loglog_{asctime()}.{EXPORT_FORMAT}")
    plt.show()

    log_e_hill_climbing(filename=f"./results/lahc_hist_size/log_e_{best[0][7:-5]}.csv", show=False)
    if save:
        plt.savefig(f"./results/plots/lahc_hist_best_log_e_{asctime()}.{EXPORT_FORMAT}")
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

    best: None | tuple[str, float] = None
    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        name = filename[15:-5]

        best_log_e = float(content["log_e"])
        if not best or best_log_e > best[1]:
            best = (name, best_log_e)
        data.setdefault(name, {}).setdefault("log_e", []).append(best_log_e)
        data.setdefault(name, {}).setdefault("iccs", []).append(
            int(content["iccs_checked"])
        )
    print(f"best: {best[0]} with {best[1]}")  # type: ignore

    sorted_data = [x for x in data.items()]
    sorted_data = sort_datapoints(sorted_data)
    colors = MapToColor(
        ["n_", "m", "c"], plt.cm.tab10.colors[0:3], ["None", "Merge", "Constructive"]  # type: ignore
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
    add_runtime_axis("construct_beam")
    if save:
        plt.savefig(f"./results/plots/greedy_{asctime()}.{EXPORT_FORMAT}")
    plt.show()


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
    add_runtime_axis("construct_beam")
    if save:
        plt.savefig(f"./results/plots/tabu_{asctime()}.{EXPORT_FORMAT}")
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
    add_runtime_axis("construct_beam")
    if save:
        plt.savefig(f"./results/plots/tabu_trajectory_{asctime()}.{EXPORT_FORMAT}")
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
    add_runtime_axis("construct_beam")
    if save:
        plt.savefig(f"./results/plots/tabu2_{asctime()}.{EXPORT_FORMAT}")
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
    add_runtime_axis("construct_beam")
    if save:
        plt.savefig(f"./results/plots/tabu2_trajectory_{asctime()}.{EXPORT_FORMAT}")
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
    add_runtime_axis("construct_beam")
    if save:
        plt.savefig(f"./results/plots/tabu_only_{asctime()}.{EXPORT_FORMAT}")
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
    add_runtime_axis("construct_beam")
    if save:
        plt.savefig(f"./results/plots/tabu_only_trajectory_{asctime()}.{EXPORT_FORMAT}")
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
    add_runtime_axis("construct_beam")
    if save:
        plt.savefig(f"./results/plots/merge_beam_{asctime()}.{EXPORT_FORMAT}")
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
    add_runtime_axis("construct_beam")
    if save:
        plt.savefig(f"./results/plots/construct_beam_{asctime()}.{EXPORT_FORMAT}")
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
    assert best

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
    add_runtime_axis("simulated_annealing")
    if save:
        plt.savefig(f"./results/plots/adaptive_annealing_{asctime()}.{EXPORT_FORMAT}")
    plt.show()

    print(f"best: {best[0]} with {best[1]}")  

    log_e_annealing(filename=f"./results/adaptive_annealing/log_e_{best[0]}.csv", show=False)
    if save:
        plt.savefig(f"./results/plots/adaptive_annealing_best_log_e_{asctime()}.{EXPORT_FORMAT}")
    plt.show()


STEPS = {"-1": 8192, "-0": 16384, "+1": 32768, "+2": 65536}


def shuffle_tag(name: str) -> str | None:
    """Extracts '-1', '-0', '+1' or '+2' from a group name."""
    match = re.search(r"[-+]\d", name)
    return match.group(0) if match and match.group(0) in STEPS else None


def draw_group_boxes(
    boxes: dict[str, tuple[list[float], list[float]]], pad: float = 0.02
) -> None:
    """Dashed box + label around all points of each shuffle group."""
    ax = plt.gca()
    all_x = [x for xs, _ in boxes.values() for x in xs]
    all_y = [y for _, ys in boxes.values() for y in ys]
    pad_x = (max(all_x) - min(all_x)) * pad
    pad_y = (max(all_y) - min(all_y)) * pad

    for tag, (xs, ys) in sorted(boxes.items(), key=lambda item: STEPS[item[0]]):
        x0, x1 = min(xs) - pad_x, max(xs) + pad_x
        y0, y1 = min(ys) - pad_y, max(ys) + pad_y
        ax.add_patch(
            Rectangle(
                (x0, y0),
                x1 - x0,
                y1 - y0,
                fill=False,
                linestyle="--",
                edgecolor="k",
                linewidth=1.2,
                zorder=0,
            )
        )
        ax.text(
            (x0 + x1) / 2,
            y1 + pad_y * 0.5,
            f"{STEPS[tag]:,} steps",
            ha="center",
            va="bottom",
            fontsize="small",
        )
    ax.margins(y=0.1)  # leave room for the labels
    ax.autoscale_view()


@app.command()
def parallel_tempering(save: bool = False) -> None:
    folder = "./results/parallel_tempering"
    data: dict[str, dict[str, list[float | int]]] = {}

    best: None | tuple[str, float, str] = None
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
            best = (filename[7:-5], best_log_e, filename)
    assert best

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

    boxes: dict[str, tuple[list[float], list[float]]] = {}
    for name, points in sorted_data:
        log_e = np.array(points["log_e"])
        iccs = np.array(points["iccs"])

        tag = shuffle_tag(name)
        if tag is not None:
            xs, ys = boxes.setdefault(tag, ([], []))
            xs.extend(iccs.tolist())
            ys.extend(log_e.tolist())

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
    draw_group_boxes(boxes)
    colors.plot_legend()
    shapes.plot_legend()

    plt.legend()
    plt.xlabel("ICCs calculated")
    plt.ylabel("log(E)")
    plt.title("Parallel tempering swap amount")
    add_runtime_axis("parallel_tempering_timing")
    if save:
        plt.savefig(f"./results/plots/parallel_tempering_{asctime()}.{EXPORT_FORMAT}")
    plt.show()

    par_temp_pools(
        filepath_ids=f"./results/parallel_tempering/mcm_id_{best[0]}.csv",
        filepath_log_e=f"./results/parallel_tempering/log_e_{best[0]}.csv",
        show=False,
    )
    if save:
        plt.savefig(
            f"./results/plots/parallel_tempering_best_pools_{asctime()}.{EXPORT_FORMAT}"
        )
    plt.show()
    par_temp_log_e(
        filepath_ids=f"./results/parallel_tempering/mcm_id_{best[0]}.csv",
        filepath_log_e=f"./results/parallel_tempering/log_e_{best[0]}.csv",
        show=False,
    )
    if save:
        plt.savefig(
            f"./results/plots/parallel_tempering_best_log_e_{asctime()}.{EXPORT_FORMAT}"
        )
    plt.show()


# dataset folder -> (filename-to-group-name function, folder used for the runtime fit)
PERCENTILE_DATASETS: dict[str, tuple[Any, str]] = {
    "greedy": (lambda f: f[15:-5], "construct_beam"),
    "simulated_annealing_shapes": (lambda f: f[7:-14], "simulated_annealing"),
    "simulated_annealing_longer": (
        lambda f: str(find_floats(f)[1]),
        "simulated_annealing",
    ),
    "lahc_hist_size": (
        lambda f: str(int(find_floats(f)[1])) + ("_a" if "_a_" in f else ""),
        "simulated_annealing",
    ),
    "adaptive_annealing": (lambda f: f[7:-14], "simulated_annealing"),
    "parallel_tempering": (lambda f: f[7:-14], "parallel_tempering_timing"),
}


def load_time_fit(folder: str):
    """Same linear fit (time_seconds vs ICCs) as prove_time_dependence."""
    time_seconds = []
    iccs_covered = []
    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        time_seconds.append(float(content["time_seconds"]))
        iccs_covered.append(int(content["iccs_checked"]))
    return linregress(iccs_covered, time_seconds)


@cache
def cached_time_fit(time_folder: str):
    return load_time_fit(f"./results/{time_folder}")


def add_runtime_axis(time_folder: str) -> None:
    """Adds a top x-axis with the extrapolated runtime to the current axes."""
    ax = plt.gca()
    fit = cached_time_fit(time_folder)
    secax = ax.secondary_xaxis(
        "top",
        functions=(
            lambda iccs: fit.slope * iccs + fit.intercept,
            lambda seconds: (seconds - fit.intercept) / fit.slope,
        ),
    )
    secax.set_xlabel("extrapolated runtime (s)")
    secax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.0f}"))

    # full y numbers instead of "-2.38" + "1e6" offset text, which collides with the top axis
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.0f}"))
    ax.yaxis.get_offset_text().set_visible(False)


def collect_percentile_points(
    folder: str, namer
) -> list[tuple[str, np.ndarray, np.ndarray]]:
    """Returns, per group in `folder`, the (iccs, log_e) 25/50/75 percentiles."""
    data: dict[str, dict[str, list[float | int]]] = {}
    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        name = namer(filename)
        data.setdefault(name, {}).setdefault("log_e", []).append(
            float(content["log_e"])
        )
        data.setdefault(name, {}).setdefault("iccs", []).append(
            int(content["iccs_checked"])
        )

    return [
        (
            name,
            np.percentile(points["iccs"], [25, 50, 75]),
            np.percentile(points["log_e"], [25, 50, 75]),
        )
        for name, points in data.items()
    ]


def plot_percentiles_all_datasets(
    save: bool = False,
    ylim: tuple[float, float] | None = None,
    suffix: str = "",
    runtime: bool = False,
) -> None:
    base = "./results"

    # greedy keeps its own colors (tab10[0:3]) and shapes, as in the greedy plot
    greedy_colors = MapToColor(
        ["n_", "m", "c"],
        plt.cm.tab10.colors[0:3],  # type: ignore
        ["Greedy: None", "Greedy: Merge", "Greedy: Constructive"],
    )
    greedy_shapes = MapToShape(
        ["t", "l", "nr"],
        ["^", "s", "o"],
        ["Greedy: Tabu", "Greedy: Local", "Greedy: None"],
    )

    # other datasets use the remaining tab10 colors so they don't clash with greedy
    other_colors = iter(plt.cm.tab10.colors[3:])  # type: ignore
    dataset_colors = {
        dataset: (None if dataset == "greedy" else next(other_colors))
        for dataset in PERCENTILE_DATASETS
    }

    fits: dict[str, Any] = {}
    greedy_plotted = False
    for dataset, (namer, time_folder) in PERCENTILE_DATASETS.items():
        is_greedy = dataset == "greedy"
        folder = f"{base}/{dataset}"
        if not os.path.isdir(folder):
            print(f"skipping {dataset}: folder not found")
            continue

        if runtime:
            if time_folder not in fits:
                time_path = f"{base}/{time_folder}"
                if not os.path.isdir(time_path):
                    print(f"skipping {dataset}: timing folder {time_folder} not found")
                    continue
                fits[time_folder] = load_time_fit(time_path)
            fit = fits[time_folder]

        for i, (name, iccs_p, log_e_p) in enumerate(
            collect_percentile_points(folder, namer)
        ):
            if is_greedy:
                point_color = greedy_colors(name)
                point_shape = greedy_shapes(name)
                label = None  # greedy is covered by the color/shape legend below
                greedy_plotted = True
            else:
                point_color = dataset_colors[dataset]
                point_shape = "o"
                label = dataset if i == 0 else None

            if runtime:
                x_p = np.sort(fit.slope * iccs_p + fit.intercept)
            else:
                x_p = iccs_p

            plt.errorbar(
                x_p[1],
                log_e_p[1],
                xerr=[[x_p[1] - x_p[0]], [x_p[2] - x_p[1]]],
                yerr=[[log_e_p[1] - log_e_p[0]], [log_e_p[2] - log_e_p[1]]],
                fmt=point_shape,
                color=point_color,
                ms=9 if is_greedy else 7,
                capsize=3,
                lw=2 if is_greedy else 1.5,
                alpha=1.0 if is_greedy else 0.85,
                label=label,
            )

    if greedy_plotted:
        greedy_colors.plot_legend()
        greedy_shapes.plot_legend()

    plt.legend(ncol=2, fontsize="small")
    plt.xlabel("extrapolated runtime (s)" if runtime else "ICCs calculated")
    plt.ylabel("log(E)")
    plt.title("Percentile points of all datasets")
    if ylim is not None:
        plt.ylim(ylim)
    if save:
        file_name = "percentiles_runtime" if runtime else "percentiles_all"
        plt.savefig(f"./results/plots/{file_name}{suffix}_{asctime()}.{EXPORT_FORMAT}")
    plt.show()


ZOOM_YLIM = (-2_440_000, -2_380_000)


@app.command()
def percentiles_all(save: bool = False) -> None:
    plot_percentiles_all_datasets(save)


@app.command()
def percentiles_all_zoomed(save: bool = False) -> None:
    plot_percentiles_all_datasets(save, ylim=ZOOM_YLIM, suffix="_zoomed")


@app.command()
def percentiles_runtime(save: bool = False, zoomed: bool = False) -> None:
    plot_percentiles_all_datasets(
        save,
        ylim=ZOOM_YLIM if zoomed else None,
        suffix="_zoomed" if zoomed else "",
        runtime=True,
    )


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

    percentiles_all(save)
    percentiles_all_zoomed(save)
    percentiles_runtime(save)
    percentiles_runtime(save, zoomed=True)


if __name__ == "__main__":
    app()
