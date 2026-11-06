"""Plot experiment results with shared loading, styling and summary helpers.

Existing Typer commands are preserved. Combined-plot commands:
    python plot_agregate.py combined-percentiles --save
    python plot_agregate.py combined-percentiles-limited --save

The combined plot contains one marginal median/IQR summary per configuration
from every percentile-based experiment, plus all raw results/greedy points.
The limited version retains the ICCs x-axis and limits the log(E) y-axis.
Filename parsing deliberately retains the original experiment conventions.
Requires Python 3.10+, matplotlib, numpy, scipy and typer.
"""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from datetime import datetime
from itertools import cycle
from pathlib import Path
from typing import Any, Callable, Iterable, Literal

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import typer
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from scipy.stats import linregress

plt.rcParams.update({
    "figure.constrained_layout.use": True,
    "figure.figsize": (10, 8),
})
app = typer.Typer(help=__doc__)
RESULTS_DIR = Path("results")
COMBINED_LOG_E_LIMITS = (-2_470_000, -2_380_000)
Key = str | float | int
TAB10 = mpl.colormaps["tab10"].colors
FLOAT_PATTERN = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


@dataclass(frozen=True)
class Samples:
    iccs: np.ndarray
    log_e: np.ndarray


@dataclass(frozen=True)
class StyleMap:
    """Ordered substring rules shared by both color and marker mappings."""

    rules: tuple[tuple[str, Any, str], ...]

    def __call__(self, name: Key) -> Any:
        for key, value, _ in self.rules:
            if key in str(name):
                return value
        raise ValueError(f"No style rule matches configuration {name!r}")

    def handles(self, kind: Literal["color", "marker"], prefix: str = "") -> list[Line2D]:
        return [
            Line2D(
                [], [], linestyle="none",
                marker="s" if kind == "color" else value,
                color=value if kind == "color" else "black",
                markerfacecolor=value if kind == "color" else "white",
                markersize=9, label=f"{prefix}{label}",
            )
            for _, value, label in self.rules
        ]


def find_floats(text: str) -> list[float]:
    return [float(value) for value in FLOAT_PATTERN.findall(text)]


def numeric_sort_key(name: Key) -> tuple[bool, float, str]:
    numbers = find_floats(str(name))
    return (not bool(numbers), numbers[0] if numbers else 0.0, str(name))


def filename_key(filename: str) -> str:
    return filename[7:-5]


def repeated_run_key(filename: str) -> str:
    return filename[7:-14]


def parameter_key(filename: str) -> float:
    numbers = find_floats(filename)
    if len(numbers) < 2:
        raise ValueError(f"Expected a parameter as the second number in {filename!r}")
    return numbers[1]


def annealing_marker(name: Key) -> str:
    return next((marker for key, marker in (("exp", "^"), ("lin", "s"), ("cus", "x"))
                 if key in str(name)), "o")


@dataclass(frozen=True)
class Experiment:
    folder: str
    title: str
    output: str
    key: Callable[[str], Key] = filename_key
    sort: Literal["numeric", "alphabetical"] = "numeric"
    percentiles: bool = False
    palette: str = "tab10"
    color_map: StyleMap | None = None
    marker_map: StyleMap | None = None
    marker: Callable[[Key], str] | None = None
    legend_columns: int = 1
    legend: bool = True


GREEDY_COLORS = StyleMap(tuple(zip(
    ("n", "m", "c"), TAB10[:3], ("None", "Merge", "Constructive"),
)))
GREEDY_MARKERS = StyleMap((("t", "^", "Tabu"), ("l", "s", "Local"), ("nr", "o", "None")))
ADAPTIVE_COLORS = StyleMap(tuple((key, color, key) for key, color in zip(
    ("0.004", "0.003", "0.002", "0.0015", "0.001", "0.0007", "0.0005"), TAB10,
)))
ADAPTIVE_MARKERS = StyleMap(tuple(zip(
    ("0.95", "0.9", "0.85", "0.8", "0.75"), ("o", "s", "^", "*", "P"),
    ("0.95", "0.90", "0.85", "0.80", "0.75"),
)))
PARALLEL_COLORS = StyleMap(tuple((key, color, key) for key, color in zip(
    ("32", "64", "128", "256", "512", "1024"), TAB10,
)))
PARALLEL_MARKERS = StyleMap(tuple(zip(
    ("e", "g", "i", "l"), ("o", "s", "^", "*"),
    ("exponential", "geometric", "inverse linear", "linear"),
)))

EXPERIMENTS = {
    "greedy": Experiment("greedy", "Greedy Strategy Comparison", "greedy",
                         key=lambda name: name[15:-5], color_map=GREEDY_COLORS,
                         marker_map=GREEDY_MARKERS),
    "tabu": Experiment("tabu_divergence", "Tabu list size", "tabu", legend_columns=2),
    "tabu_2": Experiment("tabu_divergence_lower", "Tabu list size (lower range)", "tabu2",
                         legend_columns=2),
    "tabu_only": Experiment("pure_tabu", "Tabu only strategy", "tabu_only",
                            palette="gist_rainbow", legend=False),
    "merge_beam": Experiment("merge_beam", "Merge beam size", "merge_beam",
                             palette="gist_rainbow", legend=False),
    "construct_beam": Experiment("construct_beam", "Constructive beam size", "construct_beam",
                                 palette="gist_rainbow", legend=False),
    "sa_shapes": Experiment("simulated_annealing_shapes", "Annealing shapes comparison",
                            "annealing_shapes", key=repeated_run_key, sort="alphabetical",
                            percentiles=True, marker=annealing_marker, legend_columns=2),
    "sa_exp_decay": Experiment("simulated_annealing_longer", "Exponential decay speed",
                               "annealing_exp_decay", key=parameter_key, percentiles=True),
    "lahc_hist_size": Experiment("lahc_hist_size", "Late Acceptance Hill Climbing",
                                 "lahc_hist_size", key=lambda name: int(parameter_key(name)),
                                 percentiles=True),
    "adaptive": Experiment("adaptive_annealing", "Adaptive annealing", "adaptive_annealing",
                           key=repeated_run_key, sort="alphabetical", percentiles=True,
                           color_map=ADAPTIVE_COLORS, marker_map=ADAPTIVE_MARKERS),
    "parallel_tempering": Experiment("parallel_tempering", "Parallel tempering swap amount",
                                     "parallel_tempering", key=repeated_run_key,
                                     sort="alphabetical", percentiles=True,
                                     color_map=PARALLEL_COLORS, marker_map=PARALLEL_MARKERS),
}


# Loading and statistics are independent of figure creation.
def result_files(folder: str | Path, suffix: str) -> list[Path]:
    path = Path(folder)
    if not path.is_dir():
        raise FileNotFoundError(f"Results directory does not exist: {path}")
    return sorted(file for file in path.glob(f"*{suffix}") if file.is_file())


def json_records(folder: str | Path) -> Iterable[tuple[Path, dict[str, Any]]]:
    files = result_files(folder, ".json")
    if not files:
        raise ValueError(f"No JSON results in {folder}")
    for path in files:
        with path.open(encoding="utf-8") as file:
            yield path, json.load(file)


def load_experiment(spec: Experiment) -> dict[Key, Samples]:
    grouped: dict[Key, list[tuple[int, float]]] = {}
    for path, record in json_records(RESULTS_DIR / spec.folder):
        try:
            name = spec.key(path.name)
            pair = (int(record["iccs_checked"]), float(record["log_e"]))
            if not np.isfinite(pair).all():
                raise ValueError("Non-finite ICC count or log(E)")
        except (KeyError, ValueError, TypeError, IndexError) as error:
            raise ValueError(f"Invalid result in {path}: {error}") from error
        grouped.setdefault(name, []).append(pair)
    names = sorted(grouped, key=numeric_sort_key if spec.sort == "numeric" else str)
    return {
        name: Samples(
            iccs=np.array([pair[0] for pair in grouped[name]], dtype=np.int64),
            log_e=np.array([pair[1] for pair in grouped[name]], dtype=float),
        )
        for name in names
    }


def percentile_summary(samples: Samples) -> tuple[np.ndarray, np.ndarray]:
    """Return marginal medians and asymmetric 25th/75th percentile distances.

    Medians are computed independently for ICCs and log(E), as in the original
    plots. They need not correspond to the same observed run. Bars are IQRs,
    not confidence intervals. Shape of errors is (coordinate, lower/upper).
    """
    quartiles = np.percentile(np.vstack((samples.iccs, samples.log_e)), [25, 50, 75], axis=1)
    median = quartiles[1]
    errors = np.stack((median - quartiles[0], quartiles[2] - median), axis=1)
    return median, errors


def best_result(data: dict[Key, Samples]) -> tuple[Key, float]:
    return max(((name, float(samples.log_e.max())) for name, samples in data.items()),
               key=lambda item: item[1])


def csv_columns(path: Path, columns: tuple[int, ...], *, exclude_constructive: bool = False
                ) -> tuple[np.ndarray, ...]:
    values: list[list[float]] = [[] for _ in columns]
    with path.open(newline="", encoding="utf-8") as file:
        reader = csv.reader(file)
        next(reader, None)
        for line_number, row in enumerate(reader, start=2):
            if not row or (exclude_constructive and row[0] == "Constructive"):
                continue
            try:
                parsed = [float(row[column]) for column in columns]
            except (IndexError, ValueError) as error:
                raise ValueError(f"Invalid CSV row in {path}, line {line_number}") from error
            for target, value in zip(values, parsed):
                target.append(value)
    return tuple(np.asarray(value, dtype=float) for value in values)


# Shared rendering helpers always use explicit axes and figures.
def styled_groups(spec: Experiment, data: dict[Key, Samples]
                  ) -> Iterable[tuple[Key, Samples, Any, str]]:
    if spec.palette == "tab10":
        palette = cycle(TAB10)  # Never truncate configurations after ten colors.
    else:
        palette = iter(mpl.colormaps[spec.palette](np.linspace(0, 1, len(data))))
    for name, samples in data.items():
        color = spec.color_map(name) if spec.color_map else next(palette)
        marker = spec.marker_map(name) if spec.marker_map else (
            spec.marker(name) if spec.marker else "o")
        yield name, samples, color, marker


def draw_raw(ax: Axes, samples: Samples, color: Any, marker: str, *,
             label: str | None = None, summary_overlay: bool = False) -> None:
    ax.plot(samples.iccs, samples.log_e, linestyle="none", marker=marker,
            color=color, label=label, alpha=0.75 if summary_overlay else 1.0,
            markersize=6 if summary_overlay else 9,
            markeredgecolor=color if summary_overlay else "black", linewidth=2, zorder=3)


def draw_percentiles(ax: Axes, samples: Samples, color: Any, marker: str, *,
                     label: str | None = None) -> np.ndarray:
    median, errors = percentile_summary(samples)
    ax.errorbar(median[0], median[1], xerr=errors[0].reshape(2, 1),
                yerr=errors[1].reshape(2, 1), fmt=marker, color=color, label=label,
                markersize=9, markeredgecolor="black", capsize=4, linewidth=2, zorder=4)
    return median


def add_legend(ax: Axes, spec: Experiment) -> None:
    handles = []
    if spec.color_map:
        handles.extend(spec.color_map.handles("color"))
    if spec.marker_map:
        handles.extend(spec.marker_map.handles("marker"))
    if spec.legend:
        if handles:
            ax.legend(handles=handles, ncol=spec.legend_columns)
        elif ax.get_legend_handles_labels()[0]:
            ax.legend(ncol=spec.legend_columns)


def finish_figure(fig: Figure, output: str, save: bool) -> None:
    try:
        if save:
            folder = RESULTS_DIR / "plots"
            folder.mkdir(parents=True, exist_ok=True)
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
            fig.savefig(folder / f"{output}_{timestamp}.svg")
        plt.show()
    finally:
        plt.close(fig)


def label_axes(ax: Axes, title: str, xlabel: str = "ICCs calculated",
               ylabel: str = "log(E)") -> None:
    ax.set(title=title, xlabel=xlabel, ylabel=ylabel)


def plot_experiment(name: str, save: bool = False) -> tuple[Experiment, dict[Key, Samples]]:
    spec = EXPERIMENTS[name]
    data = load_experiment(spec)
    fig, ax = plt.subplots()
    for key, samples, color, marker in styled_groups(spec, data):
        label = None if spec.color_map or spec.marker_map else str(key)
        draw_raw(ax, samples, color, marker, label=label, summary_overlay=spec.percentiles)
        if spec.percentiles:
            draw_percentiles(ax, samples, color, marker)
    add_legend(ax, spec)
    label_axes(ax, spec.title)
    finish_figure(fig, spec.output, save)
    best_name, value = best_result(data)
    print(f"best: {best_name} with {value}")
    return spec, data


def plot_trajectories(spec: Experiment, save: bool, *, exclude_constructive: bool = True,
                      best_name: Key | None = None,
                      ylim: tuple[float, float] = (-2430000, -2380000)) -> None:
    paths = result_files(RESULTS_DIR / spec.folder, ".csv")
    if not paths:
        typer.echo(f"Skipping trajectories: no CSV files in {spec.folder}", err=True)
        return
    paths.sort(key=lambda path: numeric_sort_key(path.name[8:-4]))
    palette = cycle(tuple(mpl.colormaps["tab20c"].colors) + tuple(mpl.colormaps["tab20b"].colors))
    fig, ax = plt.subplots()
    for path in paths:
        color = next(palette)
        name = path.name[6:-4]
        (values,) = csv_columns(path, (2,), exclude_constructive=exclude_constructive)
        is_best = best_name is not None and name == str(best_name)
        label = f"best: {name}" if is_best else (name if best_name is None else None)
        ax.plot(np.arange(len(values)), values, color="black" if is_best else color, label=label)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(ncol=4 if best_name is None else 1)
    label_axes(ax, "Tabu only trajectories" if spec.folder == "pure_tabu"
               else "Tabu list size trajectories")
    ax.set_ylim(*ylim)
    finish_figure(fig, f"{spec.output}_trajectory", save)


def prove_time_dependence(save: bool, folder: str | Path, name: str) -> None:
    records = [record for _, record in json_records(folder)]
    iccs = np.array([int(record["iccs_checked"]) for record in records])
    seconds = np.array([float(record["time_seconds"]) for record in records])
    if len(np.unique(iccs)) < 2 or not np.isfinite(seconds).all():
        raise ValueError("Time regression needs finite times and at least two distinct ICC counts")
    fit = linregress(iccs, seconds)
    x_range = np.array([iccs.min() * 0.995, iccs.max() * 1.005])
    fig, ax = plt.subplots()
    ax.plot(iccs, seconds, "o", label="datapoints")
    ax.plot(x_range, fit.slope * x_range + fit.intercept, "--", label="linear fit")
    label_axes(ax, f"{name} time dependence", ylabel="time (s)")
    ax.set_xlim(*x_range)
    ax.legend()
    finish_figure(fig, f"{name}_time_dependence", save)


@app.command()
def prove_annealing_time_dependence(save: bool = False) -> None:
    prove_time_dependence(save, RESULTS_DIR / "simulated_annealing", "Annealing")


@app.command()
def prove_parallel_time_dependence(save: bool = False) -> None:
    prove_time_dependence(save, RESULTS_DIR / "parallel_tempering_timing", "Parallel tempering")


@app.command()
def prove_greedy_time_dependence(save: bool = False) -> None:
    prove_time_dependence(save, RESULTS_DIR / "construct_beam", "Greedy")


@app.command()
def plot_sa(save: bool = False) -> None:
    """Plot all log(E) trajectories and the longest run's temperature schedule.

    A short first CSV must not determine the red curve's duration. Temperatures
    are neither padded nor extrapolated; invalid log-scale values leave gaps.
    """
    paths = result_files(RESULTS_DIR / "simulated_annealing", ".csv")
    if not paths:
        raise ValueError("No simulated-annealing CSV trajectories found")
    fig, ax = plt.subplots()
    temperature_ax = ax.twinx()
    longest_temperature = np.array([], dtype=float)
    temperature_path: Path | None = None
    for path in paths:
        log_e, temperature = csv_columns(path, (0, 1))
        ax.plot(np.arange(len(log_e)), log_e, color="#1f77b4", alpha=0.3)
        if len(temperature) > len(longest_temperature):
            longest_temperature = temperature
            temperature_path = path
    if temperature_path is None:
        plt.close(fig)
        raise ValueError("Simulated-annealing CSV trajectories contain no data rows")

    valid = np.isfinite(longest_temperature) & (longest_temperature > 0)
    if not valid.all():
        typer.echo(f"Omitting {(~valid).sum()} invalid log-scale temperatures in "
                   f"{temperature_path.name}; the curve will contain gaps", err=True)
    if valid.any():
        temperature_ax.set_yscale("log")
        temperature_ax.plot(np.arange(len(longest_temperature)),
                            np.where(valid, longest_temperature, np.nan),
                            color="red", linewidth=2, label="T (longest run)")
        temperature_ax.legend()
    else:
        typer.echo(f"No positive finite temperatures in {temperature_path.name}", err=True)
        temperature_ax.text(0.98, 0.98, "No valid temperature data", ha="right", va="top",
                            transform=temperature_ax.transAxes, color="red")
        temperature_ax.set_yticks([])
    label_axes(ax, "Log(E) of Simulated Annealing", xlabel="step")
    temperature_ax.set_ylabel("Temperature (log scale)", color="red")
    temperature_ax.tick_params(axis="y", colors="red")
    ax.set_xlim(-50, max(len(longest_temperature), 1))
    finish_figure(fig, "simulated_annealing", save)


@app.command()
def greedy(save: bool = False) -> None:
    spec, data = plot_experiment("greedy", save)
    fig, axes = plt.subplots(3, 1)
    for name, samples, _, marker in styled_groups(spec, data):
        if "t" not in str(name):
            continue
        index = 0 if "m" in str(name) else (1 if "c" in str(name) else 2)
        draw_raw(axes[index], samples, None, marker, label=str(name))
    for ax, title in zip(axes, ("Merge", "Constructive", "None")):
        label_axes(ax, title)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(ncol=3)
    fig.suptitle("Greedy tabu refinements")
    finish_figure(fig, "greedy_split", save)


@app.command()
def tabu(save: bool = False) -> None:
    spec, _ = plot_experiment("tabu", save)
    plot_trajectories(spec, save)


@app.command()
def tabu_2(save: bool = False) -> None:
    spec, _ = plot_experiment("tabu_2", save)
    plot_trajectories(spec, save)


@app.command()
def tabu_only(save: bool = False) -> None:
    spec, data = plot_experiment("tabu_only", save)
    plot_trajectories(spec, save, exclude_constructive=False, best_name=best_result(data)[0],
                      ylim=(-2450000, -2380000))


@app.command()
def merge_beam(save: bool = False) -> None:
    plot_experiment("merge_beam", save)


@app.command()
def construct_beam(save: bool = False) -> None:
    plot_experiment("construct_beam", save)


@app.command()
def sa_shapes(save: bool = False) -> None:
    plot_experiment("sa_shapes", save)


@app.command()
def sa_exp_decay(save: bool = False) -> None:
    plot_experiment("sa_exp_decay", save)


@app.command()
def adaptive(save: bool = False) -> None:
    plot_experiment("adaptive", save)


@app.command()
def parallel_tempering(save: bool = False) -> None:
    plot_experiment("parallel_tempering", save)


@app.command()
def lahc_hist_size(save: bool = False) -> None:
    spec, data = plot_experiment("lahc_hist_size", save)
    fig, ax = plt.subplots()
    for name, samples, color, _ in styled_groups(spec, data):
        ranked = np.sort(samples.log_e)
        survival = 1.0 - np.arange(len(ranked)) / len(ranked)
        positive = -ranked > 0
        if not positive.all():
            typer.echo(f"Omitting nonpositive -log(E) values for LAHC {name} on log axes", err=True)
        ax.plot(-ranked[positive], survival[positive], "o-", color=color, label=str(name))
    ax.set(xscale="log", yscale="log")
    label_axes(ax, "LAHC empirical survival", xlabel="-log(E)", ylabel="Survival probability")
    ax.legend()
    finish_figure(fig, "lahc_loglog", save)


def plot_combined_percentiles(save: bool = False, annotate: bool = False, *,
                              ylim: tuple[float, float] | None = None) -> None:
    """Shared renderer for full and limited combined plots.

    Missing/empty experiment folders are reported and skipped. Limits only
    change the view: no samples are filtered and no summaries are recomputed.
    """
    available: dict[str, dict[Key, Samples]] = {}
    for name, spec in EXPERIMENTS.items():
        if not (spec.percentiles or name == "greedy"):
            continue
        folder = RESULTS_DIR / spec.folder
        if not folder.is_dir() or not any(folder.glob("*.json")):
            typer.echo(f"Skipping {name}: no JSON results in {folder}", err=True)
            continue
        available[name] = load_experiment(spec)
    if not available:
        raise typer.BadParameter("No percentile or greedy result files were found")

    fig, ax = plt.subplots(figsize=(12, 9))
    handles: list[Line2D] = []
    markers = cycle(("D", "v", "P", "X", "*", "h"))
    colors = cycle(TAB10[3:])  # Reserve original greedy colors (indices 0, 1, 2).
    for name, spec in EXPERIMENTS.items():
        if not spec.percentiles:
            continue
        color, marker = next(colors), next(markers)
        if name not in available:
            continue
        for configuration, samples in available[name].items():
            median = draw_percentiles(ax, samples, color, marker)
            if annotate:
                ax.annotate(f"{name}: {configuration}", median, xytext=(5, 5),
                            textcoords="offset points", fontsize=8,
                            annotation_clip=True, clip_on=True)
        handles.append(Line2D([], [], linestyle="none", marker=marker, color=color,
                              markeredgecolor="black", markersize=9, label=spec.title))
    if "greedy" in available:
        for _, samples, color, marker in styled_groups(EXPERIMENTS["greedy"], available["greedy"]):
            draw_raw(ax, samples, color, marker)  # No aggregation, jitter, scaling or filtering.
        handles.extend(GREEDY_COLORS.handles("color", "Greedy strategy: "))
        handles.extend(GREEDY_MARKERS.handles("marker", "Greedy refinement: "))
    title = "All percentile summaries with raw greedy points"
    if ylim is not None:
        ax.set_ylim(*ylim)
        title += " (limited log(E) range)"
    label_axes(ax, title)
    ax.legend(handles=handles, loc="best", ncol=2, fontsize=9)
    ax.text(0.01, 0.01, "Summaries: marginal medians with 25th–75th percentile bars; greedy: raw runs",
            transform=ax.transAxes, fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8})
    output = "combined_percentiles" if ylim is None else "combined_percentiles_limited"
    finish_figure(fig, output, save)


@app.command()
def combined_percentiles(save: bool = False, annotate: bool = False) -> None:
    """Combine every configuration's median/IQR with unchanged raw greedy points.

    Use --annotate to label individual configurations; the default legend
    identifies experiment families.
    """
    plot_combined_percentiles(save=save, annotate=annotate)


@app.command()
def combined_percentiles_limited(save: bool = False, annotate: bool = False) -> None:
    """Combined plot with the log(E) y-axis limited to -2470000 through -2380000."""
    plot_combined_percentiles(save=save, annotate=annotate, ylim=COMBINED_LOG_E_LIMITS)


@app.command()
def all(save: bool = False) -> None:
    """Run the original plot suite and both combined-percentiles plots."""
    commands = (
        prove_greedy_time_dependence, prove_annealing_time_dependence,
        prove_parallel_time_dependence, greedy, tabu, tabu_only, merge_beam,
        construct_beam, plot_sa, sa_exp_decay, sa_shapes, lahc_hist_size,
        adaptive, parallel_tempering, combined_percentiles, combined_percentiles_limited,
    )
    for command in commands:
        command(save=save)


if __name__ == "__main__":
    app()
