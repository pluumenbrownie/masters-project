import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator
import matplotlib as mpl
import numpy as np
import typer
from time import asctime
import json
import os 
from scipy.stats import linregress # type: ignore

plt.rcParams["figure.constrained_layout.use"] = True

app = typer.Typer()


@app.command()
def prove_time_dependence(save: bool = False):
    folder = "./results/simulated_annealing"
    time_seconds = []
    iccs_covered = []

    for filename in filter(lambda name: name[-4:] == "json", os.listdir(folder)):
        with open(folder + "/" + filename, "r") as file:
            content = json.load(file)
        time_seconds.append(float(content["time_seconds"]))
        iccs_covered.append(int(content["iccs_checked"]))

    x_range = [20000, 22500]
    fit = linregress(iccs_covered, time_seconds)
    fitline = lambda x: fit.slope * x + fit.intercept
    
    plt.plot(iccs_covered, time_seconds, "o", label="SA data")
    plt.plot(x_range, list(map(fitline, x_range)), "--", label="Linear fit")

    plt.xlabel("ICCs calculated")
    plt.ylabel("time (s)")
    plt.title("Runtime over ICCs calculated")
    plt.legend()
    plt.xlim(*x_range)
    plt.show()


@app.command()
def plot_sa(save: bool = True):
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


if __name__ == "__main__":
    app()
