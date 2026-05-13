#!/usr/bin/env python
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MaxNLocator


PLOTS_DIR = Path("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/plots_MVA/run3")
MODEL_FILE = Path("./model_Za_BDT_run3.pkl")


def savefig_and_show(pdf_name, dpi=300):
    out = PLOTS_DIR / pdf_name
    out.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.gcf()
    fig.savefig(str(out), format="pdf", bbox_inches="tight", dpi=dpi)
    plt.close(fig)


def plot_loss_vs_n_estimators(model, pdf_name="loss_vs_nEstimators.pdf"):
    try:
        evals_result = model.evals_result()
    except Exception as exc:
        raise RuntimeError(
            "This model does not contain XGBoost eval history. "
            "Train once with eval_set=... before replotting logloss."
        ) from exc

    if not evals_result:
        raise RuntimeError("This model has an empty XGBoost eval history.")

    dataset_labels = {
        "validation_0": "Train",
        "validation_1": "Test",
    }
    dataset_colors = {
        "validation_0": "darkorange",
        "validation_1": "darkgreen",
    }

    fig, ax = plt.subplots(figsize=(8, 6))
    plotted = False
    for dataset_name, metrics in evals_result.items():
        if "logloss" not in metrics:
            continue

        losses = np.asarray(metrics["logloss"], dtype=float)
        n_estimators = np.arange(1, len(losses) + 1)
        label = dataset_labels.get(dataset_name, dataset_name)
        color = dataset_colors.get(dataset_name)
        ax.plot(n_estimators, losses, color=color, lw=2, label=label)
        plotted = True

        if len(losses) > 0 and np.any(np.isfinite(losses)):
            best_idx = int(np.nanargmin(losses))
            ax.scatter([best_idx + 1], [losses[best_idx]], color=color, s=45, zorder=3)
            print("{} minimum logloss: nEstimators {}, logloss {:.6g}".format(
                label,
                best_idx + 1,
                losses[best_idx],
            ))

    if not plotted:
        raise RuntimeError("The model eval history does not contain logloss.")

    ax.tick_params(axis='x', labelsize=18, direction='in', top=True, right=True, length=11)
    ax.tick_params(axis='y', labelsize=18, direction='in', top=True, right=True, length=11)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(which='minor', direction='in', top=True, right=True, length=5)
    fig.subplots_adjust(left=0.17, right=0.96, top=0.91, bottom=0.17)
    ax.set_xlabel('nEstimators', fontsize=20, labelpad=12)
    ax.set_ylabel('Log Loss', fontsize=20, labelpad=12)
    ax.legend(loc="upper right", fontsize=17, frameon=False)
    savefig_and_show(pdf_name)
    print("Saved:", PLOTS_DIR / pdf_name)


if __name__ == "__main__":
    model = pickle.load(open(MODEL_FILE, "rb"))
    plot_loss_vs_n_estimators(model)
