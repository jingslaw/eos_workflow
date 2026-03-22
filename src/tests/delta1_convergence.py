import json
import numpy as np
import matplotlib.pyplot as plt
from eos_workflow.delta_metric import load_ae_birch_murnaghan


def plot_abs_delta_vs_ecut(
    filename,
    y_range=None,
    y_horizon_lines=None,
):
    """
    Plot |delta(Ecut) - delta(Ecut_max)| vs cutoff energy.

    Units: meV / atom
    """

    # =====================
    # Load data
    # =====================
    with open(filename, "r") as f:
        data = json.load(f)

    element = data.get("element", None)
    config = data.get("configuration", None)
    unit = data.get("delta/natoms unit", "meV/atom")
    ref_v0_b0_b1 = data.get("reference_ae_V0_B0_B1", None)

    if element is None or config is None:
        print(f"some key information is missing in element: {element}, configuration:{config}")
    if ref_v0_b0_b1 is None:
        bm_fit = load_ae_birch_murnaghan(element, config)
        ref_v0_b0_b1 = np.array([bm_fit["min_volume"], bm_fit["bulk_modulus_ev_ang3"], bm_fit["bulk_deriv"],])

    ecut_keys = sorted(
        [k for k in data if k.startswith("ecut-")],
        key=lambda x: int(x.split("-")[1]),
    )

    ecut = np.array([int(k.split("-")[1]) for k in ecut_keys])

    delta = np.array([
        data[k]["delta/natoms"] if data[k]["delta/natoms"] is not None else np.nan
        for k in ecut_keys
    ], dtype=float)

    eos_ok = np.array([data[k]["eos_is_converged"] for k in ecut_keys], dtype=bool)

    v0_b0_b1 = np.array([
        data[k]["v0_b0_b1"] if data[k]["v0_b0_b1"] is not None else [np.nan, np.nan, np.nan]
        for k in ecut_keys
    ], dtype=float)

    v0 = v0_b0_b1[:, 0]
    b0 = v0_b0_b1[:, 1]

    v0_ref, b0_ref = ref_v0_b0_b1[:2]

    delta1 = delta * (v0_ref * b0_ref) / (v0 * b0)

    valid = eos_ok & np.isfinite(delta1)
    delta1[~valid] = np.nan

    # =====================
    # Reference (highest Ecut)
    # =====================
    delta_ref = delta1[-1]
    abs_diff = np.abs(delta1 - delta_ref)

    # =====================
    # Convergence thresholds (meV / atom)
    # =====================
    thresholds = {
        "low": 2,
        "normal": 1,
        "high": 0.5,
    }

    conv_ecut = {}
    for label, thr in thresholds.items():
        mask = abs_diff < thr

        stable = np.logical_and.accumulate(mask[::-1])[::-1]
        idx = np.where(stable)[0]

        conv_ecut[label] = ecut[idx[0]] if len(idx) > 0 else None

    # =====================
    # Plot
    # =====================
    fig, ax = plt.subplots(figsize=(6, 4))

    ax.plot(ecut, abs_diff, "-o", label=r"$|\Delta \delta|$")

    # Mark EOS failures
    for i, ok in enumerate(eos_ok):
        if not ok:
            ax.plot(
                ecut[i], abs_diff[i], "x",
                color="r", markersize=8, markeredgewidth=2
            )

    # =====================
    # Horizontal reference lines
    # =====================
    if y_horizon_lines:
        for y0 in y_horizon_lines:
            ax.axhline(y=y0, linestyle="--", linewidth=1, color="k")

    # =====================
    # Axes & labels
    # =====================
    ax.set_xlabel("Ecut (Ha)")
    ax.set_ylabel(rf"$|\Delta \delta|$ ({unit})")
    ax.grid(True)

    if y_range:
        ax.set_ylim(y_range)

    ymax = ax.get_ylim()[1]

    # =====================
    # Vertical convergence lines
    # =====================
    colors = {
        "low": "C1",
        "normal": "C2",
        "high": "C3",
    }

    for label, ec in conv_ecut.items():
        if ec is None:
            continue

        ax.axvline(
            ec,
            linestyle="--",
            linewidth=1.5,
            color=colors[label],
            label=f"{label} conv",
        )

        ax.text(
            ec,
            ymax * 0.98,
            f"{label} {ec} Ha",
            rotation=90,
            va="top",
            ha="right",
            fontsize=9,
            color=colors[label],
            clip_on=True,
        )

    title = f"{element} - {config} (Normalised Δ$_1$ difference convergence)"
    ax.set_title(title)
    ax.legend()

    plt.tight_layout()
    plt.savefig(f"{element}-{config}-abs-delta-convergence.pdf")
    plt.show()

    return fig


if __name__ == "__main__":
    plot_abs_delta_vs_ecut(
        "eos_converge_results.json",
        y_range=[0, 4]
    )
