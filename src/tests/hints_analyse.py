import json
import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt


def plot_energy_vs_ecut(
    filename,
    error_type='abs_error',
    y_range=None,
    y_horizon_lines=None,
):
    """
    Plot total energy convergence vs cutoff energy with automatic meV scaling.

    Parameters
    ----------
    filename : str
        Path to JSON data file (hints.txt).
    error_type : str or None
        None        : plot total energy directly (eV)
        'error'     : E(Ecut) - E(ref)   → plotted in meV
        'abs_error' : |E(Ecut) - E(ref)| → plotted in meV
        'rel_error' : |E(Ecut) - E(ref)| / |E(ref)|
    y_range : tuple or None
        (ymin, ymax) for y-axis limits (same unit as plotted data).
    y_horizon_lines : list[float] or None
        Horizontal reference lines (same unit as plotted data).
    """

    # =====================
    # Load data
    # =====================
    with open(filename, "r") as f:
        data = json.load(f)

    energy_unit = data.get("energy_unit", "eV")
    element = data.get("element", "")
    config = data.get("configuration", "")

    ecut_keys = sorted(
        [k for k in data if k.startswith("ecut-")],
        key=lambda x: int(x.split("-")[1]),
    )

    ecut = np.array([int(k.split("-")[1]) for k in ecut_keys])
    energy = np.array([data[k]["energy"] for k in ecut_keys])
    states = [data[k]["task_state"] for k in ecut_keys]

    # =====================
    # Reference energy
    # =====================
    ref_energy = energy[-1]
    delta = energy - ref_energy
    abs_delta = np.abs(delta)
    rel_delta = abs_delta / abs(ref_energy)

    # =====================
    # Decide what to plot
    # =====================
    scale = 1.0
    unit = energy_unit

    if error_type == "error":
        y = delta
        scale = 1e3   # eV → meV
        unit = "meV"
        ylabel = f"ΔE ({unit})"
    elif error_type == "abs_error":
        y = abs_delta
        scale = 1e3
        unit = "meV"
        ylabel = f"|ΔE| ({unit})"
    elif error_type == "rel_error":
        y = rel_delta
        ylabel = "Relative Error"
    else:
        y = energy
        ylabel = f"Total Energy ({energy_unit})"

    y = y * scale

    # =====================
    # Convergence thresholds (meV)
    # =====================
    thresholds_mev = {
        "low": 10.0,
        "normal": 5.0,
        "high": 2.0,
    }

    conv_ecut = {}
    for label, thr in thresholds_mev.items():
        idx = np.where(abs_delta * 1e3 < thr)[0]
        conv_ecut[label] = ecut[idx[0]] if len(idx) > 0 else None

    # =====================
    # Plot
    # =====================
    fig, ax = plt.subplots(figsize=(6, 4))

    ax.plot(ecut, y, "-o", label="Total Energy")

    for i, st in enumerate(states):
        if st != "successful":
            ax.plot(
                ecut[i], y[i], "x",
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
    ax.set_ylabel(ylabel)
    ax.grid(True)
    # =====================
    # Apply y-range FIRST
    # =====================
    if y_range:
        ax.set_ylim(y_range)

    ymax = ax.get_ylim()[1]

    # =====================
    # Vertical convergence lines (with legend labels)
    # =====================
    colors = {
        "low": "C1",
        "normal": "C2",
        "high": "C3",
    }

    labels_done = set()

    for label, ec in conv_ecut.items():
        if ec is None:
            continue

        legend_label = label if label not in labels_done else None
        labels_done.add(label)

        ax.axvline(
            ec,
            linestyle="--",
            linewidth=1.5,
            color=colors[label],
            label=legend_label,
        )

        # place text safely inside axis
        ax.text(
            ec,
            ymax * 0.98,
            f"{label} {conv_ecut[label]} Ha",
            rotation=90,
            va="top",
            ha="right",
            fontsize=9,
            color=colors[label],
            clip_on=True,
        )

    title = f"{element} - {config} (Ecut convergence)"
    if error_type:
        title += f" [{error_type}]"
    ax.set_title(title)

    ax.legend()

    plt.tight_layout()
    plt.savefig(
        f"{element}-{config}-energy-convergency-{error_type or 'energy'}.pdf"
    )
    plt.show()

    return fig


def plot_energy_vs_ecut_plotly(
    filename,
    error_type="abs_error",
    y_range=None,
    y_horizon_lines=None,
    html_out=None,
):
    """
    Plot total energy convergence vs cutoff energy using Plotly (HTML output).

    Parameters
    ----------
    filename : str
        Path to JSON data file (hints.txt).
    error_type : str or None
        None        : plot total energy directly (eV)
        'error'     : E(Ecut) - E(ref)   → meV
        'abs_error' : |E(Ecut) - E(ref)| → meV
        'rel_error' : |E(Ecut) - E(ref)| / |E(ref)|
    y_range : tuple or None
        (ymin, ymax) for y-axis limits.
    y_horizon_lines : list[float] or None
        Horizontal reference lines.
    html_out : str or None
        Output HTML filename.
    """

    # =====================
    # Load data
    # =====================
    with open(filename, "r") as f:
        data = json.load(f)

    energy_unit = data.get("energy_unit", "eV")
    element = data.get("element", "")
    config = data.get("configuration", "")

    ecut_keys = sorted(
        [k for k in data if k.startswith("ecut-")],
        key=lambda x: int(x.split("-")[1]),
    )

    ecut = np.array([int(k.split("-")[1]) for k in ecut_keys])
    energy = np.array([data[k]["energy"] for k in ecut_keys])
    states = [data[k]["task_state"] for k in ecut_keys]

    # =====================
    # Reference energy
    # =====================
    ref_energy = energy[-1]
    delta = energy - ref_energy
    abs_delta = np.abs(delta)
    rel_delta = abs_delta / abs(ref_energy)

    # =====================
    # Decide what to plot
    # =====================
    scale = 1.0
    unit = energy_unit

    if error_type == "error":
        y = delta
        scale = 1e3
        unit = "meV"
        ylabel = f"ΔE ({unit})"
    elif error_type == "abs_error":
        y = abs_delta
        scale = 1e3
        unit = "meV"
        ylabel = f"|ΔE| ({unit})"
    elif error_type == "rel_error":
        y = rel_delta
        ylabel = "Relative Error"
    else:
        y = energy
        ylabel = f"Total Energy ({energy_unit})"

    y = y * scale

    # =====================
    # Convergence thresholds (meV)
    # =====================
    thresholds_mev = {
        "low": 10.0,
        "normal": 5.0,
        "high": 2.0,
    }

    conv_ecut = {}
    for label, thr in thresholds_mev.items():
        idx = np.where(abs_delta * 1e3 < thr)[0]
        conv_ecut[label] = ecut[idx[0]] if len(idx) > 0 else None

    # =====================
    # Plotly figure
    # =====================
    fig = go.Figure()

    # main line
    fig.add_trace(
        go.Scatter(
            x=ecut,
            y=y,
            mode="lines+markers",
            name="Total Energy",
        )
    )

    # failed points
    failed_idx = [i for i, st in enumerate(states) if st != "successful"]
    if failed_idx:
        fig.add_trace(
            go.Scatter(
                x=ecut[failed_idx],
                y=y[failed_idx],
                mode="markers",
                marker=dict(symbol="x", size=10, color="red"),
                name="Failed",
            )
        )

    # =====================
    # Horizontal reference lines
    # =====================
    if y_horizon_lines:
        for y0 in y_horizon_lines:
            fig.add_hline(
                y=y0,
                line=dict(dash="dash", width=1),
            )

    # =====================
    # Vertical convergence lines
    # =====================
    colors = {
        "low": "orange",
        "normal": "green",
        "high": "purple",
    }

    ymax = y_range[1] if y_range else max(y)

    for label, ec in conv_ecut.items():
        if ec is None:
            continue

        fig.add_vline(
            x=ec,
            line=dict(
                dash="dash",
                width=2,
                color=colors[label],
            ),
            annotation_text="", # f"{label} {ec} Ha",
            # annotation_position="top right",
            annotation_font_color=colors[label],
        )

    for label, ec in conv_ecut.items():
        if ec is None:
            continue

        # vertical line as a trace (legend-aware)
        fig.add_trace(
            go.Scatter(
                x=[ec, ec],
                y=[min(y), max(y)],
                mode="lines",
                line=dict(
                    dash="dash",
                    width=2,
                    color=colors[label],
                ),
                name=f"{label} {ec} Ha",
                showlegend=True,
            )
        )


    # =====================
    # Layout
    # =====================
    title = f"{element} - {config} (Ecut convergence)"
    if error_type:
        title += f" [{error_type}]"

    fig.update_layout(
        title=title,
        xaxis_title="Ecut (Ha)",
        yaxis_title=ylabel,
        template="plotly_white",
        width=700,
        height=450,
    )

    if y_range:
        fig.update_yaxes(range=y_range)

    # =====================
    # Save HTML
    # =====================
    if html_out is None:
        html_out = f"{element}-{config}-energy-convergency-{error_type or 'energy'}.html"

    fig.write_html(html_out, include_plotlyjs="cdn")
    fig.show()

    return fig


if __name__ == "__main__":
    filename = "hints.txt"
    plot_energy_vs_ecut_plotly(
        filename,
        y_range=(0, 15),
    )
