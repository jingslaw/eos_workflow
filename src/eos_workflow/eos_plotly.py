import json
import numpy as np
import csv
from math import ceil

import plotly.graph_objects as go
from plotly.subplots import make_subplots
from eos_workflow.delta_metric import birch_murnaghan_function, load_ae_birch_murnaghan


# ---------------------------
# Plot function (replacement of eos_plot)
# ---------------------------
def eos_plot_plotly(
    fig,
    row,
    col,
    ref_ae_data,
    energy0,
    volumes,
    energies,
    line_data=None,
    nu=None,
):
    dense_volume_max = max(volumes)
    dense_volume_min = min(volumes)
    dense_volumes = np.linspace(dense_volume_min, dense_volume_max, 100)

    ref_v0, ref_b0, ref_b1 = ref_ae_data

    ae_eos_fit_energy = birch_murnaghan_function(
        v=dense_volumes,
        e0=energy0,
        v0=ref_v0,
        b0=ref_b0,
        b1=ref_b1,
    )

    show_legend = (row == 1 and col == 1)

    # ---- RAW data points ----
    fig.add_trace(
        go.Scatter(
            x=volumes,
            y=energies,
            mode="markers",
            name="RAW equation of state",
            marker=dict(size=10, color="blue"),
            showlegend=show_legend,
        ),
        row=row,
        col=col,
    )

    # ---- AE curve ----
    fig.add_trace(
        go.Scatter(
            x=dense_volumes,
            y=ae_eos_fit_energy,
            mode="lines",
            name="AE",
            line=dict(dash="dash", color="red"),
            showlegend=show_legend,
        ),
        row=row,
        col=col,
    )

    # ---- Pseudo curve + filled area ----
    if line_data:
        v0, b0, b1 = line_data

        psp_eos_fit_energy = birch_murnaghan_function(
            v=dense_volumes,
            e0=energy0,
            v0=v0,
            b0=b0,
            b1=b1,
        )

        fig.add_trace(
            go.Scatter(
                x=dense_volumes,
                y=psp_eos_fit_energy,
                mode="lines",
                name="Pseudo",
                line=dict(color="blue"),
                showlegend=show_legend,
            ),
            row=row,
            col=col,
        )

        # fill between AE and Pseudo
        fig.add_trace(
            go.Scatter(
                x=np.concatenate([dense_volumes, dense_volumes[::-1]]),
                y=np.concatenate([ae_eos_fit_energy, psp_eos_fit_energy[::-1]]),
                fill="toself",
                fillcolor="rgba(255,0,0,0.3)",
                line=dict(width=0),
                showlegend=False,
                name="Difference",
            ),
            row=row,
            col=col,
        )

        # vertical line at v0
        fig.add_vline(
            x=v0,
            line=dict(color="gray", dash="dash"),
            row=row,
            col=col,
        )

    # ---- ν annotation ----
    if nu:
        x_min, x_max = min(volumes), max(volumes)
        y_min, y_max = min(energies), max(energies)

        # Place annotation at 3/4 of x-range, 90% of y-range
        ann_x = x_min + 0.75 * (x_max - x_min)
        ann_y = y_min + 0.9 * (y_max - y_min)

        fig.add_annotation(
            x=ann_x,
            y=ann_y,
            text=f"ν = {round(nu,3)}",
            showarrow=False,
            row=row,
            col=col,
        )

    # ---- Axis formatting ----
    fig.update_xaxes(
        title_text="Cell volume per formula unit (Å³)",
        row=row,
        col=col,
        showticklabels=False,
    )

    fig.update_yaxes(
        title_text="E - TS (eV)",
        row=row,
        col=col,
        showticklabels=False,
    )


# ---------------------------
# Main function (replacement of eos_inspect)
# ---------------------------
def eos_inspect_plotly(filepath="eos_fitting_results.json", save_html=True):
    with open(filepath, 'r') as fp:
        eos_results = json.load(fp)

    figs = {}

    for element, results in eos_results.items():

        nplots = len(results)
        rows = ceil(nplots / 2)
        cols = 2

        fig = make_subplots(
            rows=rows,
            cols=cols,
            subplot_titles=[f"EOS of {cfg}" for cfg in results.keys()],
            vertical_spacing=0.02,
            horizontal_spacing=0.02,
        )

        i = 0
        for config, res in results.items():
            row = i // 2 + 1
            col = i % 2 + 1

            volumes = res["volumes"]
            energies = res["energies"]

            if res["eos_is_converged"] is True:
                ref_data = res["reference_ae_V0_B0_B1"]
                line_data = res["birch_murnaghan_results"]
                energy0 = res["energy0"]
                nu = res["rel_errors_vec_length"]
            else:
                bm_fit = load_ae_birch_murnaghan(element, config)
                ref_data = [
                    bm_fit["min_volume"],
                    bm_fit["bulk_modulus_ev_ang3"],
                    bm_fit["bulk_deriv"],
                ]
                energy0 = min(energies)
                line_data = None
                nu = None

            eos_plot_plotly(
                fig,
                row,
                col,
                ref_data,
                energy0,
                volumes,
                energies,
                line_data=line_data,
                nu=nu,
            )

            i += 1

        fig.update_layout(
            # title=f"{element} EOS precision check",
            height=350 * rows,
            width=900,
            showlegend=True,
        )

        if save_html:
            fig.write_html(f"{element}_precision.html")

        figs[element] = fig

    return figs


def plot_phonon_vs_ecut_plotly(
    filename,
    acoustic_only=True,
    y_ranges=None,
    error_type=None,
    csv_output=True,
    y_horizon_lines=None,
    abs_val=True,
    outfile_html=True,
    outfile_image=False,
):
    """
    Plot phonon frequencies vs cutoff energy using Plotly and optionally export errors as CSV.
    """

    # ===== Load data =====
    with open(filename, "r") as f:
        data = json.load(f)

    qpts = data["qpt_list"]
    ecut_keys = sorted(
        [k for k in data.keys() if k.startswith("ecut-")],
        key=lambda x: int(x.split("-")[1])
    )

    # Map q-points to collected data
    qpt_results = {tuple(q): {"ecut": [], "freqs": [], "state": []} for q in qpts}

    for ecut_key in ecut_keys:
        ecut = int(ecut_key.split("-")[1])
        for entry in data[ecut_key]:
            qpt = tuple(entry["phonon wavevector"])
            freqs = entry["phonon frequencies (cm^-1)"]
            state = entry["calculation state"]

            qpt_results[qpt]["ecut"].append(ecut)
            qpt_results[qpt]["freqs"].append(freqs)
            qpt_results[qpt]["state"].append(state)

    # Reference cutoff
    ref_key = max(ecut_keys, key=lambda x: int(x.split("-")[1]))
    ref_freqs = {
        tuple(entry["phonon wavevector"]): entry["phonon frequencies (cm^-1)"]
        for entry in data[ref_key]
    }

    # ===== Prepare plotly subplots =====
    fig = make_subplots(
        rows=len(qpts),
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.07,
        subplot_titles=[f"q-point {tuple(q)}" for q in qpts]
    )

    csv_rows = []

    for row_idx, qpt in enumerate(qpts, start=1):

        result = qpt_results[tuple(qpt)]
        ecut = result["ecut"]
        freqs = result["freqs"]
        states = result["state"]

        num_modes = len(freqs[0])
        if acoustic_only:
            num_modes = min(3, num_modes)

        for mode in range(num_modes):

            y = [f[mode] for f in freqs]
            ref = ref_freqs[tuple(qpt)][mode]

            if abs_val:
                y = [abs(i) for i in y]
                ref = abs(ref)

            # Errors
            errors, abs_errors, rel_errors = [], [], []

            for val in y:
                delta = val - ref
                errors.append(delta)
                abs_errors.append(abs(delta))
                rel_errors.append(np.nan if abs(ref) < 1e-6 else abs(delta) / abs(ref))

            # Save CSV
            for i, ec in enumerate(ecut):
                csv_rows.append({
                    "qpt": qpt,
                    "mode": mode + 1,
                    "ecut": ec,
                    "freq": y[i],
                    "ref_freq": ref,
                    "error": errors[i],
                    "abs_error": abs_errors[i],
                    "rel_error": rel_errors[i],
                    "state": states[i],
                })

            # Decide what to plot
            if error_type == "error":
                y_plot = errors
                ylabel = "Error (cm⁻¹)"
            elif error_type == "abs_error":
                y_plot = abs_errors
                ylabel = "Absolute Error (cm⁻¹)"
            elif error_type == "rel_error":
                y_plot = rel_errors
                ylabel = "Relative Error"
            else:
                y_plot = y
                ylabel = "Frequency (cm⁻¹)"

            # ===== Line =====
            fig.add_trace(
                go.Scatter(
                    x=ecut,
                    y=y_plot,
                    mode="lines",
                    name=f"Acoustic {mode+1}",
                    legendgroup=f"{mode}",
                    showlegend=(row_idx == 1)
                ),
                row=row_idx, col=1
            )

            # ===== Markers =====
            for i in range(len(ecut)):

                if states[i] == "successful":
                    symbol = "circle"
                    color = f"rgba(0,0,0,0.7)"
                else:
                    symbol = "x"
                    color = "red"

                fig.add_trace(
                    go.Scatter(
                        x=[ecut[i]],
                        y=[y_plot[i]],
                        mode="markers",
                        marker=dict(symbol=symbol, size=8, color=color),
                        showlegend=False,
                        hovertemplate=(
                            f"Mode: {mode+1}<br>"
                            f"Ecut: {ecut[i]} Ha<br>"
                            f"Value: {y_plot[i]:.4f}<br>"
                            f"State: {states[i]}"
                        )
                    ),
                    row=row_idx, col=1
                )

        # ===== Horizontal line if requested =====
        if y_horizon_lines and tuple(qpt) in y_horizon_lines:
            y0 = y_horizon_lines[tuple(qpt)]
            fig.add_hline(
                y=y0,
                line_dash="dash",
                line_color="black",
                row=row_idx,
                col=1
            )

        # ===== y-range =====
        if y_ranges and tuple(qpt) in y_ranges:
            ymin, ymax = y_ranges[tuple(qpt)]
            fig.update_yaxes(range=[ymin, ymax], row=row_idx, col=1)

        fig.update_yaxes(title_text=ylabel, row=row_idx, col=1)

    fig.update_xaxes(title_text="Ecut (Ha)", row=len(qpts), col=1)

    element = data.get("element", "")
    config = data.get("configuration", "")
    title = f"{element} - {config}"
    if error_type:
        title += f" ({error_type})"

    fig.update_layout(
        title=dict(text=title, x=0.5),
        height=350 * len(qpts),
        width=800,
        template="plotly_white"
    )

    # ===== Output =====
    base = f"{element}-{config}-convergency-{error_type or 'freqs'}"

    if outfile_html:
        fig.write_html(base + ".html")

    if outfile_image:
        fig.write_image(base + ".pdf")  # needs: pip install kaleido

    fig.show()

    # ===== Save CSV =====
    if csv_output:
        csv_file = f"{element}-{config}-convergency-data.csv"
        with open(csv_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=[
                "qpt", "mode", "ecut", "freq", "ref_freq",
                "error", "abs_error", "rel_error", "state"
            ])
            writer.writeheader()
            for row in csv_rows:
                writer.writerow(row)

        print(f" Saved results to {csv_file}")

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
    filename = "phonon.txt"
    y_range = {(0.0, 0.0, 0.0): (0, 10), (0.5, 0.5, 0.5): (0, 10)}
    y_lines = {(0.0, 0.0, 0.0): 1.0, (0.5, 0.5, 0.5): 1.0}
    plot_phonon_vs_ecut_plotly(
        filename,
        acoustic_only=True,
        y_ranges=y_range,
        error_type="abs_error",
        y_horizon_lines=y_lines,
        abs_val=False,
    )
