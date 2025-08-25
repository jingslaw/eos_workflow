import json
import matplotlib.pyplot as plt
import numpy as np
import csv


def plot_phonon_vs_ecut(filename, acoustic_only=True, y_ranges=None, error_type=None, csv_output=True):
    """
    Plot phonon frequencies vs cutoff energy and optionally export errors as CSV.

    Parameters
    ----------
    filename : str
        Path to JSON data file.
    acoustic_only : bool
        Whether to plot only acoustic modes (default: True).
    y_ranges : dict[tuple, (ymin, ymax)]
        Optional dict mapping q-point tuples -> (ymin, ymax) for y-axis limits.
    error_type : str or None
        None: plot frequencies directly.
        'error': plot f(Ecut) - f(Ecut_ref).
        'abs_error': plot |f(Ecut) - f(Ecut_ref)|.
        'rel_error': plot relative error wrt f(Ecut_ref).
    csv_output : bool
        If True, write a CSV file with frequencies and errors.
    """

    # Load data
    with open(filename, "r") as f:
        data = json.load(f)

    qpts = data["qpt_list"]
    ecut_keys = [k for k in data.keys() if k.startswith("ecut-")]
    ecut_keys_sorted = sorted(ecut_keys, key=lambda x: int(x.split("-")[1]))

    # Map q-points to collected data
    qpt_results = {tuple(q): {"ecut": [], "freqs": [], "state": []} for q in qpts}

    for ecut_key in ecut_keys_sorted:
        ecut = int(ecut_key.split("-")[1])
        for entry in data[ecut_key]:
            qpt = tuple(entry["phonon wavevector"])
            freqs = entry["phonon frequencies (cm^-1)"]
            state = entry["calculation state"]

            qpt_results[qpt]["ecut"].append(ecut)
            qpt_results[qpt]["freqs"].append(freqs)
            qpt_results[qpt]["state"].append(state)

    # Reference cutoff (largest cutoff available, e.g. 150 Ha)
    ref_key = max(ecut_keys_sorted, key=lambda x: int(x.split("-")[1]))
    ref_freqs = {}
    for entry in data[ref_key]:
        qpt = tuple(entry["phonon wavevector"])
        ref_freqs[qpt] = entry["phonon frequencies (cm^-1)"]

    # Prepare CSV rows
    csv_rows = []

    # Plot
    fig, axes = plt.subplots(len(qpts), 1, figsize=(6, 4*len(qpts)), sharex=True)
    if len(qpts) == 1:
        axes = [axes]

    for ax, qpt in zip(axes, qpts):
        ecut = qpt_results[tuple(qpt)]["ecut"]
        freqs = qpt_results[tuple(qpt)]["freqs"]
        states = qpt_results[tuple(qpt)]["state"]

        num_modes = len(freqs[0])
        if acoustic_only:
            num_modes = min(3, num_modes)  # only 3 acoustic branches

        for mode in range(num_modes):
            y = [f[mode] for f in freqs]
            ref = ref_freqs[tuple(qpt)][mode]

            # Compute errors for CSV
            errors = []
            abs_errors = []
            rel_errors = []
            for val in y:
                delta = val - ref
                errors.append(delta)
                abs_errors.append(abs(delta))
                rel_errors.append(np.nan if abs(ref) < 1e-6 else abs(delta) / abs(ref))

            # Write CSV rows
            for i, ec in enumerate(ecut):
                csv_rows.append({
                    "qpt": qpt,
                    "mode": mode+1,
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
                ylabel = "Error (cm$^{-1}$)"
            elif error_type == "abs_error":
                y_plot = abs_errors
                ylabel = "Absolute Error (cm$^{-1}$)"
            elif error_type == "rel_error":
                y_plot = rel_errors
                ylabel = "Relative Error"
            else:
                y_plot = y
                ylabel = "Frequency (cm$^{-1}$)"

            # Draw connecting line
            ax.plot(ecut, y_plot, "-", color=f"C{mode}", label=f"Acoustic {mode+1}")

            # Overlay markers
            for i in range(len(ecut)):
                if states[i] == "successful":
                    ax.plot(ecut[i], y_plot[i], "o", color=f"C{mode}")
                else:
                    ax.plot(ecut[i], y_plot[i], "x", color="r", markersize=8, markeredgewidth=2)

        ax.set_title(f"q-point {qpt}")
        ax.set_ylabel(ylabel)
        ax.grid(True)

        # Apply user-provided y-limits
        if y_ranges and tuple(qpt) in y_ranges:
            ax.set_ylim(y_ranges[tuple(qpt)])

        # Legend
        handles, labels = ax.get_legend_handles_labels()
        if "Unconverged" not in labels:
            handles.append(plt.Line2D([0], [0], color="r", marker="x",
                                      linestyle="None", markersize=8,
                                      markeredgewidth=2))
            labels.append("Unconverged")
        ax.legend(handles, labels)

    axes[-1].set_xlabel("Ecut (Ha)")

    element = data.get("element", "")
    config = data.get("configuration", "")
    title = f"{element} - {config}"
    if error_type:
        title += f" ({error_type})"
    plt.suptitle(title, fontsize=14, fontweight="bold")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"{element}-{config}-convergency-{error_type or 'freqs'}.pdf")
    plt.show()

    # === Save CSV ===
    if csv_output:
        csv_file = f"{element}-{config}-convergency-data.csv"
        with open(csv_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=[
                "qpt", "mode", "ecut", "freq", "ref_freq", "error", "abs_error", "rel_error", "state"
            ])
            writer.writeheader()
            for row in csv_rows:
                writer.writerow(row)
        print(f"Saved results to {csv_file}")


if __name__ == "__main__":
    # For Diamond only acoustic modes
    plot_phonon_vs_ecut("phonon.txt", acoustic_only=True)
