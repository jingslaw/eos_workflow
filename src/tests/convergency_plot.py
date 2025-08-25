import json
import matplotlib.pyplot as plt


def plot_phonon_vs_ecut(filename, acoustic_only=True):
    # Load data
    with open(filename, "r") as f:
        data = json.load(f)

    qpts = data["qpt_list"]
    ecut_keys = [k for k in data.keys() if k.startswith("ecut-")]

    # Map q-points to collected data
    qpt_results = {tuple(q): {"ecut": [], "freqs": [], "state": []} for q in qpts}

    for ecut_key in sorted(ecut_keys, key=lambda x: int(x.split("-")[1])):
        ecut = int(ecut_key.split("-")[1])
        for entry in data[ecut_key]:
            qpt = tuple(entry["phonon wavevector"])
            freqs = entry["phonon frequencies (cm^-1)"]
            state = entry["calculation state"]

            qpt_results[qpt]["ecut"].append(ecut)
            qpt_results[qpt]["freqs"].append(freqs)
            qpt_results[qpt]["state"].append(state)

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
            num_modes = min(3, num_modes)  # only take 3 acoustic branches

        for mode in range(num_modes):
            y = [f[mode] for f in freqs]

            # Line + circle markers for all points
            ax.plot(
                ecut, y, "-o", color=f"C{mode}", label=f"Acoustic {mode+1}"
            )

            # Overlay red X at unconverged points
            for i in range(len(ecut)):
                if states[i] != "successful":
                    ax.plot(ecut[i], y[i], "rx", markersize=8, markeredgewidth=2)

        ax.set_title(f"q-point {qpt}")
        ax.set_ylabel("Frequency (cm$^{-1}$)")
        ax.grid(True)

        # Build legend: one per branch + one for unconverged
        handles, labels = ax.get_legend_handles_labels()
        # Add a dummy handle for unconverged marker if not already added
        if "Unconverged" not in labels:
            handles.append(plt.Line2D([0], [0], color="r", marker="x",
                                      linestyle="None", markersize=8,
                                      markeredgewidth=2))
            labels.append("Unconverged")
        ax.legend(handles, labels)

    axes[-1].set_xlabel("Ecut (Ha)")

    # === Add global label with element + configuration ===
    element = data.get("element", "")
    config = data.get("configuration", "")
    plt.suptitle(f"{element} - {config}", fontsize=14, fontweight="bold")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"{element}-{config}-convergency.pdf")


if __name__ == "__main__":
    # For Diamond only acoustic modes
    plot_phonon_vs_ecut("phonon.txt", acoustic_only=True)
