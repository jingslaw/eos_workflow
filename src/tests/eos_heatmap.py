import matplotlib.pyplot as plt
import numpy as np
import os
import json

from eos_workflow.utilities import LANTHANIDE_ELEMENTS, ALL_ELEMENTS
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

# ============================================================
# Base colormap (unchanged)
# ============================================================
cmap = LinearSegmentedColormap.from_list(
    "ref_cmap",
    [
        # (0.0, "#0203BF"),
        # (0.33 / 1.65, "#ffff66"),
        # (1.0, "#cc0000"),

        (0.0, "#555998"),
        (0.1, "#6B71AD"),
        (0.33, "#EEE992"),
        (1.00, "#F28C53"),

        # (0.00, "#006837"),  # dark green
        # (0.05, "#1a9850"),  # green
        # (0.10, "#66bd63"),  # light green

        # (0.17, "#d9ef8b"),  # yellow-green
        # (0.27, "#fee08b"),  # yellow-orange
        # (0.33, "#fdae61"),  # light orange

        # (0.45, "#f46d43"),  # strong orange
        # (0.65, "#e34a33"),  # orange-red
        # (1.00, "#b30000"),  # dark red
    ],
)

def truncate_cmap(cmap_in, minval=0.0, maxval=1.0, n=256, name="trunc"):
    """Return a truncated version of a colormap."""
    colors = cmap_in(np.linspace(minval, maxval, n))
    return LinearSegmentedColormap.from_list(name, colors)

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
    "text.latex.preamble": r"\usepackage{sfmath}",
})

configuration = ["BCC", "FCC", "SC", "Diamond", "X2O", "XO", "X2O3", "XO2", "X2O5", "XO3"]
config_labels = ["BCC", "FCC", "SC", "Diamond", "X$_2$O", "XO",
                 "X$_2$O$_3$", "XO$_2$", "X$_2$O$_5$", "XO$_3$"]

file_path = "/home/wjing/project/eos/Dojov1.0"

# ============================================================
# Load data
# ============================================================
eos_results = []
for i in range(57, 72):
    element = ALL_ELEMENTS[i - 1]
    json_path = os.path.join(file_path, f"{i}_{element}.json")

    with open(json_path, "r") as fp:
        result = json.load(fp)

    eos_result = result[element]
    tmp = [eos_result[c]["rel_errors_vec_length"] for c in configuration]

    unaries = sum(tmp[:4]) / 4
    oxides = sum(tmp[-6:]) / 6
    together = sum(tmp) / 10

    tmp += [unaries, oxides, together]
    eos_results.append(tmp)

eos_results = np.array(eos_results)

config_labels += [r"$\bar{\nu}_{\rm u}$", r"$\bar{\nu}_{\rm o}$", r"$\bar{\nu}_{10}$"]

# ============================================================
# Plot
# ============================================================
fig, ax = plt.subplots(figsize=(6.0, 6.5), dpi=150, constrained_layout=True)

# --- Keep GRID mapping unchanged (0 -> 1.65)
# vmax_phys = 1.65
#vmax_phys = 1.00
#norm_img = Normalize(vmin=0.0, vmax=vmax_phys, clip=True)
im = ax.imshow(eos_results, cmap=cmap)

# ============================================================
# Colorbar that is 0.00 -> 1.00, with "over" values shown in red
# ============================================================
# We want 1.00 on the bar to correspond to the SAME color as data=1.00 in the grid.
# In the grid, data=1.00 corresponds to colormap position 1.00/1.65.
# cut = 1.0 / vmax_phys  # ~0.606
cut = 1.0

# Use only the [0, cut] portion of the original colormap so 1.00 becomes orange
cmap_cb = truncate_cmap(cmap, 0.0, cut, name="ref_cmap_upto_1")
cmap_cb.set_over("#cc0000")  # red cap for values > 1.00

norm_cb = Normalize(vmin=0.0, vmax=1.0, clip=False)
sm = ScalarMappable(norm=norm_cb, cmap=cmap_cb)
sm.set_array([])

cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
cbar.set_ticks([0.00, 0.10, 0.33, 1.00])
cbar.set_ticklabels(["0.00", "0.10", "0.33", "1.00"])
cbar.ax.minorticks_off()

for label in cbar.ax.get_yticklabels():
    label.set_fontfamily("sans-serif")
    label.set_fontname("Helvetica")

# ============================================================
# Axes ticks and labels
# ============================================================
ax.set_xticks(
    range(len(config_labels)),
    labels=config_labels,
    rotation=45,
    ha="right",
    rotation_mode="anchor",
    fontsize=14,
)
ax.set_yticks(range(len(LANTHANIDE_ELEMENTS)), labels=LANTHANIDE_ELEMENTS, fontsize=14)

# Cell annotations
for i in range(len(LANTHANIDE_ELEMENTS)):
    for j in range(len(config_labels)):
        nu = eos_results[i, j]
        ax.text(
            j, i, f"{nu:.2f}",
            ha="center", va="center",
            color="white" if nu < 0.1 else "black",
        )

ax.set_title(r"EOS discrepancy: $\nu$", fontsize=14)

plt.savefig("eos-nu.pdf", bbox_inches="tight", pad_inches=0.01)
plt.show()

