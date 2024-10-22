import json
import click
import numpy as np
import matplotlib.pyplot as plt
from math import ceil
from eos_workflow.delta_metric import birch_murnaghan_function, load_ae_birch_murnaghan


# referred from aiida_sssp_workflow.cli.inspect and considered birch-murnaghan fitting failed case of PS
def eos_plot(
    ax, ref_ae_data, energy0, volumes, energies, line_data=None, nu=None, title="EOS", fontsize=8
):
    """plot EOS result on ax"""
    dense_volume_max = max(volumes)
    dense_volume_min = min(volumes)

    dense_volumes = np.linspace(dense_volume_min, dense_volume_max, 100)

    ref_v0, ref_b0, ref_b1 = ref_ae_data
    ae_eos_fit_energy = birch_murnaghan_function(
        v=dense_volumes,
        e0=energy0,  # in future update E0 from referece json, where ACWF has E0 stored.
        v0=ref_v0,
        b0=ref_b0,
        b1=ref_b1,
    )
    # Plot AE EOS & raw volume-energy points
    ax.tick_params(axis="y", labelsize=6, rotation=45)
    ax.plot(volumes, energies, "ob", label="RAW equation of state")
    ax.plot(dense_volumes, ae_eos_fit_energy, "--r", label="AE")

    if line_data:
        v0, b0, b1 = line_data
        psp_eos_fit_energy = birch_murnaghan_function(
            v=dense_volumes,
            e0=energy0,
            v0=v0,
            b0=b0,
            b1=b1,
        )
        ax.axvline(v0, linestyle="--", color="gray")
        ax.plot(dense_volumes, psp_eos_fit_energy, "-b", label=f"Pseudo")
        ax.fill_between(
            dense_volumes,
            ae_eos_fit_energy,
            psp_eos_fit_energy,
            alpha=0.5,
            color="red",
        )

    if nu:
        center_x = (max(volumes) + min(volumes)) / 2
        center_y = (max(energies) + min(energies)) / 2

        # write text of nu value in close middle
        nu = round(nu, 3)
        ax.text(center_x, center_y, f"$\\nu$={nu}")

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.0), fontsize=fontsize)
    ax.set_xlabel("Cell volume per formula unit ($\\AA^3$)", fontsize=fontsize)
    ax.set_ylabel("$E - TS$ (eV)", fontsize=fontsize)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=fontsize)


def convergence_plot(
    ax,
    convergence_data,
    converged_xy,
    y_thresholds_range,
    y_label,
    title,
):
    xs = convergence_data["xs"]
    ys = convergence_data["ys"]

    # xlim a bit larger than the range of xs
    xlims = (min(xs) * 0.9, max(xs) * 1.02)

    ax.plot(xs, ys, "-x")
    ax.scatter(
        *converged_xy,
        marker="s",
        s=100,
        label=f"Converge at {converged_xy[0]} Ha",
        facecolors="none",
        edgecolors="red",
    )
    ax.fill_between(
        x=xlims,
        y1=y_thresholds_range[0],
        y2=y_thresholds_range[1],
        alpha=0.3,
        color="green",
    )
    ax.legend(loc="upper right", fontsize=8)

    # twice the range of ylimits if the y_thresholds_range just cover the y range
    if y_thresholds_range[1] > max(ys) or y_thresholds_range[1] < max(ys) * 4:
        y_max = y_thresholds_range[1] * 2
        y_min = -0.05 * y_max
        ax.set_ylim(bottom=y_min, top=y_max)

    # change ticks size
    ax.tick_params(axis="both", labelsize=6)

    ax.set_xlim(*xlims)
    ax.set_ylabel(y_label, fontsize=8)
    ax.set_xlabel("Cutoff (Ry)", fontsize=8)
    ax.set_title(title, fontsize=8)


def eos_inspect(filepath="eos_fitting_results.json"):
    with open(filepath, 'r') as fp:
        eos_results = json.load(fp)
    for element, results in eos_results.items():
        rows = ceil(len(results) / 2)
        fig, axs = plt.subplots(rows, 2, figsize=(8.27, 11.69), dpi=100)
        i = 0
        for config, res in results.items():
            volumes = res["volumes"]
            energies = res["energies"]
            title = f"EOS of {config}"
            ax = axs.flat[i]
            if res["eos_is_converged"] is True:
                ref_data = res["reference_ae_V0_B0_B1"]
                data = res["birch_murnaghan_results"]
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
                data = None
                nu = None
            eos_plot(ax, ref_data, energy0, volumes, energies, data, nu, title=title)
            i += 1
        # fig to pdf
        fig.tight_layout()
        fig.savefig(f"{element}_precision.pdf", bbox_inches="tight")


def converge_inspect(filepath="eos_converge_results.json"):
    with open(filepath, 'r') as fp:
        eos_results = json.load(fp)
    fig, ax = plt.subplots(1, 1, figsize=(11.69, 8.27), dpi=100)
    y_thresholds_range = [0.0, 0.2]
    y_label = 'Δ meV/atom'
    title = "eos convergence wrt wavefunction cutoff"

    element = eos_results["element"]
    configuration = eos_results["configuration"]
    ecut_list = eos_results['ecut']
    delta_list = eos_results['delta/natoms']
    conv_data = {"xs": [], "ys": []}
    for i in range(len(ecut_list)):
        if delta_list[i] != "NaN":
            conv_data["xs"].append(ecut_list[i])
            conv_data["ys"].append(delta_list[i])
    wfc_scan_healthy = len(conv_data["xs"]) / len(ecut_list)

    if wfc_scan_healthy != 1:
        color = "red"
    else:
        color = "green"
    click.secho(
        f"Convergence scan healthy check for {property}: "
        f"wavefunction scan = {round(wfc_scan_healthy * 100, 2)}%", fg=color,
    )
    if len(conv_data["xs"]) == 0:
        click.secho("All eos calculations are failed", fg="red")
        ax.text(0.5, 0.5, "All eos calculations are failed", color="red")
        ax.set_axis_off()
    else:
        reference = conv_data["ys"][-1]
        conv_data["ys"] = [x - reference for x in conv_data["ys"]]
        index = 0
        for delta in reversed(conv_data["ys"]):
            if abs(delta) > y_thresholds_range[1]:
                index = conv_data["ys"].index(delta)
                break
        if abs(conv_data["ys"][0]) < y_thresholds_range[1]:
            pass
        else:
            index += 1
        _x = conv_data["xs"][index]
        _y = conv_data["ys"][index]
        converged_xy = (_x, _y)
        convergence_plot(ax, conv_data, converged_xy, y_thresholds_range, y_label=y_label, title=title)
    # fig to pdf
    fig.tight_layout()
    fig.savefig(f"{element}-{configuration}-converge.pdf", bbox_inches="tight")
