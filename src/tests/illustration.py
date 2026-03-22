import json
from math import ceil
import numpy as np
import matplotlib.pyplot as plt
# from eos_workflow.eos_inspect import eos_plot
from eos_workflow.delta_metric import load_ae_birch_murnaghan, birch_murnaghan_function
from eos_workflow.delta_metric import birch_murnaghan_fit, metric_analyze


def eos_plot(
    ax, ref_ae_data, energy0, volumes, energies, line_data=None, nu=None, title="EOS", fontsize=17, position=None
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
        if v0 > dense_volume_max:
            ax.axvline(dense_volume_max, linestyle="--", color="gray")
        elif v0 < dense_volume_min:
            ax.axvline(dense_volume_min, linestyle="--", color="gray")
        else:
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
        nu = round(nu, 3)
        ax.text(
            0.65, 0.5,
            f"$\\nu$={nu}",
            transform=ax.transAxes,
            fontsize=fontsize,
        )

    if position == "B":
        ax.legend(loc="upper center", bbox_to_anchor=(0.3, 0.5), fontsize=fontsize)
    elif position == "C":
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.0), fontsize=fontsize)
    else:
        ax.legend(loc="upper center", bbox_to_anchor=(0.25, 1.0), fontsize=fontsize)

    ax.set_xlabel("Cell volume per formula unit ($\\AA^3$)", fontsize=fontsize)
    ax.set_ylabel("$E - TS$ (eV)", fontsize=fontsize)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=fontsize)


def eos_inspect(filepath="eos_fitting_results.json"):
    LABEL_POSSISION = ['T', 'B', 'C']
    with open(filepath, 'r') as fp:
        eos_results = json.load(fp)
    for element, results in eos_results.items():
        # rows = ceil(len(results) / 2)
        rows = len(results)
        fig, axs = plt.subplots(rows, 1, figsize=(8.27, 11.69), dpi=100)
        i = 0
        for config, res in results.items():
            config = config.split('-')[0]
            volumes = res["volumes"]
            energies = res["energies"]
            title = f"EOS of {config}"
            # TODO: now this can only handle config number > 1
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
                fitting_results = birch_murnaghan_fit(res)

                v0 = fitting_results["volume0"]
                b0 = fitting_results["bulk_modulus0"]
                b1 = fitting_results["bulk_deriv0"]
                num = res["num_of_atoms"]

                eos_results = metric_analyze(element, config, v0, b0, b1, num)
                data = eos_results["birch_murnaghan_results"]
                nu = eos_results["rel_errors_vec_length"]
            eos_plot(ax, ref_data, energy0, volumes, energies, data, nu, title=title, position=LABEL_POSSISION[i])
            i += 1
        # fig to pdf
        fig.tight_layout()
        fig.savefig(f"{element}_precision.pdf", bbox_inches="tight")
        return fig


if __name__ == "__main__":
    eos_inspect()
