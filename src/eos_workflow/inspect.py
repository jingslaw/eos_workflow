import json
import numpy as np
import matplotlib.pyplot as plt
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


def eos_inspect(filepath="eos_fitting_results.json"):
    with open(filepath, 'r') as fp:
        eos_results = json.load(fp)
    for element, results in eos_results.items():
        rows = len(results) // 2
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
