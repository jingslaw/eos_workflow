import json
import numpy as np
from scipy.integrate import quad
from pathlib import Path

LANTHANIDE_ELEMENTS = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
OXIDE_CONFIGURATIONS = ["XO", "XO2", "XO3", "X2O", "X2O3", "X2O5"]
UNARIE_CONFIGURATIONS = ["BCC", "FCC", "SC", "Diamond"]


def birch_murnaghan_fit(volume_energy):
    # energy (eV) -- volume (A^3)
    volumes = np.array(volume_energy["volumes"])
    energies = np.array(volume_energy["energies"])
    num_of_atoms = volume_energy["num_of_atoms"]
    fitdata = np.polyfit(volumes ** (-2.0 / 3.0), energies, 3, full=True)

    ssr = fitdata[1]
    sst = np.sum((energies - np.average(energies)) ** 2.0)
    if len(ssr) == 0:
        residuals0 = 1000
    else:
        residuals0 = ssr[0] / sst
    # adjusted_residuals0 = (1 - residuals0) * (num_of_atoms - 1) / (num_of_atoms - 4 - 1)

    deriv0 = np.poly1d(fitdata[0])
    deriv1 = np.polyder(deriv0, 1)
    deriv2 = np.polyder(deriv1, 1)
    deriv3 = np.polyder(deriv2, 1)

    volume0 = 0
    x = 0
    for x in np.roots(deriv1):
        if x > 0 and deriv2(x) > 0:
            volume0 = x ** (-3.0 / 2.0)
            break

    if volume0 == 0:
        print("no positive volume0 found for the birch murnaghan fitting")

    if not isinstance(volume0, float):
        print(f"get a complex volume0={volume0}, birch murnaghan fitting failed")
        volume0 = 0

    if volume0 != 0:
        # using equations below to obtain energy0 (E0), bulk_modulus0 (B0), and bulk_deriv0 (B1)
        # x0 = V0 ** (-2.0/3.0)
        # deriv0(x0) = E0
        # deriv1(x0) = 0
        # deriv2(x0) = 9.0 / 4.0 * B0 * V0 ** (7.0 / 3.0)
        # deriv3(x0) = 27.0 / 8.0 * B0 * V0 ** 3 * (B1 - 4)
        # So we have:
        bulk_modulus0 = 4.0 / 9.0 * deriv2(x) * x ** (7.0 / 2.0)
        bulk_deriv0 = 4.0 + 2.0 / 3.0 * x * deriv3(x) / deriv2(x)
        energy0 = deriv0(x)
    else:
        bulk_modulus0 = -1.0
        bulk_deriv0 = -1.0
        energy0 = -1.0
        residuals0 = 10000

    echarge = 1.60217733e-19
    return {
            "volume0": np.round(volume0, 7),
            "energy0": np.round(energy0, 7),
            "num_of_atoms": int(num_of_atoms),
            "bulk_modulus0": np.round(bulk_modulus0, 7),
            "bulk_modulus0_GPa": np.round(bulk_modulus0 * echarge * 1.0e21, 7),
            "bulk_deriv0": np.round(bulk_deriv0, 7),
            "residuals0": np.round(residuals0, 7),
            "volume0_unit": "A^3",
            "bulk_modulus0_unit": "eV A^3",
            "bulk_modulus0_GPa_unit": "GPa",
    }


def eos_coefficients(v0, b0, b1):
    # EOS, Total energy E is the function of Volume V according to Birch-Murnaghan equation of state
    # if we extend E(V) as E(V) = E0 + a3 * x**3 + a2 * x**2 + a1 * x + a0, where x = V ** (-2/3)
    # coefficients a3, a2, a1, a0 can be represented by v0, b0, and b1 as:

    a3 = 9.0 * v0 ** 3.0 * b0 / 16.0 * (b1 - 4.0)
    a2 = 9.0 * v0 ** (7.0 / 3.0) * b0 / 16.0 * (14.0 - 3.0 * b1)
    a1 = 9.0 * v0 ** (5.0 / 3.0) * b0 / 16.0 * (3.0 * b1 - 16.0)
    a0 = 9.0 * v0 * b0 / 16.0 * (6.0 - b1)
    return [a3, a2, a1, a0]


def calc_delta(v0f, b0f, b1f, v0w, b0w, b1w, useasymm=False):
    # delta is calculated by Equation 1 in reference DOI: 10.1126/science.aad3000

    vref = 30.0
    bref = 100.0 * 10.0 ** 9.0 / 1.602176565e-19 / 10.0 ** 30.0

    if useasymm:
        vi = 0.94 * v0w
        vf = 1.06 * v0w
    else:
        vi = 0.94 * (v0w + v0f) / 2.0
        vf = 1.06 * (v0w + v0f) / 2.0

    calc_eos = np.poly1d(eos_coefficients(v0f, b0f, b1f))
    referred_eos = np.poly1d(eos_coefficients(v0w, b0w, b1w))

    def integrand_delta(v):
        return (calc_eos(v ** (-2 / 3)) - referred_eos(v ** (-2 / 3))) ** 2

    def integrand_average_eos(v):
        return ((calc_eos(v ** (-2 / 3)) + referred_eos(v ** (-2 / 3))) / 2) ** 2

    delta = quad(integrand_delta, vi, vf)[0]
    delta = 1000 * np.sqrt(delta / (vf - vi))
    average_delta = quad(integrand_average_eos, vi, vf)[0]
    average_delta = 1000 * np.sqrt(average_delta / (vf - vi))
    deltarel = delta / average_delta * 100

    if useasymm:
        delta1 = delta / (v0w * b0w) * (vref * bref)
    else:
        delta1 = delta / (v0w + v0f) / (b0w + b0f) * 4.0 * (vref * bref)

    return delta, deltarel, delta1


def rel_errors_vec_length(v0w, b0w, b1w, v0f, b0f, b1f, prefact=100, weight_b0=1/20, weight_b1=1/400):
    """
    Returns the length of the vector formed by the relative error of V0, B0, B1
    THE SIGNATURE OF THIS FUNCTION HAS BEEN CHOSEN TO MATCH THE ONE OF ALL THE OTHER FUNCTIONS
    RETURNING A QUANTITY THAT IS USEFUL FOR COMPARISON, THIS SIMPLIFIES THE CODE LATER.
    Even though config_string is not usd
    """
    v0err = 2 * (v0w - v0f) / (v0w + v0f)
    b0err = 2 * (b0w - b0f) / (b0w + b0f)
    b1err = 2 * (b1w - b1f) / (b1w + b1f)
    leng = np.sqrt(v0err ** 2 + (weight_b0 * b0err) ** 2 + (weight_b1 * b1err) ** 2)
    return leng * prefact


def metric_analyze(element, configuration, v0, b0, b1, natoms):
    """
    The calcfunction calculate the metric factor.
    return delta factor with unit (eV/atom)

    The configuration can be one of:
    - RE: for Rare-earth element
    - GS: for BM fit results from Corttiner's paper
    - One of OXIDES:
        - XO
        - XO2
        - XO3
        - X2O
        - X2O3
        - X2O5
    - One of UNARIES:
        - BCC
        - FCC
        - SC
        - Diamond
    - For actinides and Ar, Fr, Ra: using FCC from ACWF dataset

    conf_key is key in json file for configurations of every element.
    """
    if configuration == "RE":
        assert element in LANTHANIDE_ELEMENTS

        ref_json = "WIEN2K_LANN.json"
        conf_key = f"{element}N"
    elif configuration == "GS":
        ref_json = "WIEN2K_GS.json"
        conf_key = f"{element}"
    elif configuration in UNARIE_CONFIGURATIONS:
        ref_json = "AE-average-unaries.json"
        conf_key = f"{element}-X/{configuration}"
    elif configuration in OXIDE_CONFIGURATIONS:
        ref_json = "AE-average-oxides.json"
        conf_key = f"{element}-{configuration}"
    else:
        print(f"{configuration} is not one of RE, GS, BCC, FCC, SC, Diamond, X2O, XO, X2O3, XO2, X2O5, XO3")
        exit(-1)

    package_dir = Path(__file__).resolve().parent
    ae_eos_path = package_dir / "statics/AE_EOS" / ref_json
    if ae_eos_path.exists():
        # Open and read the JSON file
        with ae_eos_path.open("rb") as handle:
            data = json.load(handle)
    else:
        raise FileNotFoundError(f"File not found: {ae_eos_path}")

    # import_path = (importlib.resources.path("abinit_fireworks.AE_EOS", ref_json))
    # with import_path as path, open(path, "rb") as handle:
    #    data = json.load(handle)

    bm_fit = data["BM_fit_data"][conf_key]
    ref_v0, ref_b0, ref_b1 = (
        bm_fit["min_volume"],
        bm_fit["bulk_modulus_ev_ang3"],
        bm_fit["bulk_deriv"],
    )

    results = {
        "birch_murnaghan_results": [v0, b0, b1],
        "reference_ae_V0_B0_B1": [ref_v0, ref_b0, ref_b1],
        "V0_B0_B1_units_info": "eV/A^3 for B0",
    }
    # Delta computation
    try:
        delta, deltarel, delta1 = calc_delta(v0, b0, b1, ref_v0, ref_b0, ref_b1)
    except Exception:
        pass
    else:
        results.update(
            {
                "delta": delta,
                "delta1": delta1,
                "delta_unit": "meV/atom",
                "delta_relative": deltarel,
                "delta_relative_unit": "%",
                "natoms": natoms,
                "delta/natoms": delta / natoms,
            }
        )

    # The nu_measure is a measure of the relative error of the fit parameters
    try:
        nu_measure = rel_errors_vec_length(ref_v0, ref_b0, ref_b1, v0, b0, b1)
    except Exception:
        pass
    else:
        results.update({"rel_errors_vec_length": nu_measure})

    return results


if __name__ == "__main__":
    result = metric_analyze('Nd', 'XO', 29.6867157, 1.9282279, -20.1279188, 2)
    # response = run_locally(result)
    pass
