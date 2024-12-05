import json

import numpy as np

from eos_workflow.eos_workflow import eos_workflow
from eos_workflow.sets import EosDoc
from eos_workflow.utilities import ACWF_CONFIGURATIONS, ELEMENTS_INCLUDE_F_ELECTRONS
from eos_workflow.delta_metric import birch_murnaghan_fit
from jobflow import Flow, job


def sum_up(workflow_outputs, num_of_volumes=None) -> EosDoc:
    eos_result = EosDoc()
    warning_lines = []
    failed_wfs = []
    all_missing_outputs = {}
    completely_off = []
    all_eos_data = {}
    all_stress_data = {}
    all_BM_fit_data = {}
    num_atoms_in_sim_cell = {}

    if num_of_volumes is None:
        num_of_volumes = 0
        for element, result in workflow_outputs.items():
            for configuration, res in result.items():
                if len(res["volumes"]) > num_of_volumes:
                    num_of_volumes = len(res["volumes"])
        if num_of_volumes == 0:
            print("WARNING: All results are empty")
            return eos_result

    for element, result in workflow_outputs.items():
        for configuration, res in result.items():
            # Initialize to None if the outputs are not there
            eos_data = None
            stress_data = None
            BM_fit_data = None
            num_atoms = res["num_of_atoms"]
            if res["eos_is_converged"] is False:
                if len(res["volumes"]) / float(num_of_volumes) < 0.8:
                    failed_wfs.append({
                        'element': element,
                        'configuration': configuration,
                        'process_state': 'finished',
                        'exit_status': 801,
                    })
                    # Return the info collected so far, eos_data, stress_data, BM_fit_data are still None
                    all_eos_data[f'{element}-{configuration}'] = eos_data
                    num_atoms_in_sim_cell[f'{element}-{configuration}'] = num_atoms
                    all_stress_data[f'{element}-{configuration}'] = stress_data
                    all_BM_fit_data[f'{element}-{configuration}'] = BM_fit_data
                    continue
                else:
                    tmp = num_of_volumes - len(res["volumes"])
                    if tmp > 0:
                        all_missing_outputs[f'{element}-{configuration}'] = tmp
                        warning_lines.append(
                            "  WARNING! MISSING OUTPUTS: {0}".format(all_missing_outputs[f'{element}-{configuration}'])
                        )
            eos_data = list(zip(res["volumes"], res["energies"]))
            eos_data = sorted(eos_data, key=lambda x: x[0])
            energies = [e for _, e in eos_data]

            min_loc = np.array(energies).argmin()
            if min_loc == 0:
                # Side is whether the minimum occurs on the left side (small volumes) or right side (large volumes)
                completely_off.append({'element': element, 'configuration': configuration, 'side': 'left'})
            elif min_loc == len(res["energies"]) - 1:
                completely_off.append({'element': element, 'configuration': configuration, 'side': 'right'})

            try:
                stress_data = list(zip(res["volumes"], res["stresses"]))
                stress_data = sorted(stress_data, key=lambda x: x[0])
            except:
                stress_data = None

            try:
                fitting_results = birch_murnaghan_fit(res)
                residuals0 = fitting_results["residuals0"]
                if residuals0 >= 1000:
                    warning_lines.append(f"WARNING! Unable to fit for {element} {configuration}")
                else:
                    BM_fit_data = {
                        'min_volume': fitting_results["volume0"],
                        'E0': fitting_results['energy0'],
                        'bulk_modulus_ev_ang3': fitting_results["bulk_modulus0"],
                        'bulk_deriv': fitting_results["bulk_deriv0"],
                        'residuals': residuals0
                    }
                    if residuals0 > 1.e-3:
                        warning_lines.append(
                            f"WARNING! High fit residuals: {residuals0} for {element} {configuration}")

            except ValueError:
                warning_lines.append(f"WARNING! Unable to fit for {element}-{configuration}")

            all_eos_data[f'{element}-{configuration}'] = eos_data
            num_atoms_in_sim_cell[f'{element}-{configuration}'] = num_atoms
            all_stress_data[f'{element}-{configuration}'] = stress_data
            all_BM_fit_data[f'{element}-{configuration}'] = BM_fit_data

    eos_result.failed_wfs = failed_wfs
    eos_result.all_missing_outputs = all_missing_outputs
    eos_result.completely_off = completely_off
    eos_result.all_eos_data = all_eos_data
    eos_result.all_stress_data = all_stress_data
    eos_result.num_atoms_in_sim_cell = num_atoms_in_sim_cell
    eos_result.warning_lines = warning_lines
    eos_result.all_BM_fit_data = all_BM_fit_data
    return eos_result


@job
def export_result(workflow_results, file_name="eos_fitting_results.json"):
    with open(file_name, 'w') as fp:
        json.dump(workflow_results, fp, indent=4)
    eos_result = sum_up(workflow_results)
    return eos_result


def eos_workflows(
    element: str,
    ecut: float,
    pseudos: str = "ONCVPSP-PBE-SR-PDv0.4:standard",
    configurations: str | list[str] | None = None,
    volume_scaling_list: list[float] | None = None,
    precision: str = 'standard'
):
    if configurations is None:
        configurations = ACWF_CONFIGURATIONS
    if isinstance(configurations, str):
        configurations = [configurations]
    workflows = []
    outputs = {}
    for configuration in configurations:
        workflow = eos_workflow(
            element,
            configuration,
            ecut,
            pseudos,
            volume_scaling_list=volume_scaling_list,
            precision=precision
        )
        workflows.append(workflow)
        outputs.update({configuration: workflow.output})
    flow = Flow(workflows, output={element: outputs})
    return Flow([flow, export_result(flow.output)], name=f"eos-{element}")


@job
def converge_results(element, configuration, ecuts, outputs, print_raw=True):
    nu = []
    delta = []
    for output in outputs:
        if output["eos_is_converged"]:
            nu.append(output["rel_errors_vec_length"])
            delta.append(output["delta/natoms"])
        else:
            nu.append("NaN")
            delta.append("NaN")
    results = {
        "element": element,
        "configuration": configuration,
        "ecut": ecuts,
        "nu": nu,
        "delta/natoms": delta,
        "delta/natoms unit": "meV/natoms"
    }
    if print_raw is True:
        results.update({"raw": outputs})
    return results


def eos_converge_workflows(
    element: str,
    configuration: str,
    ecuts: list[float] | None = None,
    pseudos: str = "ONCVPSP-PBE-SR-PDv0.4:standard",
    volume_scaling_list: list[float] | None = None,
    precision: str = 'standard'
):
    if ecuts is None:
        if element in ELEMENTS_INCLUDE_F_ELECTRONS:
            ecuts = [50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 120, 150]
        else:
            ecuts = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 45, 50, 60, 75, 100]
    workflows = []
    outputs = []
    for ecut in ecuts:
        workflow = eos_workflow(
            element,
            configuration,
            ecut,
            pseudos,
            volume_scaling_list=volume_scaling_list,
            precision=precision
        )
        workflows.append(workflow)
        outputs.append(workflow.output)
    flow = Flow(workflows, output=outputs)
    result = converge_results(element, configuration, ecuts, flow.output)
    return Flow([flow, result, export_result(result.output, file_name="eos_converge_results.json")],
                name=f"eos-converge-{element}")
