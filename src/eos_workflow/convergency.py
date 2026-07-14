import json
import logging
import numpy as np
from dataclasses import dataclass, field
from pathlib import Path

from atomate2.abinit.sets.base import AbinitInputGenerator
from atomate2.abinit.sets.core import StaticSetGenerator
from atomate2.abinit.jobs.core import StaticMaker
from atomate2.abinit.jobs.base import BaseAbinitMaker
from atomate2.abinit.schemas.calculation import TaskState
from atomate2.abinit.schemas.task import AbinitTaskDoc
from atomate2.abinit.utils.history import JobHistory
from jobflow import Flow, job, Response

from eos_workflow.utilities import ELEMENTS_INCLUDE_F_ELECTRONS
from eos_workflow.sets import get_standard_structure, eos_kpoints_generation, eos_input_generation
from eos_workflow.workflows import eos_workflow, export_result

logger = logging.getLogger(__name__)


XC = ['LDA', 'PBE', 'PBEsol']


@job
def write_to_file(results: str, fname="output.txt"):
    with open(fname, 'w') as fp:
        json.dump(results, fp, indent=4)
    return results


@job
def hints_check(hints_outputs):
    hints_results = {"energy_unit": "eV"}
    for key, values in hints_outputs.items():
        if isinstance(values, str):
            hints_results[key] = values
        else:
            hints_results[key] = {
                "task_state": values.state,
                "energy": values.output.energy,
            }
    return hints_results


@dataclass
class HintsMaker(StaticMaker):
    calc_type: str = "scf"
    name: str = "hints"
    input_set_generator: AbinitInputGenerator = field(default_factory=StaticSetGenerator)

    def get_response(
        self,
        task_document: AbinitTaskDoc,
        history: JobHistory,
        max_restarts: int = 5,
        prev_outputs: str | tuple | list | Path | None = None,
    ) -> Response:
        """Rewrite StaticMaker.get_response(). Switch off the unconverged error when run_number > max_restarts"""
        if task_document.state == TaskState.SUCCESS:
            return Response(
                output=task_document,
            )

        if history.run_number > max_restarts:
            return Response(
                output=task_document,
            )

        logger.info("Getting restart job.")

        new_job = self.make(
            structure=task_document.structure,
            restart_from=task_document.dir_name,
            prev_outputs=prev_outputs,
            history=history,
        )

        return Response(
            output=task_document,
            replace=new_job,
        )


def convergency_workflows(
    element: str,
    configuration: str = 'FCC',
    ecuts: list[float] | None = None,
    pseudos: str = "ONCVPSP-PBE-SR-PDv0.4:standard",
    precision: str = 'hints',
    prev_results: dict | None = None,
    xc: str | None = None,
):
    if xc is None:
        xc = pseudos.split('-')[1]
        if xc.upper() in XC:
            pass
        else:
            print(f'ERROR: UNKNOWN XC format: {xc}')
            exit(-1)

    structure = get_standard_structure(element, configuration, xc=xc)
    if ecuts is None:
        if element in ELEMENTS_INCLUDE_F_ELECTRONS:
            ecuts = [50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 150]
        else:
            ecuts = [15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100]
    hints_jobs = []
    results = {
        "element": element,
        "configuration": configuration,
        "pseudos": pseudos,
    }
    for ecut in ecuts:
        abinit_settings = eos_input_generation(element, configuration, ecut, pseudos, precision=precision)
        kpoints_settings = eos_kpoints_generation(structure, precision=precision)

        hints_maker = HintsMaker()
        hints_maker.input_set_generator.user_abinit_settings = abinit_settings
        hints_maker.input_set_generator.user_kpoints_settings = kpoints_settings
        hints_maker.input_set_generator.pseudos = pseudos
        hints_job = hints_maker.make(structure=structure)
        hints_job.name = f"hints-{element}-{configuration}-{ecut}"
        hints_jobs.append(hints_job)

        results[f"ecut-{ecut}"] = hints_job.output

    if prev_results:
        merge_jobs = merge_results(results, prev_results)
        hints_jobs.append(merge_jobs)
        check_job = hints_check(merge_jobs.output)
    else:
        check_job = hints_check(results)
    hints_jobs.append(check_job)
    write_job = write_to_file(check_job.output, fname="hints.txt")
    hints_jobs.append(write_job)
    return Flow(hints_jobs, name=f"Hints-{element}-{configuration}-convergency", output=results)


@job
def get_hints(flow_outputs, thresholds):
    data = {"energy_unit": "eV"}
    for key, values in flow_outputs.items():
        if isinstance(values, str):
            data[key] = values
        else:
            data[key] = {
                "task_state": values.state,
                "energy": values.output.energy,
            }

    ecut_keys = sorted(
        [k for k in data if k.startswith("ecut-")],
        key=lambda x: int(x.split("-")[1]),
    )

    ecut = np.array([int(k.split("-")[1]) for k in ecut_keys])
    energy = np.array([data[k]["energy"] for k in ecut_keys])

    ref_energy = energy[-1]
    delta = energy - ref_energy
    abs_delta = np.abs(delta) * 1000

    conv_ecut = {}
    for label, thr in thresholds.items():
        conv_idx = None
        for i in range(len(abs_delta)):
            if np.all(abs_delta[i:] < thr):
                conv_idx = i
                break
        conv_ecut[label] = ecut[conv_idx] if conv_idx is not None else None
    return conv_ecut


@job
def generate_insert_ecuts(ecuts_old, hints):
    ecuts_new = []
    ecuts_old = sorted(ecuts_old)
    for label, hint in hints.items():
        idx = ecuts_old.index(hint)
        if idx > 0:
            tmp = [i for i in range(ecuts_old[idx - 1], ecuts_old[idx])]
            tmp.pop(0)
            ecuts_new += tmp
    ecuts_new = sorted(set(ecuts_new))
    if len(ecuts_new) == 0:
        ecuts_new = [max(5, ecuts_old[0]-5)]
    return ecuts_new


@job
def merge_results(dict1, dict2):
    dict1.update(dict2)
    return dict1


@job
def decorate_convergency_workflows(
    element: str,
    configuration: str = 'FCC',
    ecuts: list[float] | None = None,
    pseudos: str = "ONCVPSP-PBE-SR-PDv0.4:standard",
    precision: str = 'hints',
    prev_results: dict | None = None,
    xc: str | None = None,
):
    if xc is None:
        xc = pseudos.split('-')[1]
        if xc.upper() in XC:
            pass
        else:
            print(f'ERROR: UNKNOWN XC format: {xc}')
            exit(-1)

    flow = convergency_workflows(
        element=element,
        configuration=configuration,
        ecuts=ecuts,
        pseudos=pseudos,
        precision=precision,
        prev_results=prev_results,
        xc=xc,
    )
    return Response(addition=flow, output=flow.output)


def etot_convergency_workflows(
    element: str,
    configuration: str = 'FCC',
    ecuts: list[float] | None = None,
    pseudos: str = "ONCVPSP-PBE-SR-PDv0.4:standard",
    precision: str = 'hints',
    thresholds: dict | None = None,
    xc: str | None = None,
):
    if xc is None:
        xc = pseudos.split('-')[1]
        if xc.upper() in XC:
            pass
        else:
            print(f'ERROR: UNKNOWN XC format: {xc}')
            exit(-1)

    if ecuts is None:
        if element in ELEMENTS_INCLUDE_F_ELECTRONS:
            ecuts = [50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 150]
        else:
            ecuts = [15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100]

    if thresholds is None:
        # threshoulds in meV
        thresholds = {
            "low": 10.0,
            "normal": 5.0,
            "high": 2.0,
        }

    flow0 = convergency_workflows(
        element=element,
        configuration=configuration,
        pseudos=pseudos,
        ecuts=ecuts,
        precision=precision,
        xc=xc
    )

    hints_job = get_hints(flow0.output, thresholds)
    insert_job = generate_insert_ecuts(ecuts, hints_job.output)
    flow1 = decorate_convergency_workflows(
        element=element,
        configuration=configuration,
        pseudos=pseudos,
        ecuts=insert_job.output,
        precision=precision,
        prev_results=flow0.output,
        xc=xc
    )

    return Flow([flow0, hints_job, insert_job, flow1], name=f"Hints-{element}-{configuration}-convergency", output=flow1.output)


@job
def converge_results(element, configuration, ecuts, pseudos, outputs, print_raw=False):
    results = {
        "element": element,
        "configuration": configuration,
        "pseudos": pseudos,
    }
    ref = None
    i = 0
    for output in outputs:
        tmp = {
            "eos_is_converged": False,
            'nu': None,
            'delta/natoms': None,
            'delta1': None,
            'v0_b0_b1': None,
            'energies': None,
        }
        if output["eos_is_converged"]:
            tmp['nu'] = output["rel_errors_vec_length"]
            tmp['delta/natoms'] = output["delta/natoms"]
            tmp['eos_is_converged'] = output["eos_is_converged"]
            tmp['v0_b0_b1'] = output["birch_murnaghan_results"]
            tmp['delta1'] = output['delta1']
            tmp['energies'] = output['energies']
        results.update({f'ecut-{ecuts[i]}': tmp})
        i += 1
    try:
        results.update({
            "reference_ae_V0_B0_B1": output["reference_ae_V0_B0_B1"],
            "delta/natoms unit": "meV/natoms",
            "num_of_atoms": output["num_of_atoms"],
            "volumes": output["volumes"],
            "scaling_factors": output["scaling_factors"],
            "volume_unit": output["volume_unit"],
            "energy_unit": output["energy_unit"],
        })
    except KeyError:
        print('ERROR in the common part for the convergency results.')
    if print_raw is True:
        results.update({"raw": outputs})
    return results


@job
def merge_results(dict1, dict2):
    dict1.update(dict2)
    return dict1


def eos_converge_workflows(
    element: str,
    configuration: str,
    ecuts: list[float] | None = None,
    pseudos: str = "ONCVPSP-PBE-SR-PDv0.4:standard",
    volume_scaling_list: list[float] | None = None,
    precision: str = 'standard',
    prev_results: dict | None = None,
    xc: str | None = None,
):
    if xc is None:
        xc = pseudos.split('-')[1]
        if xc.upper() in XC:
            pass
        else:
            print(f'ERROR: UNKNOWN XC format: {xc}')
            exit(-1)

    if ecuts is None:
        if element in ELEMENTS_INCLUDE_F_ELECTRONS:
            ecuts = [50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 150]
        else:
            ecuts = [15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100]
    workflows = []
    outputs = []
    for ecut in ecuts:
        workflow = eos_workflow(
            element,
            configuration,
            ecut,
            pseudos,
            volume_scaling_list=volume_scaling_list,
            precision=precision,
            xc=xc
        )
        workflows.append(workflow)
        outputs.append(workflow.output)
    flow = Flow(workflows, output=outputs)
    result_job = converge_results(element, configuration, ecuts, pseudos, flow.output)
    if prev_results:
        merge_job = merge_results(result_job.output, prev_results)
        return Flow(
            [flow, result_job, merge_job, export_result(merge_job.output, file_name="eos_converge_results.json")],
            name=f"eos-converge-{element}", output=merge_job.output)
    return Flow([flow, result_job, export_result(result_job.output, file_name="eos_converge_results.json")],
                name=f"eos-converge-{element}", output=result_job.output)


@job
def decorate_eos_convergency_workflow(
    element: str,
    configuration: str,
    ecuts: list[float] | None = None,
    pseudos: str = "ONCVPSP-PBE-SR-PDv0.4:standard",
    volume_scaling_list: list[float] | None = None,
    precision: str = 'standard',
    prev_results: dict | None = None,
    xc: str | None = None,
):
    if xc is None:
        xc = pseudos.split('-')[1]
        if xc.upper() in XC:
            pass
        else:
            print(f'ERROR: UNKNOWN XC format: {xc}')
            exit(-1)

    flow = eos_converge_workflows(
        element=element,
        configuration=configuration,
        ecuts=ecuts,
        pseudos=pseudos,
        volume_scaling_list=volume_scaling_list,
        precision=precision,
        prev_results=prev_results,
        xc=xc
    )
    return Response(addition=flow, output=flow.output)


@job
def get_eos_hints(data, thresholds):
    ecut_keys = sorted(
        [k for k in data if k.startswith("ecut-")],
        key=lambda x: float(x.split("-")[1]),
    )
    ecut = []
    delta1 = []
    v0_b0_b1 = data["reference_ae_V0_B0_B1"]
    for k in ecut_keys:
        if not isinstance(data[k]["delta/natoms"], float):
            continue
        ecut.append(float(k.split("-")[1]))
        delta1.append(data[k]["delta/natoms"])

    ecut = np.array(ecut)
    delta1 = np.array(delta1)

    v0 = v0_b0_b1[0]
    b0 = v0_b0_b1[1]
    # delta1 (delta_prime) = v0_ref * b0_ref / (v0_AE * b0_AE), where v0_ref = 30 bohr^3, b0_ref = 100 Gpa.
    # so, 2.7747 is the multiple of v0_ref * b0_ref with a unit transform from bohr^3 * Gpa to eV
    delta1 = delta1 * 2.7747 / (v0 * b0)

    ref_delta1 = delta1[-1]
    delta = delta1 - ref_delta1
    abs_delta = np.abs(delta)

    conv_ecut = {}
    for label, thr in thresholds.items():
        conv_idx = None
        for i in range(len(abs_delta)):
            if np.all(abs_delta[i:] < thr):
                conv_idx = i
                break
        conv_ecut[label] = ecut[conv_idx] if conv_idx is not None else None
    print(conv_ecut)
    return conv_ecut


@job
def generate_insert_ecuts(ecuts_old, hints):
    ecuts_new = []
    ecuts_old = sorted(ecuts_old)
    for label, hint in hints.items():
        idx = ecuts_old.index(hint)
        if idx > 0:
            tmp = [i for i in range(ecuts_old[idx - 1], ecuts_old[idx])]
            tmp.pop(0)
            ecuts_new += tmp
    ecuts_new = sorted(set(ecuts_new))
    if len(ecuts_new) == 0:
        ecuts_new = [max(5, ecuts_old[0]-5)]
    print(ecuts_new)
    return ecuts_new


def eos_convergency_workflows(
    element: str,
    configuration: str,
    ecuts: list[float] | None = None,
    pseudos: str = "ONCVPSP-PBE-SR-PDv0.4:standard",
    volume_scaling_list: list[float] | None = None,
    precision: str = 'standard',
    thresholds: dict | None = None,
    xc: str | None = None
):
    if xc is None:
        xc = pseudos.split('-')[1]
        if xc.upper() in XC:
            pass
        else:
            print(f'ERROR: UNKNOWN XC format: {xc}')
            exit(-1)

    if ecuts is None:
        if element in ELEMENTS_INCLUDE_F_ELECTRONS:
            ecuts = [50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 150]
        else:
            ecuts = [15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100]

    flow0 = eos_converge_workflows(
        element=element,
        configuration=configuration,
        ecuts=ecuts,
        pseudos=pseudos,
        volume_scaling_list=volume_scaling_list,
        precision=precision,
        xc=xc,
    )

    if thresholds is None:
        # threshoulds in meV
        thresholds = {
            "low": 0.5,
            "normal": 0.3,
            "high": 0.1,
        }

    hints_job = get_eos_hints(flow0.output, thresholds)
    insert_job = generate_insert_ecuts(ecuts, hints_job.output)

    flow1 = decorate_eos_convergency_workflow(
        element=element,
        configuration=configuration,
        pseudos=pseudos,
        ecuts=insert_job.output,
        volume_scaling_list=volume_scaling_list,
        precision=precision,
        prev_results=flow0.output,
        xc=xc,
    )

    return Flow([flow0, hints_job, insert_job, flow1], name=f"eos-converge-{element}", output=flow1.output)


