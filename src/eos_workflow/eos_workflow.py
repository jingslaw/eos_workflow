import os
import json
import logging
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path

from atomate2.abinit.jobs.core import StaticMaker
from atomate2.abinit.sets.base import AbinitInputGenerator
from atomate2.abinit.schemas.calculation import TaskState
from atomate2.abinit.schemas.task import AbinitTaskDoc
from atomate2.abinit.utils.history import JobHistory
from pymatgen.core import Structure
from pymatgen.io.abinit.pseudos import PseudoTable
from jobflow import Flow, job, Response

from eos_workflow.sets import get_standard_structure, eos_input_generation, eos_kpoints_generation, EosSetGenerator
from eos_workflow.delta_metric import birch_murnaghan_fit, metric_analyze

logger = logging.getLogger(__name__)


def float2bson(eta: float):
    return str(eta).replace('.', '_')


def bson2float(eta: str):
    return float(eta.replace('_', '.'))


@dataclass
class EosMaker(StaticMaker):
    calc_type: str = "scf"
    name: str = "eos calculation, scaling_factor=1.00"
    input_set_generator: AbinitInputGenerator = field(default_factory=EosSetGenerator)

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


def parallel_eos_calculation(
        element: str,
        configuration: str,
        structure: Structure,
        pseudos: str | list[str] | PseudoTable | None = "ONCVPSP-PBE-SR-PDv0.4:standard",
        abinit_settings=None,
        kpoints_settings=None,
        volume_scaling_list=None
):
    if kpoints_settings is None:
        kpoints_settings = field(default_factory=dict)
    if abinit_settings is None:
        abinit_settings = field(default_factory=dict)
    if volume_scaling_list is None:
        volume_scaling_list = [0.94, 0.96, 0.98, 1.02, 1.04, 1.06]
    elif 1.00 in volume_scaling_list:
        index = volume_scaling_list.index(1.00)
        volume_scaling_list.pop(index)
    eos_jobs = []
    eos_maker = EosMaker()
    volume0 = structure.volume
    eos_maker.input_set_generator.user_abinit_settings = abinit_settings
    eos_maker.input_set_generator.user_kpoints_settings = kpoints_settings
    eos_maker.input_set_generator.pseudos = pseudos
    eos_job = eos_maker.make(structure=structure)
    eos_job.name = f"eos-{element}-{configuration}-1.00"
    eos_jobs.append(eos_job)
    eos_outputs = {float2bson(1.00): eos_job.output}
    for eta in volume_scaling_list:
        structure = structure.copy()
        structure = structure.scale_lattice(volume0 * eta)
        eos_job = eos_maker.make(structure=structure)
        eos_job.name = f"eos-{element}-{configuration}-{eta}"
        eos_jobs.append(eos_job)
        eos_outputs.update({float2bson(eta): eos_job.output})
    flow = Flow(eos_jobs, output=eos_outputs)
    return flow


def get_series_neighbors(eta0, volume_scaling_list=None):
    if volume_scaling_list is None:
        volume_scaling_list = [0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.06]
    if eta0 not in volume_scaling_list:
        return [], []
    index = volume_scaling_list.index(eta0)
    left_neighbors = volume_scaling_list[:index][::-1]
    right_neighbors = volume_scaling_list[index + 1:]
    return left_neighbors, right_neighbors


def series_eos_calculations(
        element: str,
        configuration: str,
        eta: float,
        minimum: AbinitTaskDoc,
        pseudos: str | list[str] | PseudoTable | None = "ONCVPSP-PBE-SR-PDv0.4:standard",
        abinit_settings=None,
        kpoints_settings=None,
        volume_scaling_list=None
):
    if kpoints_settings is None:
        kpoints_settings = field(default_factory=dict)
    if abinit_settings is None:
        abinit_settings = field(default_factory=dict)
    if volume_scaling_list is None:
        volume_scaling_list = [0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.06]
    left, right = get_series_neighbors(eta, volume_scaling_list=volume_scaling_list)
    eos_jobs = []
    structure = minimum.output.structure
    path = minimum.dir_name
    volume0 = structure.volume / eta
    path0 = path
    eos_maker = EosMaker()
    eos_maker.input_set_generator.user_abinit_settings = abinit_settings
    eos_maker.input_set_generator.user_kpoints_settings = kpoints_settings
    eos_maker.input_set_generator.pseudos = pseudos
    eos_outputs = {float2bson(eta): minimum}
    for eta in left:
        structure = structure.copy()
        structure = structure.scale_lattice(volume0 * eta)
        eos_job = eos_maker.make(structure=structure, prev_outputs=path)
        eos_job.name = f"eos-{element}-{configuration}-{eta}"
        path = eos_job.output.dir_name
        eos_jobs.append(eos_job)
        eos_outputs.update({float2bson(eta): eos_job.output})
    path = path0
    for eta in right:
        structure = structure.copy()
        structure = structure.scale_lattice(volume0 * eta)
        eos_job = eos_maker.make(structure=structure, prev_outputs=path)
        eos_job.name = f"eos-{element}-{configuration}-{eta}"
        path = eos_job.output.dir_name
        eos_jobs.append(eos_job)
        eos_outputs.update({float2bson(eta): eos_job.output})
    return Flow(eos_jobs, output=eos_outputs)


@job
def eos_check(eos_jobs_outputs, self_clean=True):
    volumes = []
    energies = []
    scaling_factors = []
    flag = True
    eos_log = ''
    for eta, task_doc in eos_jobs_outputs.items():
        if task_doc.state == TaskState.SUCCESS:
            volumes.append(task_doc.output.structure.volume)
            energies.append(task_doc.output.energy)
            scaling_factors.append(bson2float(eta))
        else:
            flag = False
            eos_log += "\nSome eos calculations are failed.\n"
    if len(energies) == 0:
        print("ERROR: all eos calculations are failed.")
        exit(-1)
    minimum = min(energies)
    index = energies.index(minimum)
    eta = scaling_factors[index]
    minimum_eos_result = eos_jobs_outputs[float2bson(eta)]
    num = minimum_eos_result.output.structure.num_sites
    volume_energy = {
        "num_of_atoms": num,
        "volume_unit": "A^3",
        "energy_unit": "eV",
        "volumes": volumes,
        "energies": energies,
        "scaling_factors": scaling_factors,
        "minimum_eos_result": minimum_eos_result,
        "eta_min": eta,
        "eos_is_converged": flag,
        "eos_failed_problems": eos_log,
    }
    if flag is True:
        fitting_results = birch_murnaghan_fit(volume_energy)
        residuals0 = fitting_results["residuals0"]
        if residuals0 > 0.005:
            # < 0.005 corresponds to the coefficient of determination: R^2 > 0.995, which means fitting is pretty good.
            eos_log += "\nWARNING: The behavior of eos curve is bad.\n"
            print("\nWARNING: The behavior of eos curve is bad.\n")
            flag = False
            volume_energy.update({"eos_is_converged": flag, "eos_failed_problems": eos_log})
    tmp = deepcopy(volume_energy)
    tmp.pop("minimum_eos_result")
    print(tmp)
    if self_clean:
        for eta, task_doc in eos_jobs_outputs.items():
            calc_path = task_doc.dir_name
            try:
                uuid, times = calc_path.split('_')
                times = int(times)

                # Determine range of iterations
                if task_doc.state == TaskState.SUCCESS and task_doc.output.energy == minimum:
                    iterations = range(1, times)
                else:
                    iterations = range(1, times + 1)

                # Remove files for the specified iterations
                for i in iterations:
                    path = f"{uuid}_{i}"
                    wfk_path = os.path.join(path, 'outdata/out_WFK')
                    if os.path.exists(wfk_path):
                        os.remove(wfk_path)
            except ValueError:
                if task_doc.state == TaskState.SUCCESS and task_doc.output.energy == minimum:
                    continue
                else:
                    wfk_path = os.path.join(calc_path, 'outdata/out_WFK')
                    if os.path.exists(wfk_path):
                        os.remove(wfk_path)

    return volume_energy


@job
def series_vs_parallel_results(series_results, parallel_results):
    for i in range(len(series_results["energies"])):
        eta = series_results["scaling_factors"][i]
        energy = series_results["energies"][i]
        if eta in parallel_results["scaling_factors"]:
            j = parallel_results["scaling_factors"].index(eta)
            energy_p = parallel_results["energies"][j]
            if (energy - energy_p) > 0.0001:
                # 0.0001 could be replaced to the value of toldfe, which is better need further test for all pps.
                series_results["eos_failed_problems"] += \
                    (f"\nEtot at scaling_factor={eta} in series eos calculation"
                     f"is larger than corresponding Etot at parallel calculation."
                     f"The start point for series eos calculations might be"
                     f"meta-stable states above the ground state.\n")
                series_results.update({"eos_is_converged": False})
    return series_results


@job
def wfk_calculations(
    element,
    configuration,
    volume_energy_result,
    pseudos,
    abinit_settings,
    kpoints_settings,
    volume_scaling_list
):
    if volume_energy_result["eos_is_converged"]:
        return volume_energy_result
    else:
        minimum = volume_energy_result["minimum_eos_result"]
        eta = volume_energy_result["eta_min"]
        series_eos_job = series_eos_calculations(
            element,
            configuration,
            eta,
            minimum,
            pseudos,
            abinit_settings=abinit_settings,
            kpoints_settings=kpoints_settings,
            volume_scaling_list=volume_scaling_list
        )
        check_job = eos_check(series_eos_job.output)
        compare_job = series_vs_parallel_results(check_job.output, volume_energy_result)
        series_eos_flow = [series_eos_job, check_job, compare_job]
        flow = Flow(series_eos_flow, output=compare_job.output)
        return Response(replace=flow)


@job
def export_single_workflow(element, configuration, result, fname="eos_fitting_results.json"):
    with open(fname, 'w') as fp:
        tmp = {element: {configuration: result}}
        json.dump(tmp, fp, indent=4)
    return tmp

@job
def eos_delta_calculation(
    element: str,
    configuration: str,
    volume_energy_results: dict,
):
    flag = volume_energy_results["eos_is_converged"]
    result = deepcopy(volume_energy_results)
    result.pop("minimum_eos_result")
    print(result)
    if flag is True:
        fitting_results = birch_murnaghan_fit(volume_energy_results)
        v0 = fitting_results["volume0"]
        b0 = fitting_results["bulk_modulus0"]
        b1 = fitting_results["bulk_deriv0"]
        num = volume_energy_results["num_of_atoms"]
        e0 = fitting_results["energy0"]
        result.update({"energy0": e0})
        eos_results = metric_analyze(element, configuration, v0, b0, b1, num)
        print(eos_results)
        result.update(eos_results)
    return result


@job
def get_remote_pseudo_path(pseudos_names):
    pseudos_path = os.getenv('EOS_PSEUDO_PATH')
    pseudos = [os.path.join(pseudos_path, pseudo) for pseudo in pseudos_names]
    return pseudos


def generate_abinit_inputs(element, configuration, ecut, pseudos, structure, precision="standard"):
    user_abinit_settings = eos_input_generation(element, configuration, ecut, pseudos, precision=precision)
    user_kpoints_settings = eos_kpoints_generation(structure, precision=precision)
    abinit_inputs = {
        "user_abinit_settings": user_abinit_settings,
        "user_kpoints_settings": user_kpoints_settings
    }
    return abinit_inputs


def eos_workflow(element, configuration, ecut, pseudos, volume_scaling_list=None, precision='standard'):
    """
    :param element: str, the name of element
    :param configuration: str, it should be one of the config in
    ["XO", "XO2", "XO3", "X2O", "X2O3", "X2O5", "BCC", "FCC", "SC", "Diamond"]
    :param ecut: float, Ecut (unit: Ha)
    :param pseudos: a str like "ONCVPSP-PBE-SR-PDv0.4:standard"
    :param volume_scaling_list: the scaling factors, [0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.06] by default
    :param precision: tag "standard" provides the abinit inputs as same as aiida-common workflow,
    tag "test" provides a test to this jobflow in a lower precision
    :return:
    """
    eos_jobflow = []
    if volume_scaling_list is None:
        volume_scaling_list = [0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.06]

    structure = get_standard_structure(element, configuration)
    inputs_job = generate_abinit_inputs(element, configuration, ecut, pseudos, structure, precision=precision)

    eos_jobs = parallel_eos_calculation(
        element, configuration,
        structure, pseudos,
        abinit_settings=inputs_job["user_abinit_settings"],
        kpoints_settings=inputs_job["user_kpoints_settings"],
        volume_scaling_list=volume_scaling_list.copy()
    )
    eos_jobflow.append(eos_jobs)
    check_job = eos_check(eos_jobs.output)
    eos_jobflow.append(check_job)
    wfk_job = wfk_calculations(
        element,
        configuration,
        check_job.output,
        pseudos,
        abinit_settings=inputs_job["user_abinit_settings"],
        kpoints_settings=inputs_job["user_kpoints_settings"],
        volume_scaling_list=volume_scaling_list.copy()
    )
    eos_jobflow.append(wfk_job)
    delta_job = eos_delta_calculation(element, configuration, wfk_job.output)
    eos_jobflow.append(delta_job)
    write_job = export_single_workflow(element, configuration, delta_job.output)
    eos_jobflow.append(write_job)
    workflow = Flow(eos_jobflow, output=delta_job.output)
    return workflow


if __name__ == "__main__":
    print("test")
