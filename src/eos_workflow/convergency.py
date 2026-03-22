import json
import logging
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

logger = logging.getLogger(__name__)


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
):
    structure = get_standard_structure(element, configuration)
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

    check_job = hints_check(results)
    hints_jobs.append(check_job)
    write_job = write_to_file(check_job.output, fname="hints.txt")
    hints_jobs.append(write_job)
    return Flow(hints_jobs, name=f"Hints-{element}-{configuration}-convergency", output=results)


def convergency_hints_workflows(
    element: str,
    configuration: str = 'FCC',
    ecuts: list[float] | None = None,
    pseudos: str = "ONCVPSP-PBE-SR-PDv0.4:standard",
    precision: str = 'hints',
):
    structure = get_standard_structure(element, configuration)
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

    check_job = hints_check(results)
    hints_jobs.append(check_job)
    write_job = write_to_file(check_job.output, fname="hints.txt")
    hints_jobs.append(write_job)
    return Flow(hints_jobs, name=f"Hints-{element}-{configuration}-convergency", output=results)
