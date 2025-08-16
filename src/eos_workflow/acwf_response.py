from atomate2.abinit.jobs.base import BaseAbinitMaker, abinit_job
from atomate2.abinit.schemas.task import AbinitTaskDoc
from atomate2.abinit.schemas.calculation import TaskState
from atomate2.abinit.utils.history import JobHistory
from atomate2.abinit.sets.base import AbinitInputGenerator
from atomate2.abinit.sets.response import PhononSetGenerator

from abipy.flowtk.events import (
    AbinitCriticalWarning,
    ScfConvergenceWarning,
)

from pymatgen.core.structure import Structure
from dataclasses import dataclass, field
from jobflow import Response, Job
from pathlib import Path
import logging

from typing import ClassVar
from collections.abc import Sequence

logger = logging.getLogger(__name__)


@dataclass
class AcwfBaseAbinitMaker(BaseAbinitMaker):
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


@dataclass
class AcwfResponseMaker(AcwfBaseAbinitMaker):
    """Maker for a Response Function ABINIT calculation job.

    Parameters
    ----------
    calc_type : str
        The type of RF.
    name : str
        The job name.
    """

    calc_type: str = "RF"
    name: str = "RF calculation"
    task_document_kwargs: dict = field(default_factory=dict)
    input_set_generator: AbinitInputGenerator
    stop_jobflow_on_failure: bool = True

    CRITICAL_EVENTS: ClassVar[Sequence[AbinitCriticalWarning]] = (
        ScfConvergenceWarning,
    )

    @abinit_job
    def make(
        self,
        structure: Structure | None = None,
        prev_outputs: str | list[str] | None = None,
        restart_from: str | list[str] | None = None,
        history: JobHistory | None = None,
        perturbation: dict | None = None,
    ) -> Job:
        if perturbation:
            self.input_set_generator.factory_kwargs.update(
                {f"{self.calc_type.lower()}_pert": perturbation}
            )

        return super().make.original(
            self,
            structure=structure,
            prev_outputs=prev_outputs,
            restart_from=restart_from,
            history=history,
        )


@dataclass
class AcwfPhononResponseMaker(AcwfResponseMaker):
    """Maker to create a job with a Phonon ABINIT calculation.

    Parameters
    ----------
    name : str
        The job name.
    """

    calc_type: str = "Phonon"
    name: str = "Phonon calculation"
    input_set_generator: AbinitInputGenerator = field(
        default_factory=PhononSetGenerator
    )

    CRITICAL_EVENTS: ClassVar[Sequence[AbinitCriticalWarning]] = (
        ScfConvergenceWarning,
    )
