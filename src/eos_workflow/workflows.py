import json
from eos_workflow.eos_workflow import eos_workflow
from eos_workflow.utilities import ACWF_CONFIGURATIONS
from jobflow import Flow, job


@job
def export_result(workflow_results):
    with open("eos_fitting_results.json", 'w') as fp:
        json.dump(workflow_results, fp, indent=4)


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
    return [flow, export_result(flow.output)]
