import json
from eos_workflow.eos_workflow import eos_workflow
from eos_workflow.utilities import ACWF_CONFIGURATIONS
from eos_workflow.utilities import ELEMENTS_INCLUDE_F_ELECTRONS
from jobflow import Flow, job


@job
def export_result(workflow_results, file_name="eos_fitting_results.json"):
    with open(file_name, 'w') as fp:
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
    return Flow([flow, export_result(flow.output)])


@job
def converge_results(ecuts, outputs, print_raw=True):
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
    result = converge_results(ecuts, flow.output)
    return [flow, result, export_result(result.output, file_name="eos_converge_results.json")]
