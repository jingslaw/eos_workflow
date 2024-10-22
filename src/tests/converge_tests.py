from eos_workflow.workflows import eos_converge_workflows
from jobflow import run_locally
from jobflow_remote import submit_flow, set_run_config

pseudos = "ONCVPSP-PBE-SR-PDv0.5:standard"

element = 'Nd'
ecuts = [50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
configuration = "XO3"
eos = eos_converge_workflows(element, configuration, ecuts, pseudos, precision="standard")
eos = set_run_config(eos, name_filter="eos_check", worker="lucia_frontend")
eos = set_run_config(eos, name_filter="eos_delta_calculation", worker="lucia_frontend")
eos = set_run_config(eos, name_filter="export_result", worker="lucia_frontend")

# result = run_locally(eos, create_folders=True)
result = submit_flow(eos)
print(result)
