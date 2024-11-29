from atomate2.abinit.sets.base import as_pseudo_table
from abipy.flowtk.psrepos import OncvpspRepo
from eos_workflow.eos_workflow import eos_workflow
from eos_workflow.workflows import eos_workflows, eos_converge_workflows
from jobflow import run_locally
from jobflow_remote import submit_flow, set_run_config

pseudos = "ONCVPSP-PBE-SR-PDv0.4:standard"

ecut = 100
element = 'Rn'
# configurations = ['BCC']

eos = eos_workflows(element, ecut, pseudos, precision="debug")
eos = set_run_config(eos, name_filter="eos_check", worker="lucia_frontend")
eos = set_run_config(eos, name_filter="eos_delta_calculation", worker="lucia_frontend")
eos = set_run_config(eos, name_filter="export_result", worker="lucia_frontend")


# ecuts = [10, 20, 30]
# configuration = "BCC"
# eos = eos_converge_workflows(element, configuration, ecuts, pseudos, precision="test")

# result = run_locally(eos, create_folders=True)
result = submit_flow(eos)
print(result)
print(element)
