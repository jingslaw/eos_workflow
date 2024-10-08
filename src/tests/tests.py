from atomate2.abinit.sets.base import as_pseudo_table
from abipy.flowtk.psrepos import OncvpspRepo
from eos_workflow.eos_workflow import eos_workflow
from eos_workflow.workflows import eos_workflows
from jobflow import run_locally
from jobflow_remote import submit_flow

# repo = OncvpspRepo(ps_generator="ONCVPSP", xc_name="PBE", relativity_type="SR", project_name="PD", version="0.5", url="test string")
# pseudos = repo.get_pseudos("standard")
pseudos = "ONCVPSP-PBE-SR-PDv0.5:standard"

ecut = 100
element = 'Pr'
configurations = None

eos = eos_workflows(element, ecut, pseudos, configurations=configurations, precision="debug")
# result = run_locally(eos, create_folders=True)
result = submit_flow(eos)
print(result)
