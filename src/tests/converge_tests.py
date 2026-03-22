from eos_workflow.workflows import eos_converge_workflows
from jobflow import run_locally
from jobflow_remote import submit_flow, set_run_config
"""import os

os.environ["PATH"] = (
    "/home/wjing/programs/abinit-9.10.1/build_gfortran/src/98_main:"
    + os.environ.get("PATH", "")
)

os.environ["LD_LIBRARY_PATH"] = (
    "/home/wjing/local/lib:"
    + os.environ.get("LD_LIBRARY_PATH", "")
)"""

pseudos = "ONCVPSP-PBE-SR-PDv0.6:standard"

element = 'Ce'
ecuts = [40, 45]
# ecuts = [50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 150]
configuration = "FCC"
eos = eos_converge_workflows(element, configuration, ecuts, pseudos, precision="debug")
eos = set_run_config(eos, name_filter="eos_check", worker="lucia_frontend")
eos = set_run_config(eos, name_filter="eos_delta_calculation", worker="lucia_frontend")
eos = set_run_config(eos, name_filter="export_result", worker="lucia_frontend")
eos = set_run_config(eos, name_filter="converge_results", worker="lucia_frontend")
eos = set_run_config(eos, name_filter="export_single_workflow", worker="lucia_frontend")

# result = run_locally(eos, create_folders=True)
result = submit_flow(eos)
print(result)
