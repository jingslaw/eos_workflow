from eos_workflow.convergency import convergency_workflows
from jobflow_remote import submit_flow, set_run_config
from jobflow import run_locally
"""import os

os.environ["PATH"] = (
    "/home/wjing/programs/abinit-9.10.1/build_gfortran/src/98_main:"
    + os.environ.get("PATH", "")
)

os.environ["LD_LIBRARY_PATH"] = (
    "/home/wjing/local/lib:"
    + os.environ.get("LD_LIBRARY_PATH", "")
)"""


element = "H"
# ecuts = [51, 52, 53, 54, 61, 62, 63, 64]
configuration = 'FCC'
pseudo = "ONCVPSP-PBE-SR-PDv2.0:standard"
jobs = convergency_workflows(
    element,
    configuration,
    # ecuts=ecuts,
    pseudos=pseudo,
    precision="debug"
)


jobs = set_run_config(jobs, name_filter="hints_check", worker="lucia_frontend")
jobs = set_run_config(jobs, name_filter="write_to_file", worker="lucia_frontend")
res = submit_flow(jobs)
# res = run_locally(jobs, create_folders=True)
print(res)
print(element)
print(pseudo)
