from eos_workflow.convergency import convergency_workflows, etot_convergency_workflows
from jobflow_remote import submit_flow, set_run_config
from jobflow import run_locally
import os

os.environ["PATH"] = (
    "/home/wjing/programs/abinit-9.10.1/build_gfortran/src/98_main:"
    + os.environ.get("PATH", "")
)

os.environ["LD_LIBRARY_PATH"] = (
    "/home/wjing/local/lib:"
    + os.environ.get("LD_LIBRARY_PATH", "")
)


element = "Ir"
# ecuts = [25, 30, 40]
configuration = 'FCC'
pseudo = "ONCVPSP-PBEsol-FR-PDv0.6:standard"
# jobs = convergency_workflows(
jobs = etot_convergency_workflows(
    element,
    configuration,
    # ecuts=ecuts,
    pseudos=pseudo,
    precision="debug"
)


jobs = set_run_config(jobs, name_filter="hints_check", worker="lucia_frontend")
jobs = set_run_config(jobs, name_filter="write_to_file", worker="lucia_frontend")
jobs = set_run_config(jobs, name_filter="get_hints", worker="lucia_frontend")
jobs = set_run_config(jobs, name_filter="generate_insert_ecuts", worker="lucia_frontend")
jobs = set_run_config(jobs, name_filter="merge_results", worker="lucia_frontend")
res = submit_flow(jobs)
# res = run_locally(jobs, create_folders=True)
print(res)
print(element)
print(pseudo)
