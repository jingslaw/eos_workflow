from eos_workflow.phonon_convergency_workflow import phonon_convergency_workflow
from eos_workflow.sets import get_standard_structure, eos_kpoints_generation, eos_input_generation
from jobflow_remote import submit_flow, set_run_config
from jobflow import run_locally

element = 'Ce'
# ecuts = [125]
ecuts = [50, 60, 70, 80, 90, 100, 125, 150]
configuration = 'BCC'
pseudo = "ONCVPSP-PBE-SR-PDv0.6:standard"
structure = get_standard_structure(element, configuration)
jobs = phonon_convergency_workflow(
    element,
    configuration,
    pseudo,
    ecuts=ecuts,
    qpt_list=[[0.0, 0.0, 0.0]]
)
# jobs = set_run_config(jobs, name_filter="generate_perts", worker="lucia_frontend")
jobs = set_run_config(jobs, name_filter="Merge DDB", worker="lucia_frontend")
jobs = set_run_config(jobs, name_filter="parse_phonon_files", worker="lucia_frontend")
jobs = set_run_config(jobs, name_filter="write_to_file", worker="lucia_frontend")
jobs = set_run_config(jobs, name_filter="store_inputs", worker="lucia_frontend")

# res = run_locally(jobs, create_folders=True)
res = submit_flow(jobs)
print(res)
