from atomate2.abinit.sets.base import as_pseudo_table
from abipy.flowtk.psrepos import OncvpspRepo
from eos_workflow.eos_workflow import eos_workflow
from eos_workflow.workflows import eos_workflows, eos_converge_workflows
from jobflow import run_locally
from jobflow_remote import submit_flow, set_run_config, get_jobstore
from eos_workflow.inspect import collect_results
from eos_workflow.sets import EosDoc

pseudos = "ONCVPSP-PBE-SR-PDv1.0:standard"

ecut = 100
element = 'Eu'
# configurations = ['XO', 'XO2', 'XO3', 'X2O', 'X2O3', 'X2O5']
# configurations = ['BCC', 'FCC', 'SC', 'Diamond']
# small_config = ['BCC', 'Diamond', 'X2O', 'XO', 'XO2', 'XO3']
# rest = ['FCC', 'SC', 'X2O3', 'X2O5']
configurations = ['BCC', 'XO']
eos = eos_workflows(element, ecut, pseudos,
                    precision="debug",
                    volume_scaling_list=[0.94, 0.98, 1.00, 1.02, 1.06],
                    configurations=configurations
                    )
eos = set_run_config(eos, name_filter="eos_check", worker="lucia_frontend")
eos = set_run_config(eos, name_filter="eos_delta_calculation", worker="lucia_frontend")
eos = set_run_config(eos, name_filter="export_result", worker="lucia_frontend")
eos = set_run_config(eos, name_filter="export_single_workflow", worker="lucia_frontend")

"""eos = set_run_config(eos, name_filter="eos_check", worker="manneback_frontend")
eos = set_run_config(eos, name_filter="eos_delta_calculation", worker="manneback_frontend")
eos = set_run_config(eos, name_filter="export_result", worker="manneback_frontend")
eos = set_run_config(eos, name_filter="export_single_workflow", worker="manneback_frontend")"""


# ecuts = [10, 20, 30]
# configuration = "BCC"
# eos = eos_converge_workflows(element, configuration, ecuts, pseudos, precision="test")

# result = run_locally(eos, create_folders=True)
result = submit_flow(eos)
print(result)
print(element)


"""js = get_jobstore()
js.connect()
result = list(js.query({"name": "export_result"}))
data = collect_results('test', 'test', func_name="export_result")
pass"""

