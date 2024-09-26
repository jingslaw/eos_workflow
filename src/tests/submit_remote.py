from eos_workflow.eos_workflow import eos_workflow
from jobflow_remote import submit_flow, set_run_config
from jobflow import Flow, run_locally
from pymatgen.io.abinit.pseudos import PseudoParser
from atomate2.abinit.sets.base import AbinitInputGenerator, as_pseudo_table

if __name__ == "__main__":
    ecut = 50
    element = 'O'
    configuration = 'BCC'
    pseudos = "ONCVPSP-PBE-SR-PDv0.4:standard"

    eos = eos_workflow(element, configuration, ecut, pseudos, precision="standard")
    result = submit_flow(eos)
    print(result)
