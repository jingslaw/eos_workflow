from eos_workflow.eos_workflow import eos_workflow
from jobflow import run_locally

if __name__ == "__main__":
    ecut = 50
    element = 'O'
    configuration = 'BCC'
    pseudos = "ONCVPSP-PBE-SR-PDv0.4:standard"

    eos = eos_workflow(element, configuration, ecut, pseudos, precision="standard")
    result = run_locally(eos, create_folders=True)
    print(result)
