from eos_workflow.eos_workflow import eos_workflow
from jobflow import run_locally

if __name__ == "__main__":
    ecut = 50
    element = 'Li'
    configuration = 'XO2'
    pseudos = "ONCVPSP-PBE-SR-PDv0.4:standard"

    eos = eos_workflow(element, configuration, ecut, pseudos, precision="test")
    result = run_locally(eos, create_folders=True)
    print(result)
