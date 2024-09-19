from eos_workflow.eos_workflow import eos_workflow
from jobflow import run_locally

if __name__ == "__main__":
    ecut = 38
    element = 'O'
    configuration = 'BCC'
    oxygen = "O.psp8"
    pseudo_files = [oxygen]

    eos = eos_workflow(element, configuration, ecut, pseudo_files, precision="test")

    result = run_locally(eos, create_folders=True)
    print(result)
