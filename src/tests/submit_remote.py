from eos_workflow.eos_workflow import eos_workflow
from jobflow_remote import submit_flow

if __name__ == "__main__":
    ecut = 38
    element = 'O'
    configuration = 'BCC'
    oxygen = "O.psp8"
    pseudo_files = [oxygen]

    eos = eos_workflow(element, configuration, ecut, pseudo_files, precision="test")

    result = submit_flow(eos)
    print(result)
