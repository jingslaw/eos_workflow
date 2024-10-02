from eos_workflow.eos_workflow import eos_workflow
from jobflow_remote import submit_flow


if __name__ == "__main__":
    ecut = 20
    element = 'O'
    configuration = 'BCC'
    pseudos = "ONCVPSP-PBE-SR-PDv0.4:standard"

    eos = eos_workflow(element, configuration, ecut, pseudos, precision="debug")
    result = submit_flow(eos)
    # print(result)
