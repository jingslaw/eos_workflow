from eos_workflow.eos_workflow import eos_workflow
from jobflow import run_locally

if __name__ == "__main__":
    ecut = 38
    element = 'O'
    configuration = 'BCC'
    lanthanide = "/home/wjing/eos_workflow/statics/pseudos/Sm-4f.psp8"
    oxygen = "/home/wjing/eos_workflow/statics/pseudos/O.psp8"
    lithium = "O.psp8"
    # pseudo_files = [lanthanide, oxygen]
    pseudo_files = [lithium]

    eos = eos_workflow(element, configuration, ecut, pseudo_files, precision="test")
    # result = submit_flow(eos, worker='lucia')

    result = run_locally(eos, create_folders=True)
    print(result)
