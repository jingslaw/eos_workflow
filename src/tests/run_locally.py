from eos_workflow.eos_workflow import eos_workflow
from jobflow import run_locally
import os

os.environ["PATH"] = (
    "/home/wjing/programs/abinit-9.10.1/build_gfortran/src/98_main:"
    + os.environ.get("PATH", "")
)

if __name__ == "__main__":
    ecut = 20
    element = 'O'
    configuration = 'Diamond'
    pseudos = "ONCVPSP-PBE-SR-PDv0.4:standard"

    eos = eos_workflow(element, configuration, ecut, pseudos, precision="test")
    result = run_locally(eos, create_folders=True)
    print(result)
