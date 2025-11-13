from atomate2.abinit.sets.core import StaticSetGenerator
from atomate2.abinit.jobs.core import RelaxMaker, StaticMaker
from eos_workflow.sets import get_standard_structure
from jobflow import run_locally, job

@job
def print_stucture(test):
    print(test.output)
    return test.output


element = 'Ce'
ecuts = [40]
# ecuts = [50, 60, 70, 80, 90, 100, 125, 150]
configuration = 'BCC'
pseudo = "ONCVPSP-PBE-SR-PDv0.6:standard"
p = "/home/wjing/.abinit/pseudos/ONCVPSP-PBE-SR-PDv0.6/Si/Si.psp8"
structure = get_standard_structure(element, configuration)

rlx = RelaxMaker.full_relaxation()
rlx.input_set_generator.pseudos = pseudo
job = rlx.make(structure=structure)

# job2 = print_stucture(job.output)
res = run_locally([job], create_folders=True)
print(res)
