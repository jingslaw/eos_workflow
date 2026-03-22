from atomate2.abinit.jobs.core import RelaxMaker
from eos_workflow.sets import get_standard_structure, eos_kpoints_generation, eos_input_generation
from jobflow_remote import submit_flow

pseudos = "ONCVPSP-PBE-SR-PDv1.0:standard"
element = 'Nd'
configuration = 'RE'
ecut = 100

structure = get_standard_structure(element, configuration)
abinit_settings = eos_input_generation(element, configuration, ecut, pseudos, precision='LnN')
abinit_settings.pop("toldfe")
abinit_settings.update({"tolrff": 0.02})
kpoints_settings = eos_kpoints_generation(structure, precision='phonon')

relax_maker = RelaxMaker()
relax_maker.name = f'{element}N relaxation'
relax_maker.input_set_generator.user_abinit_settings = abinit_settings
relax_maker.input_set_generator.user_kpoints_settings = kpoints_settings
relax_maker.input_set_generator.pseudos = pseudos

relax_job = relax_maker.make(structure=structure)
results = submit_flow(relax_job)
print(results)
