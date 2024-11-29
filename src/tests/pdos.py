from atomate2.abinit.jobs.core import StaticMaker, LineNonSCFMaker, UniformNonSCFMaker
from atomate2.abinit.flows.core import BandStructureMaker
from jobflow import run_locally

from eos_workflow.eos_workflow import get_standard_structure
from eos_workflow.sets import eos_input_generation, eos_kpoints_generation

# lithium = "/home/wjing/eos_workflow/src/eos_workflow/statics/pseudos/Li-s.psp8"
lithium = "/home/wjing/.abinit/pseudos/ONCVPSP-PBE-SR-PDv0.5/Li/Li-sc.psp8"
oxygen = "/home/wjing/eos_workflow/src/eos_workflow/statics/pseudos/O.psp8"
ecut = 30
pseudos = [lithium, oxygen]
element = "Li"
configuration = "SC"      # "XO"
structure = get_standard_structure(element, configuration)

bands_maker = BandStructureMaker()
input_settings = eos_input_generation(element, configuration, ecut, pseudos, precision="bands")
kpoints_settings = eos_kpoints_generation(structure, precision="bands")
bands_maker.static_maker.input_set_generator.user_abinit_settings = input_settings
bands_maker.static_maker.input_set_generator.user_kpoints_settings = kpoints_settings

bandstructure_flow = BandStructureMaker().make(structure)

# run the flow
run_locally(bandstructure_flow, create_folders=True)
