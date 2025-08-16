from eos_workflow.phonon_convergency_workflow import PhononConvergencyMaker
from eos_workflow.sets import get_standard_structure, eos_kpoints_generation, eos_input_generation
from jobflow_remote import submit_flow, set_run_config
from jobflow import run_locally
element = 'Si'
ecuts = [30]
# ecuts = [60, 70, 80, 90, 100, 125, 150]
configuration = 'Diamond'
pseudo = "ONCVPSP-PBE-SR-PDv0.4:standard"
structure = get_standard_structure(element, configuration)

for ecut in ecuts:
    pmaker = PhononConvergencyMaker(
        # run_anaddb=True
        qpt_list=[[0.5, 0.5, 0.5], [0.0, 0.0, 0.0]],
    )

    abinit_settings = eos_input_generation(element, configuration, ecut, pseudo, precision='debug')
    # kpoints_settings = eos_kpoints_generation(structure, precision='phonon')

    pmaker.static_maker.input_set_generator.user_abinit_settings = abinit_settings
    # pmaker.static_maker.input_set_generator.user_kpoints_settings = kpoints_settings
    pmaker.static_maker.input_set_generator.pseudos = pseudo
    #n pmaker.phonon_maker.input_set_generator.user_abinit_settings = {"toldfe": 5e-19}

    jobs = pmaker.make(structure=structure)
    jobs = set_run_config(jobs, name_filter="generate_perts", worker="lucia_frontend")
    jobs = set_run_config(jobs, name_filter="Merge DDB", worker="lucia_frontend")
    # result = run_locally(jobs, create_folders=True)
    result = submit_flow(jobs)
    print(result)
