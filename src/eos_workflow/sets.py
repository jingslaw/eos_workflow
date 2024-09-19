from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from ase import io

from abipy.abio.input_tags import SCF
from atomate2.abinit.jobs.core import StaticMaker
from atomate2.abinit.sets.core import StaticSetGenerator
from atomate2.abinit.sets.base import AbinitInputGenerator
from pymatgen.io.abinit.abiobjects import KSampling
from pymatgen.io.ase import AseAtomsAdaptor
from math import ceil

from eos_workflow.utilities import (
    ELEMENTS_INCLUDE_F_ELECTRONS,
    ATOM_NUMBERS_IN_CONFIG,
    UNARIE_CONFIGURATIONS,
    OXIDE_CONFIGURATIONS,
    OncvpspInput
)

EOS_PREV_OUTPUTS_DEPS = (f"{SCF}:WFK",)


@dataclass
class EosSetGenerator(StaticSetGenerator):
    calc_type: str = "static"
    prev_outputs_deps: tuple = EOS_PREV_OUTPUTS_DEPS


@dataclass
class EosMaker(StaticMaker):
    calc_type: str = "scf"
    name: str = "eos calculation, scaling_factor=1.00"
    input_set_generator: AbinitInputGenerator = field(default_factory=EosSetGenerator)


def nband_calculation(element, configuration, pseudos):
    nband = 0
    numbers = ATOM_NUMBERS_IN_CONFIG[configuration]
    if isinstance(pseudos, str):
        pseudos = [pseudos]
    pps = [OncvpspInput(file_path=pseudo) for pseudo in pseudos]
    elements = [pp.basic_setting["atsym"] for pp in pps]
    ven = [pp.valence_electron_numbers for pp in pps]
    if element not in elements:
        print(f"ERROR: Pseudos list {elements} does not have pseudos of element {element} ")
        exit(-1)
    if configuration in UNARIE_CONFIGURATIONS:
        i = elements.index(element)
        nbnd = ceil(ven[i] * numbers[0] / 2)
        nband += max(ceil(nbnd * 1.2), nbnd + 4)
    else:
        if 'O' not in elements:
            print(f"ERROR: oxides {configuration} does not have pseudos of O ")
            exit(-1)
        i = elements.index(element)
        nband += ven[i] * numbers[0]
        i = elements.index('O')
        nband += ven[i] * numbers[1]
        nband = ceil(nband / 2)
    return nband


def eos_input_generation(element, configuration, ecut, pseudos, precision='standard'):
    if element in ELEMENTS_INCLUDE_F_ELECTRONS:
        if configuration in UNARIE_CONFIGURATIONS:
            nband = 20
        else:
            nband = ceil(1.5 * nband_calculation(element, configuration, pseudos))
    else:
        nband = nband_calculation(element, configuration, pseudos)
    numbers = ATOM_NUMBERS_IN_CONFIG[configuration]
    natom = sum(numbers)
    if precision == 'standard':
        eos_settings = {"ecut": ecut, "nband": nband, "tsmear": 2.25e-3, "nstep": 100, "toldfe": 5.0e-11 * natom}
    else:
        eos_settings = {"ecut": ecut, "nband": nband}
    return eos_settings


def eos_kpoints_generation(structure, precision='standard'):
    # unit of kpoint_distance is A^-1
    if precision == 'standard':
        kpoints_distance = 0.06
    else:
        kpoints_distance = 0.5
    reciprocal_lattice_length = structure.lattice.reciprocal_lattice.abc
    ngkpt = [ceil(x / kpoints_distance) for x in reciprocal_lattice_length]
    eos_kpoints_settings = KSampling.monkhorst(ngkpt, shiftk=(0.0, 0.0, 0.0))
    return eos_kpoints_settings


def get_standard_structure(element, configuration):
    # If for delta measure workflow
    base_structure_module = "eos_workflow.statics.structures"

    # uppercase configuration
    if configuration == "LAN":
        # use LAN-nitrides of Wenzovich paper
        # from typical (gs) cif folder
        p_ctx = resources.path(f"{base_structure_module}.gs", f"{element}N.cif")
    elif configuration == "GS":
        p_ctx = resources.path(f"{base_structure_module}.gs", f"{element}.cif")
    # For elements that are verified in nat.rev.phys.2024 paper, use the XSF files.
    # https://github.com/aiidateam/acwf-verification-scripts/tree/main/0-preliminary-do-not-run
    elif configuration in OXIDE_CONFIGURATIONS:
        p_ctx = resources.path(
            f"{base_structure_module}.oxides", f"{element}-{configuration}.xsf"
        )
    elif configuration in UNARIE_CONFIGURATIONS:
        # The xsf file named as E-Diamond.xsf
        if configuration == "DC":
            configuration = "Diamond"
        p_ctx = resources.path(
            f"{base_structure_module}.unaries", f"{element}-{configuration}.xsf"
        )
    else:
        raise ValueError(f"Unknown configuration {configuration}")

    with p_ctx as path_h:
        filepath = Path(path_h)

        if filepath.suffix not in [".cif", ".xsf"]:
            raise ValueError(f"Unknown file type {filepath.suffix}")
        else:
            ase_structure = io.read(filepath)

    return AseAtomsAdaptor.get_structure(ase_structure)


if __name__ == "__main__":
    lanthanide = "/home/wjing/PycharmProjects/abinit_fireworks/Nd-4f.psp8"
    oxygen = "/home/wjing/PycharmProjects/abinit_fireworks/O.psp8"
    pseudo_files = [lanthanide, oxygen]
    result = eos_input_generation('Nd', 'XO2', 100, pseudo_files)
    print(result)
