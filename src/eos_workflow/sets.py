from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from ase import io

from abipy.abio.input_tags import SCF
from atomate2.abinit.sets.core import StaticSetGenerator
from atomate2.abinit.sets.base import as_pseudo_table
from pymatgen.io.abinit.abiobjects import KSampling
from pymatgen.io.ase import AseAtomsAdaptor
from pydantic import BaseModel, Field
from math import ceil

from eos_workflow.utilities import (
    ELEMENTS_INCLUDE_F_ELECTRONS,
    ATOM_NUMBERS_IN_CONFIG,
    UNARIE_CONFIGURATIONS,
    OXIDE_CONFIGURATIONS,
)

EOS_PREV_OUTPUTS_DEPS = (f"{SCF}:WFK",)


class EosDoc(BaseModel):
    all_missing_outputs: dict = Field(None, description="Number of missing outputs at X-config eos workflow")
    warning_lines: list = Field(None, description="")
    failed_wfs: list = Field(None, description=" A list of dictionaries with information on the workchains "
                                               "that did not finish with a 0 exit code")
    completely_off: list = Field(None, description="A list of dictionaries that indicate which elements and "
                                                   "configurations have been computed completely off-centre (meaning "
                                                   "that the minimum of all computed energies is on either of the "
                                                   "two edges, i.e. for the smallest or largest volume)")
    all_eos_data: dict = Field(None, description="Dictionary with the EOS data (volumes and energies datapoints)."
                                                 " The keys are the same as the `uuid_mapping`. Values can be None")
    all_stress_data: dict = Field(None, description="")
    all_BM_fit_data: dict = Field(None, description="Birch-Murnaghan fit data. See above for the keys. "
                                                    "Can be None.")
    num_atoms_in_sim_cell: dict = Field(None, description="")


@dataclass
class EosSetGenerator(StaticSetGenerator):
    calc_type: str = "static"
    prev_outputs_deps: tuple = EOS_PREV_OUTPUTS_DEPS


def nband_calculation(element, configuration, pseudos):
    numbers = ATOM_NUMBERS_IN_CONFIG[configuration]
    pseudos = as_pseudo_table(pseudos)
    try:
        pps = pseudos.pseudo_with_symbol(element)
    except AttributeError:
        print(f"ERROR: PseudoTable does not have pseudos of element {element} ")
        exit(-1)
    nband = pps.Z_val * numbers[0] / 2
    if configuration in OXIDE_CONFIGURATIONS:
        try:
            pps = pseudos.pseudo_with_symbol('O')
        except AttributeError:
            print(f"ERROR: PseudoTable does not have pseudos of element O ")
            exit(-1)
        nband += pps.Z_val * numbers[1] / 2
    return ceil(nband)


def eos_input_generation(element, configuration, ecut, pseudos, precision='standard'):
    nband = nband_calculation(element, configuration, pseudos)
    if element in ELEMENTS_INCLUDE_F_ELECTRONS:
        if configuration in UNARIE_CONFIGURATIONS and configuration != "Diamond":
            nband = 20
        else:
            nband = ceil(1.5 * nband)
    else:
        nband = ceil(max(1.2 * nband, nband + 4))
    numbers = ATOM_NUMBERS_IN_CONFIG[configuration]
    natom = sum(numbers)
    if precision == 'standard':
        eos_settings = {
            "ecut": ecut,
            "nband": nband,
            "tsmear": 2.25e-3,
            "nstep": 100,
            "tolvrs": 5.0e-11 * natom,
            "chkprim": 0,
            "chksymbreak": 0,
            "autoparal": 1,
            "optcell": 0,
            "ionmov": 22,
            "dilatmx": 1.0,
            "ecutsm": 0.0,
            "occopt": 3,
            "nspinor": 1,
            "nsppol": 1,
            "nspden": 1,
        }
    elif precision == "debug":
        eos_settings = {
            "ecut": ecut,
            "nband": nband,
            "nstep": 100,
            "tsmear": 2.25e-3,
            "toldfe": 5.0e-11 * natom,
            "chkprim": 0,
            "chksymbreak": 0,
            # "autoparal": 1,
            "nspinor": 1,
            "nsppol": 1,
            "nspden": 1,
        }
    elif precision == "bands":
        eos_settings = {
            "ecut": ecut,
            "nband": nband,
            "nsppol": 1,
            "nspinor": 1,
            "nspden": 1,
            "occopt": 3,
            "tsmear": 0.01,
            "tolvrs": 1e-8 * natom,
        }
    elif precision == "molecule":
        charge = {'X2O': -7, 'XO': -10, 'X2O3': -9, 'XO2': -12, 'X2O5': -7, 'XO3': -6}
        eos_settings = {
            "ecut": ecut,
            "nstep": 100,
            # "tsmear": 2.25e-3,
            # "toldfe": 5.0e-11 * natom,
            "chkprim": 0,
            "chksymbreak": 0,
            # "autoparal": 1,
            "nspinor": 1,
            "nsppol": 1,
            "nspden": 1,
            "diemac": 3.0,
            'charge': charge[configuration]
        }
    else:
        eos_settings = {"ecut": ecut, "nband": nband, "nstep": 100, "nspinor": 1, "nsppol": 1, "nspden": 1}
    return eos_settings


def eos_kpoints_generation(structure, precision='standard'):
    # unit of kpoint_distance is A^-1
    if precision == 'standard':
        kpoints_distance = 0.06
    elif precision == 'debug':
        kpoints_distance = 0.06
    elif precision == 'bands':
        kpoints_distance = 0.2
    elif precision == 'molecule':
        kpoints_distance = 100
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
