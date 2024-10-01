import io
import re
import os.path

from eos_workflow.method import (
    count_nc_nv,
    atom_name,
    electron_configuration,
    find_max_l_in_valence,
    sort_nlf_states,
    nlf2num,
    ATOMIC_SYMBOLS
)

LANTHANIDE_ELEMENTS = [
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
]

ACTINIDE_ELEMENTS = [
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr",
]

ALL_ELEMENTS = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',
]

ELEMENTS_INCLUDE_F_ELECTRONS = LANTHANIDE_ELEMENTS + ACTINIDE_ELEMENTS

OXIDE_CONFIGURATIONS = ["XO", "XO2", "XO3", "X2O", "X2O3", "X2O5"]
UNARIE_CONFIGURATIONS = ["BCC", "FCC", "SC", "Diamond"]
ACWF_CONFIGURATIONS = OXIDE_CONFIGURATIONS + UNARIE_CONFIGURATIONS

ATOM_NUMBERS_IN_CONFIG = {'FCC': (1, ), 'BCC': (1, ), 'SC': (1, ), "Diamond": (2, ),
                          "X2O": (2, 1), "XO": (1, 1), "X2O3": (4, 6), "XO2": (1, 2), "X2O5": (4, 10), "XO3": (1, 3)
                          }


class Lpsp(object):
    def __init__(self, ll, rc=3.0, ep=0.0, ncon=4, nbas=9, qcut=4.5):
        self.ll = ll
        self.rc = rc
        self.ep = ep
        self.ncon = ncon
        self.nbas = nbas
        self.qcut = qcut

    def __str__(self):
        return f"{self.ll} {self.rc} {self.ep} {self.ncon} {self.nbas} {self.qcut}"


class Lproj(object):
    def __init__(self, ll, nproj=2, debl=4):
        self.ll = ll
        self.nproj = nproj
        self.debl = debl

    def __str__(self):
        return f"{self.ll} {self.nproj} {self.debl}"


class OncvpspInput(object):

    def __init__(self, Z=None, file_path=None):
        if file_path:
            lines = []
            with io.open(file_path, "rt", encoding="latin-1") as fp:
                for i, line in enumerate(fp):
                    line = line.strip()
                    lines.append(line)

            header = "# atsym  z    nc    nv    iexc   psfile"
            for i, line in enumerate(lines):
                if re.match(r'^.*atsym.*', line):
                    values = lines[i + 1].split()
                    if len(values) != len(header[1:].split()):
                        print('ERROR: line{0}: {1} needs {2} parameters but only {3} input found'
                              .format(i + 1, header, len(header[1:].split()), len(values)))
                        exit(-1)
                    keys = header[1:].split()
                    values = [int(v) if v.isdigit() else v for v in values]
                    tmp_dict = dict(zip(keys, values))
                    self._basic_setting = tmp_dict
                    break

            header = "# n l f"
            for i, line in enumerate(lines):
                if re.match(r'^.*n\s*l\s*f\s*', line):
                    tmp_list = []
                    for j in range(i + 1, len(lines)):
                        if lines[j].startswith('#'):
                            break
                        substrings = lines[j].split()
                        values = [float(substring) if '.' in substring else int(substring) for substring in substrings]
                        if len(values) != len(header[1:].split()):
                            print('ERROR: line{0}: {1} needs {2} parameters but only {3} input found'
                                  .format(i + 1, header, len(header[1:].split()), len(values)))
                            exit(-1)
                        tmp_list.append(values)
                    self._e_config = tmp_list
                    break

            header = "# lmax"
            for i, line in enumerate(lines):
                if re.match(r'^.*lmax\s*$', line):
                    values = lines[i + 1].split()
                    if len(values) != len(header[1:].split()):
                        print('ERROR: line{0}: {1} needs {2} parameters but only {3} input found'
                              .format(i + 1, header, len(header[1:].split()), len(values)))
                        exit(-1)
                    self._lmax = int(values[0])
                    break

            header = "# l rc ep ncon nbas qcut"
            for i, line in enumerate(lines):
                if re.match(r'^.*ncon.*', line):
                    tmp_list = []
                    for j in range(i + 1, len(lines)):
                        if lines[j].startswith('#'):
                            break
                        substrings = lines[j].split()
                        values = [float(substring) if '.' in substring else int(substring) for substring in substrings]
                        if len(values) != len(header[1:].split()):
                            print('ERROR: line{0}: {1} needs {2} parameters but only {3} input found'
                                  .format(i + 1, header, len(header[1:].split()), len(values)))
                            exit(-1)
                        tmp_list.append(Lpsp(*values))
                    self._lpsp = tmp_list
                    break

            header = "# lloc lpopt rc5 dvloc0"
            for i, line in enumerate(lines):
                if re.match(r'^.*lpopt.*', line):
                    values = lines[i + 1].split()
                    if len(values) != len(header[1:].split()):
                        print('ERROR: line{0}: {1} needs {2} parameters but only {3} input found'
                              .format(i + 1, header, len(header[1:].split()), len(values)))
                        exit(-1)
                    keys = header[1:].split()
                    values = [float(v) if '.' in v else int(v) for v in values]
                    tmp_dict = dict(zip(keys, values))
                    for k, v in tmp_dict.items():
                        if k in ("lloc", "lpopt"):
                            tmp_dict[k] = int(v)
                        else:
                            tmp_dict[k] = float(v)
                    self._lloc = tmp_dict
                    break

            header = "# l nproj debl"
            for i, line in enumerate(lines):
                if re.match(r'^.*l.*nproj.*', line):
                    tmp_list = []
                    for j in range(i + 1, len(lines)):
                        if lines[j].startswith('#'):
                            break
                        substrings = lines[j].split()
                        values = [float(substring) if '.' in substring else int(substring) for substring in substrings]
                        if len(values) != len(header[1:].split()):
                            print('ERROR: line{0}: {1} needs {2} parameters but only {3} input found'
                                  .format(i + 1, header, len(header[1:].split()), len(values)))
                            exit(-1)
                        tmp_list.append(Lproj(*values))
                    self._projectors = tmp_list
                    break

            header = "# icmod fcfact rcfact"
            for i, line in enumerate(lines):
                if re.match(r'^.*icmod.*', line):
                    values = lines[i + 1].split()
                    values = [float(v) if '.' in v else int(v) for v in values]
                    icmod = int(values[0])
                    if 0 <= icmod <= 2:
                        if len(values) < 2 or len(values) > 3:
                            print('ERROR: for icmod = {0}, at least 2 or 3 parameters are required'
                                  .format(icmod))
                            exit(-1)
                        tmp_dict = {'icmod': icmod, 'fcfact': values[1], 'rcfact': 0.0}
                    elif icmod == 3:
                        if len(values) != len(header[1:].split()):
                            print('ERROR: line{0}: {1} needs {2} parameters but only {3} input found'
                                  .format(i + 1, header, len(header[1:].split()), len(values)))
                            exit(-1)
                        keys = header[1:].split()
                        tmp_dict = dict(zip(keys, values))
                    elif icmod == 4:
                        tmp_dict = {'icmod': 4, 'fcfact': 0.0, 'rcfact': 0.0}
                    else:
                        print('ERROR: icmod should be in 0, 1, 2, 3, 4.')
                        exit(-1)
                    self._nlcc = tmp_dict
                    break

            header = "# epsh1 epsh2 depsh"
            for i, line in enumerate(lines):
                if re.match(r'^.*epsh.*', line):
                    values = [float(x) for x in lines[i + 1].split()]
                    if len(values) != len(header[1:].split()):
                        print('ERROR: line{0}: {1} needs {2} parameters but only {3} input found'
                              .format(i + 1, header, len(header[1:].split()), len(values)))
                        exit(-1)
                    keys = header[1:].split()
                    tmp_dict = dict(zip(keys, values))
                    self._e_range = tmp_dict
                    break

            header = "# rlmax drl"
            for i, line in enumerate(lines):
                if re.match(r'^.*rlmax.*', line):
                    values = [float(x) for x in lines[i + 1].split()]
                    if len(values) != len(header[1:].split()):
                        print('ERROR: line{0}: {1} needs {2} parameters but only {3} input found'
                              .format(i + 1, header, len(header[1:].split()), len(values)))
                        exit(-1)
                    keys = header[1:].split()
                    tmp_dict = dict(zip(keys, values))
                    self._r_range = tmp_dict
                    break

            self._test_configs = 0  # TO DO in the future

            expected_attributes = ["_basic_setting", "_e_config", "_lmax", "_lpsp", "_lloc", "_projectors", "_nlcc",
                                   "_e_range", "_r_range", "_test_configs"]
            if not all(hasattr(self, attr) for attr in expected_attributes):
                raise ValueError("Incomplete initialization: The class should have all required attributes.")

        elif Z:
            nc, nv = count_nc_nv(Z)
            self._basic_setting = {'atsym': atom_name(Z), 'z': Z, 'nc': nc, 'nv': nv, 'iexc': 4, 'psfile': 'both'}

            electron_config = tuple(electron_configuration(Z))
            ll = find_max_l_in_valence(nc, electron_config)
            self._lmax = ll
            if ll + 1 > nv:
                # such like Cu: [Ar]4s1 3d10, ll=lmax=2 and p is empty, which leads error in l from 0 to lmax in lpsp
                self._basic_setting['nc'] = nc - 2
                self._basic_setting['nv'] = nv + 2

            electron_config = sort_nlf_states(electron_config)
            self._e_config = nlf2num(electron_config)

            self._lpsp = [Lpsp(i) for i in range(ll + 1)]

            self._lloc = {'lloc': 4, 'lpopt': 5, 'rc5': 3.0, 'dvloc0': 0.0}

            self._projectors = [Lproj(i) for i in range(ll + 1)]

            self._nlcc = {'icmod': 0, 'fcfact': 0.25, 'rcfact': 0.0}

            self._e_range = {'epsh1': -12.0, 'epsh2': 12.0, 'depsh': 0.02}

            self._r_range = {'rlmax': 6.0, 'drl': 0.01}

            self._test_configs = 0
        else:
            print("ERROR: either Z or path to input file should not be None.")
            exit(-1)

    @property
    def basic_setting(self):
        return self._basic_setting

    @basic_setting.setter
    def basic_setting(self, values):
        for key, value in values.items():
            if key in self._basic_setting:
                self._basic_setting[key] = value
            else:
                raise ValueError(f"Key '{key}' does not exist in basic_setting")

    @property
    def e_config(self):
        return self._e_config

    @e_config.setter
    def e_config(self, values):
        for key, value in values.items():
            if key in self._e_config:
                self._e_config[key] = value
            else:
                raise ValueError(f"Key '{key}' does not exist in e_config")

    @property
    def lmax(self):
        return self._lmax

    @lmax.setter
    def lmax(self, value):
        self._lmax = value

    @property
    def lpsp(self):
        return self._lpsp

    @lpsp.setter
    def lpsp(self, new_lpsp_list):
        if not isinstance(new_lpsp_list, list):
            raise ValueError("lpsp must be a list of Lpsp objects")

        if not all(isinstance(obj, Lpsp) for obj in new_lpsp_list):
            raise ValueError("All elements in the lpsp list must be instances of Lpsp")

        self._lpsp = new_lpsp_list

    @property
    def lloc(self):
        return self._lloc

    @lloc.setter
    def lloc(self, values):
        for key, value in values.items():
            if key in self._lloc:
                self._lloc[key] = value
            else:
                raise ValueError(f"Key '{key}' does not exist in lloc")

    @property
    def projectors(self):
        return self._projectors

    @projectors.setter
    def projectors(self, new_lproj_list):
        if not isinstance(new_lproj_list, list):
            raise ValueError("projectors must be a list of Lproj objects")

        if not all(isinstance(obj, Lproj) for obj in new_lproj_list):
            raise ValueError("All elements in the projectors list must be instances of Lproj")

        self._projectors = new_lproj_list

    @property
    def nlcc(self):
        return self._nlcc

    @nlcc.setter
    def nlcc(self, values):
        for key, value in values.items():
            if key in self._nlcc:
                self._nlcc[key] = value
            else:
                raise ValueError(f"Key '{key}' does not exist in nlcc")

    @property
    def e_range(self):
        return self._e_range

    @e_range.setter
    def e_range(self, values):
        for key, value in values.items():
            if key in self._e_range:
                self._e_range[key] = value
            else:
                raise ValueError(f"Key '{key}' does not exist in e_range")

    @property
    def r_range(self):
        return self._r_range

    @r_range.setter
    def r_range(self, values):
        for key, value in values.items():
            if key in self._r_range:
                self._r_range[key] = value
            else:
                raise ValueError(f"Key '{key}' does not exist in r_range")

    @property
    def test_configs(self):
        return self._test_configs

    @test_configs.setter
    def test_configs(self, value):
        self._test_configs = value

    def write_to_path(self, path=None):
        string = ''
        header = '# atsym z nc nv iexc psfile'
        values = header[1:].split()
        string += header + '\n'
        for value in values:
            string += str(self.basic_setting[value]) + ' '
        string += '\n'

        header = "# n l f"
        string += header + '\n'
        for nlf in self.e_config:
            string += '{0} {1} {2}\n'.format(nlf[0], nlf[1], nlf[2])

        string += "# lmax\n" + str(self.lmax) + '\n'

        header = "# l rc ep ncon nbas qcut"
        string += header + '\n'
        for chanel in self.lpsp:
            string += str(chanel) + '\n'

        header = "# lloc lpopt rc5 dvloc0"
        string += header + '\n'
        values = header[1:].split()
        for value in values:
            string += str(self.lloc[value]) + ' '
        string += '\n'

        header = "# l nproj debl"
        string += header + '\n'
        for proj in self.projectors:
            string += str(proj) + '\n'

        header = "# icmod fcfact rcfact"
        string += header + '\n'
        values = header[1:].split()
        for value in values:
            string += str(self.nlcc[value]) + ' '
        string += '\n'

        header = "# epsh1 epsh2 depsh"
        string += header + '\n'
        values = header[1:].split()
        for value in values:
            string += str(self.e_range[value]) + ' '
        string += '\n'

        header = "# rlmax drl"
        string += header + '\n'
        values = header[1:].split()
        for value in values:
            string += str(self.r_range[value]) + ' '
        string += '\n'

        string += str(self.test_configs)

        if path is None:
            z = self.basic_setting['z']
            file_name = ATOMIC_SYMBOLS[z - 1] + '.in'
        elif os.path.isdir(path):
            z = self.basic_setting['z']
            file_name = os.path.join(path, ATOMIC_SYMBOLS[z - 1] + '.in')
        else:
            file_name = path
        with open(file_name, 'w') as file:
            file.write(string)

    def copy(self):
        from copy import deepcopy
        return deepcopy(self)

    @property
    def valence_electron_numbers(self):
        number = 0
        nv = int(self._basic_setting['nv'])
        valance_config = self._e_config[-nv:]
        for config in valance_config:
            number += config[-1]
        return number
