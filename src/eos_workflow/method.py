from typing import Optional, Iterable, Iterator
from bisect import bisect_left

INERTS_CONFIG_INDEX = (1, 3, 5, 8, 11, 15, 19)
INERTS_NUMBER = (2, 10, 18, 36, 54, 86, 118)
SUBSHELL_VALUE = (2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 10, 6)
SUBSHELL_NAME = ("1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p", "5s", "4d", "5p",
                 "6s", "4f", "5d", "6p", "7s", "5f", "6d", "7p")
INERT_ELEMS = ('He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn')
L_MAPPING = {'s': 0, 'p': 1, 'd': 2, 'f': 3}

ATOMIC_SYMBOLS = (
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
)

CONFIGS_OF_ODDS = {
    24: ('Ar', ("4s1", "3d5")),
    29: ('Ar', ("4s1", "3d10")),
    41: ('Kr', ("5s1", "4d4")),
    42: ('Kr', ("5s1", "4d5")),
    44: ('Kr', ("5s1", "4d7")),
    45: ('Kr', ("5s1", "4d8")),
    46: ('Kr', ("4d10",)),
    47: ('Kr', ("5s1", "4d10")),
    57: ('Xe', ("6s2", "5d1")),
    58: ('Xe', ("6s2", "4f1", "5d1")),
    64: ('Xe', ("6s2", "4f7", "5d1")),
    78: ('Xe', ("6s1", "4f14", "5d9")),
    79: ('Xe', ("6s1", "4f14", "5d10")),
    89: ('Rn', ("7s2", "6d1")),
    90: ('Rn', ("7s2", "6d2")),
    91: ('Rn', ("7s2", "5f2", "6d1")),
    92: ('Rn', ("7s2", "5f3", "6d1")),
    93: ('Rn', ("7s2", "5f4", "6d1")),
    96: ('Rn', ("7s2", "5f7", "6d1")),

    63: ('Xe', ("6s2", "4f6", "5d1")),
    65: ('Xe', ("6s2", "4f8", "5d1")),
    66: ('Xe', ("6s2", "4f9", "5d1")),
    70: ('Xe', ("6s2", "4f13", "5d1")),
}


def generate_electrons(index: int, copy: int) -> Iterator[str]:
    for base in SUBSHELL_VALUE[index:]:
        if base >= copy:
            if copy > 0:
                yield SUBSHELL_NAME[index] + str(copy)
            break
        yield SUBSHELL_NAME[index] + str(base)
        index += 1
        copy -= base


def make_compact_config(electrons: int) -> tuple[
    Optional[str],  # inert head
    Iterable[str],  # electron config nls
]:
    odd_config = CONFIGS_OF_ODDS.get(electrons)
    if odd_config is not None:
        return odd_config

    if electrons <= 2:
        return None, generate_electrons(index=0, copy=electrons)

    inert_index = bisect_left(INERTS_NUMBER, electrons) - 1
    inert_number = INERTS_NUMBER[inert_index]
    subshell_index = INERTS_CONFIG_INDEX[inert_index]
    inert_head = ATOMIC_SYMBOLS[inert_number - 1]
    return inert_head, generate_electrons(index=subshell_index, copy=electrons - inert_number)


def make_full_config(electrons: int) -> Iterator[int]:
    inert_elem, compact_config = make_compact_config(electrons)

    if electrons > 2:
        inert_electrons = ATOMIC_SYMBOLS.index(inert_elem) + 1
        if inert_elem in INERT_ELEMS:
            yield from make_full_config(inert_electrons)

    yield from compact_config


def electron_configuration(Z, compact=False):
    if compact:
        inert_head, config = make_compact_config(Z)
        return inert_head, tuple(config)
    else:
        return tuple(make_full_config(Z))


def atom_name(Z: int) -> str:
    if 1 <= Z <= 118:
        return ATOMIC_SYMBOLS[Z - 1]
    else:
        print('ERROR: the atomic number is out of range [1, 118]')
        exit()


def count_nc_nv(Z: int):
    inert_head, config = make_compact_config(Z)
    try:
        nc = INERTS_CONFIG_INDEX[INERT_ELEMS.index(inert_head)]
    except ValueError:
        nc = 0
    nv = len(tuple(config))
    return nc, nv


def extract_n_and_l(nlf_state: str) -> tuple[int, int]:
    """Extracts the values of n and l from the NLS state string."""
    n = int(nlf_state[0])
    ll = L_MAPPING.get(nlf_state[1].lower())
    return n, ll


def sort_nlf_states(nlf_states) -> list[str]:
    """Sorts the list of NLF states based on the values of n and l."""

    # Define the custom sorting key function
    def custom_sort_key(nlf_state: str) -> tuple[int, int]:
        return extract_n_and_l(nlf_state)

    # Sort the list of NLS states using the custom sorting key
    sorted_nlf_states = sorted(nlf_states, key=custom_sort_key)

    return sorted_nlf_states


def nlf2num(nlf_states) -> list:
    nlf_num = []
    for nlf in nlf_states:
        nlf_num.append([int(nlf[0]), L_MAPPING.get(nlf[1]), float(nlf[2:])])
    return nlf_num


def get_l_num(nlf_state: str) -> int:
    return L_MAPPING.get(nlf_state[1])


def find_max_l_in_valence(nc, econfig_list) -> int:
    l_max = 0
    for i in range(nc, len(econfig_list)):
        nlf = econfig_list[i]
        ll = L_MAPPING.get(nlf[1])
        if ll > l_max:
            l_max = ll
    return l_max


if __name__ == "__main__":
    for i in range(len(ATOMIC_SYMBOLS)):
        z = i+1
        nc, nv = count_nc_nv(z)
        result = tuple(electron_configuration(z))
        ll_max = find_max_l_in_valence(nc, result)
        print('{0}: {1}, lmax: {2}\n'.format(ATOMIC_SYMBOLS[i], result, ll_max))
