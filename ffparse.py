#!/usr/bin/env python3

from molsym import PointGroup

from itertools import accumulate
from operator import mul

import argparse

parser = argparse.ArgumentParser(description='Parse Firefly output file and print XMCQDPT2 results')
parser.add_argument('path', type=str, help='path to Firefly output file')

args = parser.parse_args()

filename = args.path

# class FireflyOutput()

# in state data structure: rel_energy relative to state # 1 ?

# implement reading state tracking
# it seems like:
# 1. xmcqdpt2 lists the CASCI states.
# 2. the CASCI states are in correct order and their CSF printout is correct
# 3. if one wants to find the MCSCF state with the CASCI state number, the ntrack output has to be used
# eg. CASCI state #7 -> MCSCF state #8
# -> this is important when states are to be addressed in the input file (wstate istate etc)

# sind dann nicht mehr nach energie sortiert, sondern nach CASCI reihenfolge der zustände. sinnvoll?
# cas_states[7], cas_states[8] = cas_states[8], cas_states[7]
# eher eine remapping-lookup-tabelle anlegen? könnte in einem FireflyOutput() unsichtbar geregelt werden
# könnte auch einfach als zusätzliche information mit ausgegeben werden, wenn die ci_states abgefragt werden

def prod(lst):
    for value in accumulate(lst, mul):
        pass
    return value

def transition_symmetry(occupations):
    # return prod([pg(mo_symmetries[cas_nmcc + i]) * int(mo_occ) for i, mo_occ in enumerate(occupations)])
    pass

# ##### GLOBAL CONSTANTS #####
a0 = 5.2917720859E-11   # Bohr radius
eV = 0.000123984        # cm-1 to eV; eV = hc 10^2 / e
Eh = 27.2114            # Hartree energy
u = 1.66053892E-27      # atomic mass unit
me = 9.10938291E-31     # electron mass


ff_point_groups = ['C1', 'CS', 'CI', 'CN', 'S2N', 'CNH', 'CNV', 'DN', 'DNH', 'DND', 'T', 'TH', 'TD', 'O', 'OH']

def read_states(f, block_end, state_end=None):
    # parses:

    # STATE #    1  ENERGY =    -626.229533327

    #      CSF      COEF    OCCUPANCY (IGNORING CORE)
    #      ---      ----    --------- --------- -----
    #        1    0.942146  222000000
    #      184   -0.106500  121001001
    #
    # etc ....

    # f file handler
    # single state is delimited by newline, or optionally state_end (str or list of str)
    # states block is delimited by block_end (str or list of str)
    # should be seek'ed to beginning of states list, else will parse entire file
    # returns dict states = {1: {'energy': -555, 'symmetry': 'b2u', 'configurations': [[1, 0.941949, '222000000'], ...]}}

    if state_end is not None and isinstance(state_end, str):
        state_end = list((state_end,))

    if isinstance(block_end, str):
        block_end = list((block_end,))

    states = dict()
    for line in f:
        if any(substr in line for substr in block_end):
            break

        if 'STATE #' in line:
            # state energy
            split_line = line.strip().split('=')
            state_energy = split_line.pop()
            state_energy = float(state_energy.strip())

            # state number
            split_line = split_line.pop()
            split_line = split_line.split('#').pop()
            state_number = int(split_line.split().pop(0))

            # skip empty lines and header
            next(f)
            next(f)
            next(f)

            # read state configurations
            state_configurations = list()
            for line in iter(f.readline, '\n'):
                if state_end is not None and any(substr in line for substr in state_end):
                    # end of CSF list of a state
                    break

                split_line = line.strip().split()
                state_csf_occupation = split_line.pop()
                state_csf_coeff = float(split_line.pop())
                state_csf = int(split_line.pop())
                state_configurations.append([state_csf, state_csf_coeff, state_csf_occupation])

            state_dict = {'energy': state_energy, 'symmetry': None, 'configurations': state_configurations}
            states[state_number] = state_dict
    return states



with open(filename) as f:
    for line in f:
        # molecular symmetry:
        if 'THE POINT GROUP OF THE MOLECULE' in line:
            # THE POINT GROUP OF THE MOLECULE IS DNH
            # THE ORDER OF THE PRINCIPAL AXIS IS     2
            split_line = line.strip().split()
            point_group = split_line.pop()
            if not point_group in ff_point_groups:
                raise RuntimeError('invalid point group')
            if 'N' in point_group:
                line = f.readline()
                split_line = line.strip().split()
                point_group_naxis = split_line.pop()
                if point_group_naxis.isdigit():
                    point_group_naxis = int(point_group_naxis)
                else:
                    raise RuntimeError('NAXIS is not an integer')
            else:
                point_group_naxis = 0
        elif 'THE POINT GROUP IS' in line:
            # THE POINT GROUP IS DNH, NAXIS= 2, ORDER= 8
            pass

        # molecular orbitals
        if 'SYMMETRIES FOR INITIAL GUESS ORBITALS FOLLOW' in line:
            next(f) # skip one line
            mo_symmetries = list()
            for line in f:
                if 'END OF INITIAL ORBITAL SELECTION' in line:
                    break

                split_line = line.strip().split()
                for mo_str in split_line:
                    _, mo = mo_str.split('=')
                    mo_symmetries.append(mo.lower())

        # number of MOs before CAS, in CAS, etc
        # how to handle DRT ?
        if 'GUGA DISTINCT ROW TABLE' in line:
            for line in f:
                if 'THE MAXIMUM ELECTRON EXCITATION WILL BE' in line:
                    # this should never execute because the for loop
                    # should be stopped below. ie could not find occupations
                    # and an error should be raised
                    raise RuntimeError('MO occupations could not be found')

                if '-CORE-' in line and '-INTERNAL-' in line and '-EXTERNAL-' in line:
                    for line in iter(f.readline, '\n'): # read until empty line
                        split_line = line.strip().split()
                        for var, value in zip(split_line[::2], split_line[1::2]):
                            if 'NDOC' in var:
                                cas_ndoc = int(value)
                            elif 'NMCC' in var:
                                cas_nmcc = int(value)
                            elif 'NVAL' in var:
                                cas_nval = int(value)
                    break

        # CASSCF states
        if 'DAVIDSON METHOD CI-MATRIX DIAGONALIZATION' in line:
            # skip two lines to NUMBER OF STATES REQUESTED
            next(f)
            next(f)
            line = f.readline()
            if 'NUMBER OF STATES REQUESTED' in line:
                split_line = line.strip().split()
                cas_nstates = split_line.pop()
            else:
                 # TODO: error handling
                pass

            # state energies and configurations
            cas_states = read_states(f, block_end='TIMING STATISTICS', state_end='END OF CI-MATRIX DIAGONALIZATION')

        # CAS-CI states
        if '-MCCI- BASED ON OPTIMIZED ORBITALS' in line:
            ci_states = read_states(f, block_end='TIMING STATISTICS')

        # XMC-QDPT2 ENERGIES
        if 'XMC-QDPT2 ENERGIES' in line:
            next(f)
            next(f)
            # {'energy': -626.228462998, 'symmetry': '', 'configurations': None}
            mp2_states = dict()
            for line in f:
                if '-------------' in line:
                    break

                split_line = line.strip().split()
                state_energy = float(split_line.pop())
                state_number = int(split_line.pop(0))
                state_dict = {'energy': state_energy, 'symmetry': '', 'configurations': None}
                mp2_states[state_number] = state_dict

        if 'EIGENVECTORS OF THE EFFECTIVE HAMILTONIAN' in line:
            next(f)

            # loop over each block of eigenvectors
            for line in f:
                if 'EIGENVALUES' in line:
                    break

                # read list of state numbers in this block
                split_line = line.strip().split()
                states_in_block = list(map(int, split_line))

                # skip to eigenvector list
                next(f)
                next(f)
                next(f)

                # read eigenvectors for all states in this block
                state_eigenvectors = [[]] * len(states_in_block)
                for line in iter(f.readline, '\n'):
                    split_line = line.strip().split()
                    state_ev_index = int(split_line.pop(0))
                    state_ev_coeffs = list(map(float, split_line))

                    for idx, state_ev_coeff in enumerate(state_ev_coeffs):
                        if state_eigenvectors[idx]:
                            state_eigenvectors[idx].append([state_ev_index, state_ev_coeff, None])
                        else:
                            state_eigenvectors[idx] = [[state_ev_index, state_ev_coeff, None]]

                # insert eigenvectors into mp2_states dictionary
                for idx, state_number in enumerate(states_in_block):
                    mp2_states[state_number]['configurations'] = state_eigenvectors[idx]


# print('point_group:', point_group)
# print('point_group_naxis:', point_group_naxis)
# print(mo_symmetries)
# print('cas_nstates:', cas_nstates)
# print(cas_states[1])
# print(ci_states[1])
# print(mp2_states[1])


pg = PointGroup(point_group, point_group_naxis)

for k, v in cas_states.items():
    configurations_sorted = sorted(v['configurations'], key=lambda x: abs(x[1]), reverse=True)
    # symmetry of all configurations should be the same, therefore take CSF with largest absolute coeff.
    csf_occupation = configurations_sorted[0][2]  # eg. '222000000' or '202020000'
    transition_sym = prod([pg(mo_symmetries[cas_nmcc + i]) * int(mo_occ) for i, mo_occ in enumerate(csf_occupation)])

    cas_states[k]['configurations'] = configurations_sorted
    cas_states[k]['symmetry'] = str(transition_sym)

for k, v in ci_states.items():
    configurations_sorted = sorted(v['configurations'], key=lambda x: abs(x[1]), reverse=True)
    # symmetry of all configurations should be the same, therefore take CSF with largest absolute coeff.
    csf_occupation = configurations_sorted[0][2]  # eg. '222000000' or '202020000'
    transition_sym = prod([pg(mo_symmetries[cas_nmcc + i]) * int(mo_occ) for i, mo_occ in enumerate(csf_occupation)])

    ci_states[k]['configurations'] = configurations_sorted
    ci_states[k]['symmetry'] = str(transition_sym)

for k, v in mp2_states.items():
    configurations_sorted = sorted(v['configurations'], key=lambda x: abs(x[1]), reverse=True)
    mp2_states[k]['configurations'] = configurations_sorted
    largest_state_contrib_number = configurations_sorted[0][0]
    mp2_states[k]['symmetry'] = ci_states[largest_state_contrib_number]['symmetry']


def get_largest_contributions(configurations, threshold=0.4):
    # returns only those configurations whose coefficient is larger than
    # <threshold> percent of the configuration with the highest coefficient
    # fkt. könnte einen generator mit yield zurückgeben
    configurations_sorted = sorted(configurations, key=lambda x: abs(x[1]), reverse=True)
    threshold_coefficient = abs(configurations_sorted[0][1] * threshold)
    thresholded_configurations = [cfg for cfg in configurations_sorted if abs(cfg[1]) >= threshold_coefficient]
    return thresholded_configurations


# number all MOs
character_counts = {character:mo_symmetries.count(character) for character in mo_symmetries}
numbered_mo_symmetries = list()
n_mos = len(mo_symmetries)
for i,c in enumerate(reversed(mo_symmetries)):
    # index i is also reversed, ie. i = 0 is the last element
    j = n_mos - i - 1
    numbered_mo_symmetries.append((character_counts[c], mo_symmetries[j]))
    character_counts[c] -= 1
numbered_mo_symmetries.reverse()


def parse_occupations(occupations):
    occ_doc = list(map(int, list(occupations[:cas_ndoc])))
    occ_doc_symmetries = mo_symmetries[cas_nmcc:cas_nmcc+cas_ndoc]

    occ_val = list(map(int, list(occupations[cas_ndoc:])))
    occ_val_symmetries = mo_symmetries[cas_nmcc+cas_ndoc:cas_nmcc+cas_ndoc+cas_nval]

    if prod(occ_doc) == 2**cas_ndoc:
        # all DOC are occupied: ground state
        # return transition HOMO -> HOMO as placeholder/indicator
        return list()

    transition_sym = prod([pg(mo_symmetries[cas_nmcc + i]) * int(mo_occ) for i, mo_occ in enumerate(occupations)])

    def parse_doc(doc, val, transitions=[]):
        if sum(val) == 0:
            return transitions

        ndoc = len(doc)
        nval = len(val)
        for i in range(ndoc):
            # loop over all DOC orbitals
            doc_symmetry = mo_symmetries[cas_nmcc+i]
            for j in range(2 - doc[i]):
                # loop over all electrons in DOC orbital
                for k in range(nval):
                    # loop over all VAL orbitals to find the excited electron
                    if val[k] > 0:
                        val_symmetry = mo_symmetries[cas_nmcc+cas_ndoc+k]
                        transitions.append((cas_nmcc+i, cas_nmcc+cas_ndoc+k))
                        all_transition_sym = prod([pg(mo_symmetries[i])*pg(mo_symmetries[f]) for (i,f) in transitions])
                        if all_transition_sym == transition_sym:
                            val[k] -= 1
                            doc[i] += 1
                            return parse_doc(doc, val, transitions)
                        else:
                            new_val = val.copy()
                            new_doc = doc.copy()
                            new_val[k] -= 1
                            new_doc[i] += 1
                            return parse_doc(new_doc, new_val, transitions)

    return parse_doc(occ_doc, occ_val)


str_header = 'XMCQDPT2 RESULTS from ' + filename
print(str_header)
print('=' * len(str_header))
print()
print('ACTIVE SPACE: NMCC:', cas_nmcc, 'NDOC:', cas_ndoc, 'NVAL:', cas_nval)
print('DOC:', ', '.join(['{0[0]} {0[1]}'.format(mo) for mo in numbered_mo_symmetries[cas_nmcc:cas_nmcc+cas_ndoc]]))
print('VAL:', ', '.join(['{0[0]} {0[1]}'.format(mo) for mo in numbered_mo_symmetries[cas_nmcc+cas_ndoc:cas_nmcc+cas_ndoc+cas_nval]]))
# print('DOC:', ', '.join(mo_symmetries[cas_nmcc:cas_nmcc+cas_ndoc]))
# print('VAL:', ', '.join(mo_symmetries[cas_nmcc+cas_ndoc:cas_nmcc+cas_ndoc+cas_nval]))
print()

ref_energy = mp2_states[1]['energy']
for k, v in mp2_states.items():
    print('STATE # {}\t\t{}\t\t{:.3f} eV'.format(k, v['symmetry'], (v['energy']-ref_energy) * Eh))
    configurations = get_largest_contributions(v['configurations'])
    for cfg in configurations:
        ci_state_number = cfg[0]
        print('\t{: 0.6f}\tCI STATE # {}'.format(cfg[1], ci_state_number))
        ci_state_configurations = get_largest_contributions(ci_states[ci_state_number]['configurations'])
        for ci_cfg in ci_state_configurations:
            ci_occ = '{0} {1}'.format(ci_cfg[2][:cas_ndoc], ci_cfg[2][cas_ndoc:])
            mo_transitions = parse_occupations(ci_cfg[2])
            print('\t\t{: 0.6f}\tCSF {:5d}\t({})'.format(ci_cfg[1], ci_cfg[0], ci_occ))
            for (i,f) in mo_transitions:
                print('\t\t\t{0[0]} {0[1]} --> {1[0]} {1[1]}'.format(numbered_mo_symmetries[i], numbered_mo_symmetries[f]))
