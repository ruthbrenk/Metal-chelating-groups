## Execute at main metal folder, the one containing pdb files folder
## This script will search for residues and atoms near 2.6 Angstrom from the selected metal
## It will not consider the mentioned molecules cited in "remove"
## Creates 1 folder (Near_timestamp) and 2 subfolders (lig_files & lig_coord_2.6A)
## Lig_files folder must be empty to generate sdfs!

### VERSION 07OCT23 BSA
### SAVES COORDINATES PER METAL OCCURRENCE

from datetime import datetime
import Bio.PDB
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser  # Creates a PDBParser object
from collections import Counter
from glob import glob
import os, os.path, sys, string
import numpy

# datetime object containing current date and time
begin = datetime.now()
dt1_string = begin.strftime('%d/%m/%Y %H:%M:%S') # dd/mm/YY H:M:S
timestamp = begin.strftime('%Y%m%d_%H%M')
print('Timestamp:', dt1_string, '\n-------------------------------\n')

p = PDBParser()
cwd = os.getcwd()
pdb_files = glob(f'{cwd}/files/*.pdb')
mf = 'near' # + '_' + timestamp

METAL = sys.argv[1]

print('There are currently', len(pdb_files), 'pdb files for this metal\n')
print(f'All compounds with heteroatoms near 2.6 A of the {METAL} will be considered')
print('Atoms C & H do not count')
print('If P near, ligand is disregarded!')

coordinates, n_value, num_dict = {}, {}, {}

# Compounds to ignore
remove = ['MN', 'MG', 'ZN', 'CL', 'CO', 'F', 'K', 'IR', 'CA', 'NA', 'OH', 'TL', 'SR', 'CD', 'FE', 'PB', 'RB', 'ALF', 'BEF', 'AF3', 'PO4', 'VO4', 'SO4', 'PO3', 'NO3', 'MGF', 'MF4', 'BF2', 'CAC', 'DOD', 'FE2', 'IRI', 'MOO', 'NCO', 'SFL', 'VN3', 'VN4', 'WO4', 'WO6', 'KQB', 'AZI', 'CMO', 'IOD', 'MH2', 'O', 'PEO', 'THJ', '8M0', 'UNX', 'LI', 'CS', 'NH4', '0BE', 'SX', 'SO3', 'FES', 'PI', 'CU', 'UNL', 'UNK', 'BR', 'H2S', 'NI', 'OXY', 'HOA'] # Mostly metals! CAC, CMO have carbons. UNX is an unknown atom or ion. UNL is an unknown ligand. UNK is an unknown amino acid.

remove_atoms = ['C', 'H']
all, removal = [], []
null_pdb = []
only = {}
all_metals, metal_per_pdb, metal_with_near, metal_wo_near = [], [], [], []
metals_same_lig = []
sdf_done = []

for fileName in sorted(pdb_files):
    structure_id = fileName.rsplit('/', 1)[1][:-4]
    structure = p.get_structure(structure_id, fileName)  # Specify PDB file
    model = structure[0]  # 'structures' may contain several proteins in this case only one
    chains = model.get_chains()
    print('\n----------------------------------------------------------------------------\n' + '<<<' + str(structure) + '>>>')
    n_metals, metal_wnear, metal_nonear = 0, 0, 0
    m_atom = []
    print_no_nb = False
    temp_lig = set()
    temp_sdf = {}

    # Print all hetero residues in chain
    def find_het():
        h_res_list = []
        for residue in model.get_residues():
            residue_id = residue.get_id()
            hetfield = residue_id[0]
            if hetfield[0] == 'H':
                h_res_list.append(residue)
        return (h_res_list)

    het_res_list = find_het() 

    for hs in het_res_list: ### only considers hetero residues
        hs_name = hs.get_resname().strip()
        if hs_name == METAL:
            n_metals += 1
            id = str(METAL) + str('{:02d}'.format(n_metals))
            for atom in hs.get_atoms():
                m_atom.append([atom, id])
            if f'{structure_id}_{id}' not in all_metals:
                all_metals.append(f'{structure_id}_{id}')

    print(f'\n {n_metals} {METAL} found in {structure_id}')

    # Find the neighbor hetero residues and atoms
    for atom in m_atom:
        id = atom[1]
        center = atom[0].get_coord()
        atoms = Bio.PDB.Selection.unfold_entities(het_res_list, 'A')
        ns = Bio.PDB.NeighborSearch(atoms)
        neighbors = ns.search(center, 2.6)  # 5.0 for distance in angstrom ## to include all 2.5xx
        #near_same_lig = False

        neighbors.remove(atom[0])

        if not neighbors:
            if structure_id not in null_pdb:
                null_pdb.append(structure_id)
            if int(id[2:]) < 2 and not print_no_nb:
                print('\nNo neighbors found')
                print_no_nb = True

            metal_nonear += 1

        residue_list = Bio.PDB.Selection.unfold_entities(neighbors, 'R')
        res_atom_dict = {}

        if neighbors:
            if len(neighbors) == 1 and neighbors[0] == atom[0]: ## If the only neighbor is the own metal ##### DO: REMOVE THIS???
                if not print_no_nb:
                    print('\nNo neighbors found')
                    print_no_nb = True
                if structure_id not in null_pdb:
                    null_pdb.append(structure_id)

            metal_wnear += 1

            for found_atom in neighbors:
                found_res = found_atom.parent
                res_atom_dict.setdefault(found_res, []).append(found_atom)

        for found_res in residue_list:
            count_found_atom = 0
            n_value.clear()
            num_dict.clear()
            make_sdf, made_sdf = False, False
            count2 = 0

            hres2 = found_res.resname
            hres1 = hres2.strip()
            chainx = res_atom_dict[found_res][0].full_id[2]

            # Separate coordinates for ligand interacting with > 1 metal
            if f'{found_res}_{chainx}' in temp_lig:
                 metals_same_lig.append(f'{structure_id}_{id}')
                 print_later = [f'! SAME LIGAND {structure_id}, {chainx}, {hres1}, {id}']
            else:
                 temp_lig.add(f'{found_res}_{chainx}')

            list_atoms = sorted(res_atom_dict[found_res])

            print(f'\n{id}', 'neighbor:', hres1, chainx, *list_atoms)

            if any(i == hres1 for i in remove): # Exclude the unwanted molecules
                if hres1 not in removal:
                    removal.append(hres1)
                print(f'* {hres1} removed *')
                try:
                    del(print_later)
                except:
                    pass
                continue

            else:
                if hres1 not in all:
                    all.append(hres1)

            for found_atom in list_atoms:
                found_atom1 = found_atom.name
                chain1 = found_atom.full_id[2]
                lig_id = found_atom.full_id[3][1]
                save_id = f'{structure_id}_{hres1}_{lig_id}_{id}'
                tu_id = (structure_id, hres1, lig_id, id)

                a_coord = found_atom.coord

                #na2 = ns.search(a_coord, 2, level='A')  # O-P single = 1.63 A, double = 1.38 A

                if any(i in found_atom1[0] for i in remove_atoms) and 'CL' not in found_atom1: # Make sure it's not a C, H or P ### try to see if it process CL01
                    print(f'* {hres1} has {found_atom1[0]} near the {METAL} *')
                    continue

                else:
                    count_found_atom += 1
                    value = str(a_coord[0]) + '  ' + str(a_coord[1]) + '  ' + str(a_coord[2])
                    n_value[count_found_atom] = value

                    num_c = numpy.array((str(a_coord[0]), str(a_coord[1]), str(a_coord[2])), dtype=float)
                    num_dict.setdefault(found_atom1, []).append(num_c) # 1st item is the coordinates
                    distance_metal = numpy.linalg.norm(center - num_c)
                    num_dict.setdefault(found_atom1, []).append(str(distance_metal)) # 2nd item is the distance to the metal

                for v in list(n_value.values()):
                    coordinates.setdefault(save_id, [])
                    if v not in coordinates[save_id]:
                        coordinates[save_id].append(v)
                        make_sdf = True

            #if save_id in sorted(coordinates.keys(), reverse = True): ## maybe it takes too long?
            if make_sdf:
                # Creates sdf files for ligands near wanted metal
                os.makedirs(f'{cwd}/{mf}/lig_files', exist_ok=True)
                sdf_check = hres1 + '_' + str(list_atoms)
                temp_sdf.setdefault(hres1, set())
                #if sdf_check in temp_sdf: ### If set of atoms is NOT unique
                if sdf_check in temp_sdf[hres1]: ### If set of atoms is NOT unique
                    try:
                        note2 = '^ Set of near atoms already save as sdf for this ligand ^'
                        if note2 not in print_later:
                            print_later.append(note2)
                    except:
                        print_later = [note2]

                #if sdf_check not in temp_sdf: ### Only if set of atoms is unique
                #if sdf_check not in reversed(temp_sdf): ### If set of atoms is NOT unique
                else: ### If set of atoms is NOT unique
                    if save_id not in sdf_done: ### Only if not made yet
                        with open(f'{cwd}/{mf}/lig_files/prova_pymol_{hres1}_{lig_id}_{chain1}.py', 'w') as f:
                            f.write(f'cmd.load("{cwd}/files/{structure_id}.pdb", object="{save_id}")\n')
                            f.write(f'cmd.remove(selection="resn hoh")\n')
                            f.write(f'cmd.select(name="lig", selection="resi {lig_id} in chain {chain1}")\n')
                            f.write(f'cmd.save(filename="{cwd}/{mf}/lig_files/{hres1}_{lig_id}_{chain1}_raw.sdf", selection="lig")\n')
                            f.write('cmd.reinitialize()')
                        os.system(f'pymol -q -c {cwd}/{mf}/lig_files/prova_pymol_{hres1}_{lig_id}_{chain1}.py')
                        with open(f'{cwd}/{mf}/lig_files/{hres1}_{lig_id}_{chain1}_raw.sdf', 'r') as raw:
                            with open(f'{cwd}/{mf}/lig_files/{save_id}.sdf', 'w') as final:
                                rl = raw.readlines()
                                print('## Created ' + rl[0][:-1] + '.sdf')
                                made_sdf = True
                                for i in rl:
                                    if i != '>  <pdb_header>\n':
                                        final.write(i)
                                    else:
                                        break

                        os.system(f'rm {cwd}/{mf}/lig_files/prova_pymol_{hres1}_{lig_id}_{chain1}.py {cwd}/{mf}/lig_files/{hres1}_{lig_id}_{chain1}_raw.sdf')

                        sdf_done.append(save_id)

                    #temp_sdf.add(sdf_check)
                    temp_sdf[hres1].add(sdf_check)

            for key, value in sorted(num_dict.items(), key=lambda e: e[1][1]):
                count2 += 1
                if count_found_atom >= 1 and made_sdf:
                    print('\n' + str(count2) + ' x ' + structure_id + ', Chain: ' + chain1 + ', Residue: ' + hres1 + ', Atom: ' + key + ' >neighbor atom<')
                    print('    ' + str(value[0]))
                    print('    ' + 'Distance to the metal: ' + str(round(float(value[1]), 4)))

                for n in range(1, count_found_atom + 1):
                    only_key = f'only_{n}'
                    only.setdefault(only_key, [])
                    if hres1 not in only[only_key]:
                        if hres1 not in remove:
                            only[only_key].append(hres1)

            try:
                for note in print_later:
                    print(note)
                del(print_later)
            except:
                pass

    metal_per_pdb.append(n_metals)
    metal_with_near.append(metal_wnear)
    metal_wo_near.append(metal_nonear)

    # Make files with coordinates
    for key in list(coordinates.keys()): ## key = {structure_id}_{hres1}_{lig_id}_{id}
        os.makedirs(f'{cwd}/{mf}/lig_coord_2.6A', exist_ok=True)
        if key in sdf_done:
        #if key == key: ## see the difference
            with open(f'{cwd}/{mf}/lig_coord_2.6A/{key}_2.6A_near_atoms.txt', 'w+') as exit_file:
                exit_file.write('\n'.join(map(str, coordinates[key])) + '\n')

with open(f'{cwd}/{mf}/ligands_only1.txt', 'w+') as st_file:
    for i in sorted(only['only_1']):
        st_file.write('%s\n' % i)

lig_done = []


for i in sdf_done:
    hres1 = i.split('_')[1]
    if hres1 not in lig_done:
        lig_done.append(hres1)

with open(f'{cwd}/{mf}/pdbs_without_near.txt', 'w') as out:
    out.write(str(sorted(null_pdb)))

print('\n                    ┌───────────────┐')
print('                    │ End of folder │')
print('                    └───────────────┘\n')

all2 = all.copy()

print('┌────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐')
print('                                    Overview')
print('\t', str(len(all2)), 'ligands seen.')

for i in lig_done:
    if i in all2:
        all2.remove(i)

actual_ligs = len(all) - len(all2)

#print('\t', str(actual_ligs), 'did not containing P.', str(lig_count2), 'had P', 'and', str(len(repeat)), 'in both.')
print('\t', str(len(all2)) + ' ligands had C or H near the metal: ' + str(all2))
print('\t', str(actual_ligs) + ' ligands in total.')
print('└────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘')

print('\n---------------------------------------------------------------------------------------------------')

mpdb = Counter(metal_per_pdb)
mnatm = Counter(metal_with_near)
mwonatm = Counter(metal_wo_near)

print(f'\n## Number of {METAL} per structure:', mpdb, '\n')
print(f'\n## Number of {METAL} WITH NEAR ATOMS per structure:', mnatm, '\n')
print(f'\n## Number of {METAL} WITHOUT NEAR ATOMS per structure:', mwonatm, '\n')

print('---------------------------------------------------------------------------------------------------')

percentage = round((int(len(null_pdb)) * 100) / int(len(pdb_files)), 2)
print(f'\n## Structures without heteroatoms within 2.6 A of the {METAL} (TOTAL: ' + str(len(null_pdb)) + '/'+ str(len(pdb_files)) + ', ' + str(percentage) + f'%). Saved as {cwd}/{mf}/pdbs_without_near.txt\n')

print('---------------------------------------------------------------------------------------------------')

print('\n## All ligands seen (TOTAL: ' + str(len(all)) + '): ' + str(sorted(all)) + '\n')
print('\n\t## Ligands REMOVED (TOTAL: ' + str(len(removal)) + '): ' + str(sorted(removal)) + '\n')
print('\n## All ligands with heteroatom within 2.6 A (TOTAL: ' + str(len(all2)) + '): ' + str(sorted(all2)) + '\n')

print('---------------------------------------------------------------------------------------------------\n')

print('## Created coordinates files for the ligands (TOTAL: ' + str(len((sdf_done))) + ').')

print('\n---------------------------------------------------------------------------------------------------\n')

print('## Ligands per number of near atoms')
for key in only.keys():
    number = key.split('_')[1]
    if number == '1':
        print('\n' + str(len(only[key])) + ' ligand(s) with ' + number + ' near atom:\t' + str(sorted(only[key])))
    else:
        print('\n' + str(len(only[key])) + ' ligand(s) with ' + number + ' near atoms:\t' + str(sorted(only[key])))

print('\n---------------------------------------------------------------------------------------------------\n')

metals_same_lig = list(set(metals_same_lig)) # makes items unique
all_metals = list(set(all_metals))

print(str(len(metals_same_lig)) + '/' + str(len(all_metals)) + ' ' + METAL + 's for the same ligand.')

end = datetime.now()
dt2_string = end.strftime('%d/%m/%Y %H:%M:%S')
print('\n-------------------------------\nTimestamp:', dt2_string)