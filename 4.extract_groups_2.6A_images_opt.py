## Based on the near atoms, it will try to make a submol
## For only 1 atom, it tries to find a ring
## If it has a ring: if near atom in it, it includes the whole ring: if not, it includes the shortest path to and the ring
## If it doesn't have a ring, makes subgraphs of 1 to 3 atoms.
## For 2 or more atoms, it makes the shortest path in between such atoms.
## It also adds the ring if any of those atoms are in it.
## For all molecules, it makes a subgraph to try to encompass "terminal" atoms. Except those already with subgraphs.
## No stereochemistry in the submols!

## BSA 04NOV23
## 2.6A radius
## Generates images if no sys.argv
## "SMARTS" for SpaceMACS

## It creates a subfolder "submols" and a dictionary of the submols as SMILES

from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from collections import Counter
from cairosvg import svg2png
from datetime import datetime
from glob import glob
import os, sys
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# datetime object containing current date and time
begin = datetime.now()
dt1_string = begin.strftime('%d/%m/%Y %H:%M:%S') # dd/mm/YY H:M:S
timestamp = begin.strftime('%Y%m%d_%H%M')
print('Timestamp:', dt1_string, '\n-------------------------------\n')

############ FROM RDKIT #########################################

### Gets fused OR NOT rings!!!

def GetRingSystems(mol, includeSpiro=False):
    ri = mol.GetRingInfo()
    systems = []
    for ring in ri.AtomRings():
        ringAts = set(ring)
        nSystems = []
        
        for system in systems:
            nInCommon = len(ringAts.intersection(system))
            #if nInCommon and (includeSpiro or nInCommon>1): # FUSED
             #   ringAts = ringAts.union(system)
            #else:
             #   nSystems.append(system)
            nSystems.append(system)  ### NOT FUSED
        nSystems.append(ringAts)
        systems = nSystems
    return systems

###########################################################

def Symbol(index):
    only_1 = False
    if type(index) is not list:
        symbol = newMol.GetAtomWithIdx(index).GetSymbol()
        only_1 = True
        return symbol
    if not only_1: # and len(index) > 1:
        result = []
        for i in index:
            symbol = newMol.GetAtomWithIdx(i).GetSymbol()
            result.append(symbol)
        return result

def ListAtomsInRings(has_rings):
    all_ring_atoms = []
    for ring in list(has_rings):
        for atom in ring:
            if atom not in all_ring_atoms:
                all_ring_atoms.append(atom)
    return all_ring_atoms

def AtomInRing(has_rings, AtomIndex):
    for ring in list(has_rings):
        if AtomIndex in list(ring):
            return list(ring)

def MakeSubgraph(lst_atoms, n_bonds):
    not_heteroatoms = ['C', 'c', 'H', 'h']
    subgraph_dict = {}

    for x in lst_atoms:
        ### Level 1 bond
        subgraph_001 = Chem.rdmolops.FindAllSubgraphsOfLengthN(newMol, 1, rootedAtAtom=x)
        for bonds in subgraph_001:
            for bond in bonds:
                atm1 = newMol.GetBondWithIdx(bond).GetBeginAtomIdx()
                atm2 = newMol.GetBondWithIdx(bond).GetEndAtomIdx()
                if atm1 != x:
                    ### Level 2 bonds
                    subgraph_021 = Chem.rdmolops.FindAllSubgraphsOfLengthN(newMol, 1, rootedAtAtom=atm1)
                    for bonds_021 in subgraph_021:
                        for bond_021 in bonds_021:
                            atm211 = newMol.GetBondWithIdx(bond_021).GetBeginAtomIdx()
                            atm212 = newMol.GetBondWithIdx(bond_021).GetEndAtomIdx()
                            if Symbol(atm211) not in not_heteroatoms:
                                subgraph_dict.setdefault('02', []).append(atm1)
                                subgraph_dict.setdefault('02', []).append(atm211)
                            if Symbol(atm212) not in not_heteroatoms:
                                subgraph_dict.setdefault('02', []).append(atm1)
                                subgraph_dict.setdefault('02', []).append(atm212)

                            ### Level 3 bonds
                            if n_bonds == 3:
                                subgraph_031 = Chem.rdmolops.FindAllSubgraphsOfLengthN(newMol, 1, rootedAtAtom=atm211)
                                for bonds_031 in subgraph_031:
                                    for bond_031 in bonds_031:
                                        atm311 = newMol.GetBondWithIdx(bond_031).GetBeginAtomIdx()
                                        atm312 = newMol.GetBondWithIdx(bond_031).GetEndAtomIdx()
                                        if Symbol(atm311) not in not_heteroatoms:
                                            subgraph_dict.setdefault('03', []).append(atm1)
                                            subgraph_dict.setdefault('03', []).append(atm211)
                                            subgraph_dict.setdefault('03', []).append(atm311)
                                        if Symbol(atm312) not in not_heteroatoms:
                                            subgraph_dict.setdefault('03', []).append(atm1)
                                            subgraph_dict.setdefault('03', []).append(atm211)
                                            subgraph_dict.setdefault('03', []).append(atm312)

                                subgraph_032 = Chem.rdmolops.FindAllSubgraphsOfLengthN(newMol, 1, rootedAtAtom=atm212)
                                for bonds_032 in subgraph_032:
                                    for bond_032 in bonds_032:
                                        atm321 = newMol.GetBondWithIdx(bond_032).GetBeginAtomIdx()
                                        atm322 = newMol.GetBondWithIdx(bond_032).GetEndAtomIdx()
                                        if Symbol(atm321) not in not_heteroatoms:
                                            subgraph_dict.setdefault('03', []).append(atm1)
                                            subgraph_dict.setdefault('03', []).append(atm212)
                                            subgraph_dict.setdefault('03', []).append(atm321)
                                        if Symbol(atm322) not in not_heteroatoms:
                                            subgraph_dict.setdefault('03', []).append(atm1)
                                            subgraph_dict.setdefault('03', []).append(atm212)
                                            subgraph_dict.setdefault('03', []).append(atm322)

                if atm2 != x:
                    ### Level 2 bonds
                    subgraph_022 = Chem.rdmolops.FindAllSubgraphsOfLengthN(newMol, 1, rootedAtAtom=atm2)
                    for bonds22 in subgraph_022:
                        for bond22 in bonds22:
                            atm221 = newMol.GetBondWithIdx(bond22).GetBeginAtomIdx()
                            atm222 = newMol.GetBondWithIdx(bond22).GetEndAtomIdx()
                            if Symbol(atm221) not in not_heteroatoms:
                                subgraph_dict.setdefault('02', []).append(atm2)
                                subgraph_dict.setdefault('02', []).append(atm221)
                            if Symbol(atm222) not in not_heteroatoms:
                                subgraph_dict.setdefault('02', []).append(atm2)
                                subgraph_dict.setdefault('02', []).append(atm222)

                            ### Level 3 bonds
                            if n_bonds == 3:
                                subgraph_033 = Chem.rdmolops.FindAllSubgraphsOfLengthN(newMol, 1, rootedAtAtom=atm221)
                                for bonds_033 in subgraph_033:
                                    for bond_033 in bonds_033:
                                        atm313 = newMol.GetBondWithIdx(bond_033).GetBeginAtomIdx()
                                        atm314 = newMol.GetBondWithIdx(bond_033).GetEndAtomIdx()
                                        if Symbol(atm313) not in not_heteroatoms:
                                            subgraph_dict.setdefault('03', []).append(atm2)
                                            subgraph_dict.setdefault('03', []).append(atm221)
                                            subgraph_dict.setdefault('03', []).append(atm313)
                                        if Symbol(atm314) not in not_heteroatoms:
                                            subgraph_dict.setdefault('03', []).append(atm2)
                                            subgraph_dict.setdefault('03', []).append(atm221)
                                            subgraph_dict.setdefault('03', []).append(atm314)

                                subgraph_034 = Chem.rdmolops.FindAllSubgraphsOfLengthN(newMol, 1, rootedAtAtom=atm222)
                                for bonds_034 in subgraph_034:
                                    for bond_034 in bonds_034:
                                        atm323 = newMol.GetBondWithIdx(bond_034).GetBeginAtomIdx()
                                        atm324 = newMol.GetBondWithIdx(bond_034).GetEndAtomIdx()
                                        if Symbol(atm323) not in not_heteroatoms:
                                            subgraph_dict.setdefault('03', []).append(atm2)
                                            subgraph_dict.setdefault('03', []).append(atm222)
                                            subgraph_dict.setdefault('03', []).append(atm323)
                                        if Symbol(atm324) not in not_heteroatoms:
                                            subgraph_dict.setdefault('03', []).append(atm2)
                                            subgraph_dict.setdefault('03', []).append(atm222)
                                            subgraph_dict.setdefault('03', []).append(atm324)

    for key in subgraph_dict.keys():
        subgraph_dict[key] = list(set(subgraph_dict[key]))
        subgraph_dict[key] = sorted(subgraph_dict[key])


    ## Combine both 02 and 03 bonds subgraphs (just in case)
    if '03' in list(subgraph_dict.keys()):
        combined = subgraph_dict['02'] + subgraph_dict['03']
        combined = list(set(combined))
        combined = sorted(combined)

        if combined != subgraph_dict['03']:
            if combined != subgraph_dict['02']:
                subgraph_dict['02-03'] = combined

    sorted_dict = dict(sorted(subgraph_dict.items()))

    return sorted_dict

def FixColon(frag): ### fix C:C to c:c
    frag_bef = frag
    positions = []
    for pos, ch in enumerate(frag):
        if ch == ':':
            position = pos
            positions.append(position)
    for position in positions:
        before = position - 1
        after = position + 1
        b_extra = before - 1
        b_extra2 = before - 2
        a_extra = after + 1

        if frag[before].isupper():
            if before == 0:
                frag = frag[before].lower() + frag[position:] 
            else:
                frag = frag[:before] + frag[before].lower() + frag[position:]
        if frag[before].isnumeric():
            frag = frag[:b_extra] + frag[b_extra].lower() + frag[before] + frag[position:]
        if frag[after].isupper():
            if after == len(frag) - 1:
                frag = frag[:position] + frag[position] + frag[after].lower()
            else:
                frag = frag[:after] + frag[after].lower() + frag[a_extra:]

    digits = []
    mb_colon = False
    for char in frag:
        if char.isdigit():
            if char not in digits: ## it also gets digits for the charge
                digits.append(char)
            mb_colon = True

    if mb_colon:
        if all(frag.count(digit) % 2 == 0 for digit in digits): ### even numbers
            frag = frag.replace(':', '')

    frag_fixed = f'{frag_bef} >>> {frag}'

    return frag, frag_fixed

def FixHs(frag):
    frag_bef = frag
    fix = {'[H]O': '[OD1]', 'O[H]': '[OD1]', 'N([H])([H])[H]': '[ND1]', '[N+]([H])([H])[H]': '[Nh3+]', '[H][N+]([H])([H])': '[Nh3+]', '[N+]([H])[H]': '[Nh2+]', '[H]N([H])': '[ND1]', 'N([H])[H]': '[ND1]', '[H][n+]': '[nh1+]', '[n+]([H])': '[nh1+]', '[H][N+]': '[Nh1+]', '[n+]2[H]': '[nh1+]2','[H]N': '[ND2]', 'N[H]': '[ND2]', '[n+]1[H]': '[nh1+]1', 'N([H])': '[ND2]', '[H]n': '[nh1]', 'N1[H]': '[ND2]1', 'n:1[H]': '[nh1]1', 'n1[H]': '[nh1]1', 'n([H])': '[nh1]', '[H]S': '[SD1]', 'S[H]': '[SD1]', 'S([H])': '[SD1]','[H][Se]': '[SeD1]', '[Se][H]': '[SeD1]', '[PH]': '[P]', '[SH]': '[SD1]', '[OH+]': '[OD1+]', '[OH2+]': '[OD1+]', '[H][O+]': '[OD1+]', '[O+][H]': '[OD1+]', '[PH+]': '[P+]', '([H])P': 'P', '[B-]([H])([H])[H]': '[BD1-]', '[O+]([H])[H]': '[OD1+]'}
    for ele in fix.keys():
        if ele in frag:
            if ele == 'n:1[H]':
                frag = frag.replace(ele, fix[ele])
                frag = frag.replace(':', '')
            else:
                frag = frag.replace(ele, fix[ele])

    frag_fixed2 = f'{frag_bef} >>> {frag}'
    return frag, frag_fixed2

broken = {}
def BrokenSmiles(name, submols_smiles): ### Compiles submols containing ':' or '.'
    broken_sym = ['.', ':']

    for symbol in broken_sym:
        if symbol in submols_smiles[-1]:
            if not broken:
                broken.setdefault(str(f'{symbol}_last'), []).append(f'{name}_0')
            if f'{symbol}_last' not in broken.keys():
                broken.setdefault(str(f'{symbol}_last'), []).append(f'{name}_0')
            if f'{name}_0' not in list(broken[f'{symbol}_last']):
                broken.setdefault(f'{symbol}_last', []).append(f'{name}_0')

        for idx, submol in enumerate(submols_smiles):
            total = len(submols_smiles) - 1
            if idx < total:
                if symbol in submol:
                    if not broken:
                        broken.setdefault(symbol, []).append(f'{name}_{idx}')
                    if symbol not in broken.keys():
                        broken.setdefault(symbol, []).append(f'{name}_{idx}')
                    if f'{name}_{idx}' not in list(broken[symbol]):
                        broken.setdefault(symbol, []).append(f'{name}_{idx}')

def Percentage(want, total):
    if isinstance(want, list):
        want = int(len(want))
    if isinstance(total, list):
        total = int(len(total))

    percentage = round((want * 100) / total, 2)
    answer = '(' + str(percentage) + '%)'

    return answer

###################### PICS ##########################

def draw_lig(ligand, mol):
    PATH_TO_OUTPUT = f'{mf}/images_{timestamp}/ligands'
    os.makedirs(f'{PATH_TO_OUTPUT}', exist_ok=True)
    Chem.RemoveStereochemistry(mol)
    Chem.RemoveHs(mol, implicitOnly=True)
    Chem.SanitizeMol(mol)

    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(800,600)
    drawer.drawOptions().useBWAtomPalette() #Will draw only in black
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg1 = drawer.GetDrawingText().replace('svg:','')
    with open(f'{PATH_TO_OUTPUT}/{ligand}_bw.svg', 'w+') as outf1:
        outf1.write(svg1)
    svg2png(url=f'{PATH_TO_OUTPUT}/{ligand}_bw.svg', write_to=f'{PATH_TO_OUTPUT}/{ligand}_bw.png')

done_submols = []
def draw_submols(ligand, mol, output, lst_atm_highlight, idx):
    PATH_TO_OUTPUT = f'{mf}/images_{timestamp}/submols'
    os.makedirs(f'{PATH_TO_OUTPUT}', exist_ok=True)
    Chem.RemoveStereochemistry(mol)
    Chem.RemoveHs(mol, implicitOnly=True)
    Chem.SanitizeMol(mol)

    submol = Chem.MolFragmentToSmiles(mol, lst_atm_highlight, isomericSmiles=False) ### no stereochemistry
    if submol not in done_submols:
        done_submols.append(submol)

    submol_n = str('{:04d}'.format(done_submols.index(submol)))

    highlight = lst_atm_highlight
    hit_bonds = []
    
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in highlight:
            if bond.GetEndAtomIdx() in highlight:
                hit_bonds.append(bond.GetIdx())
                
    color_sel1 = {}
    color_sel2 = {}
    color2 = (0.996078431372549, 0.3803921568627451, 0.0) # orange #FE6100

    for i in highlight:
        color_sel1.setdefault(i, color2)

    for i in hit_bonds:
        color_sel2.setdefault(i, color2)
    
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(800,600)
    drawer.drawOptions().useBWAtomPalette() #Will draw only in black
    drawer.DrawMolecule(mol, highlightAtoms=highlight, highlightAtomColors=color_sel1, highlightBonds=hit_bonds, highlightBondColors=color_sel2)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    with open(f'{PATH_TO_OUTPUT}/{output}_{idx+1}_submol-{submol_n}.svg', 'w+') as outf:
        outf.write(svg)
    svg2png(url=f'{PATH_TO_OUTPUT}/{output}_{idx+1}_submol-{submol_n}.svg', write_to=f'{PATH_TO_OUTPUT}/{output}_{idx+1}_submol-{submol_n}.png')
    
def draw_near(ligand, mol, output, indices):
    PATH_TO_OUTPUT = f'{mf}/images_{timestamp}/near'
    os.makedirs(f'{PATH_TO_OUTPUT}', exist_ok=True)
    Chem.RemoveStereochemistry(mol)
    Chem.RemoveHs(mol, implicitOnly=True)
    Chem.SanitizeMol(mol)

    near_atoms = Chem.MolFragmentToSmiles(mol, indices)

    color_sel1 = {}
    color1 = (0.33725490196078434, 0.7058823529411765, 0.9137254901960784) #blue #56B4E9

    for i in indices:
        color_sel1.setdefault(i, color1)
    
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(800,600)
    drawer.drawOptions().useBWAtomPalette() #Will draw only in black
    drawer.DrawMolecule(mol, highlightAtoms=indices, highlightAtomColors=color_sel1)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    with open(f'{PATH_TO_OUTPUT}/{output}_near_{near_atoms}.svg', 'w+') as outf:
        outf.write(svg)
    svg2png(url=f'{PATH_TO_OUTPUT}/{output}_near_{near_atoms}.svg', write_to=f'{PATH_TO_OUTPUT}/{output}_near_{near_atoms}.png')
    
def draw_near_submols(ligand, mol, output, lst_atm_highlight, indices, idx):
    PATH_TO_OUTPUT = f'{mf}/images_{timestamp}/near-submols'
    os.makedirs(f'{PATH_TO_OUTPUT}', exist_ok=True)
    Chem.RemoveStereochemistry(mol)
    #param = Chem.RemoveHsParameters().removeDefiningBondStereo=True
    Chem.RemoveHs(mol, implicitOnly=True)
    Chem.SanitizeMol(mol)

    near_atoms = Chem.MolFragmentToSmiles(mol, indices)

    submol = Chem.MolFragmentToSmiles(mol, lst_atm_highlight, isomericSmiles=False) ### no stereochemistry
    if submol not in done_submols:
        done_submols.append(submol)

    submol_n = str('{:04d}'.format(done_submols.index(submol)))

    hit_ats = indices.copy()
    hit_bonds = []
    
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in lst_atm_highlight:
            if bond.GetEndAtomIdx() in lst_atm_highlight:
                hit_bonds.append(bond.GetIdx())

    color_sel1 = {}
    color_sel2 = {}
    color1 = (0.33725490196078434, 0.7058823529411765, 0.9137254901960784) #blue #56B4E9
    color2 = (0.996078431372549, 0.3803921568627451, 0.0) # orange #FE6100 for the MCS

    for i in indices:
        color_sel1.setdefault(i, color1)
        
    for i in lst_atm_highlight:
        if i not in indices:
            color_sel1.setdefault(i, color2) ## N-H, O-H!
            hit_ats.append(i)
       #     if Symbol(i) == 'H': ## H blue too
        #        color_sel1.setdefault(i, color1)
         #       hit_ats.append(i)
          #  else:
           #     color_sel1.setdefault(i, color2)
            #    hit_ats.append(i)


    for i in hit_bonds:
        color_sel2.setdefault(i, color2)
    
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(800,600)
    drawer.drawOptions().useBWAtomPalette() #Will draw only in black
    drawer.DrawMolecule(mol, highlightAtoms=hit_ats, highlightAtomColors=color_sel1, highlightBonds=hit_bonds, highlightBondColors=color_sel2)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    with open(f'{PATH_TO_OUTPUT}/{output}_{idx+1}_submol-{submol_n}_near_{near_atoms}.svg', 'w+') as outf:
        outf.write(svg)
    svg2png(url=f'{PATH_TO_OUTPUT}/{output}_{idx+1}_submol-{submol_n}_near_{near_atoms}.svg', write_to=f'{PATH_TO_OUTPUT}/{output}_{idx+1}_submol-{submol_n}_near_{near_atoms}.png')
    
#################################### Grouping ###############################################

def Group(initial_dict):
    colon, all_submols, all_submols_og = [], [], []
    d, repeated = {}, {}
    
    for k, v in sorted(initial_dict.items()):
        if int(k) == 1:
            print('\n## From molecules with 1 near atom:', len(v.keys()), sorted(v.keys(), key=len, reverse=False))
        if int(k) > 1:
            print(f'\n## From molecules with {k} near atoms:', len(v.keys()), sorted(v.keys(), key=len, reverse=False))
        for smile, ligs in v.items():
            initial_dict[k][smile] = list(set(initial_dict[k][smile])) ## make it unique
            initial_dict[k][smile] = sorted(initial_dict[k][smile])
            smile_og = smile

            ### FixHs as adequate input for SpaceMACS, blocked near atoms
            if '[H]' in smile:
                fixed = FixHs(smile)
                smile = fixed[0]

            if ':' in smile_og:
                colon.append([smile_og, smile, k])
            else:
                if smile_og not in all_submols_og:
                    all_submols_og.append(smile_og)
                if smile not in all_submols:
                    all_submols.append(smile) ## For SpaceMACS

                ## make dictionary with number of molecules per SMILES
                d.setdefault(smile_og, len(initial_dict[k][smile_og]))
                
    for item in colon:
        smile_og = item[0] ### To print in the script
        smile = item[1] ### If Hs changed to SpaceMACS standard, to save in file
        k = item[2]
        search_og = smile_og.replace(':', '')
        search = smile.replace(':', '')

        if search_og not in all_submols_og:
            all_submols_og.append(smile_og)
        if search not in all_submols:
            all_submols.append(smile) ## For SpaceMACS

            ## make dictionary with number of molecules per SMILES
            d.setdefault(smile_og, len(initial_dict[k][smile_og]))

        ## If repeated, need to change value in d!
        else:
            repeated[search_og] = smile_og

            ## make dictionary with number of molecules per SMILES
            new_value = len(initial_dict[k][smile_og]) + len(initial_dict[k][search_og])
            d[search_og] = new_value
            
    all_submols.sort(key=len)
    all_submols_og.sort(key=len)
    print('\n## From ALL molecules:', len(all_submols_og), all_submols_og)
    if len(repeated) >= 1:
        print("\n## Repeated SMILES due to ':' (.x2.):", len(repeated), end = ' ')
        for k, v in repeated.items():
            if k == list(repeated.keys())[-1]:
                print(f"'{v} >>> {k}'")
            else:
                print(f"'{v} >>> {k}'", end = ', ')
            
    ## Make file with all_submols, give remark
    with open(f'{mf}/submols_{timestamp}/{metal}_all_submols_SMARTS.txt', 'w+') as exit_file:
        exit_file.write(f'{metal}_all_submols = {all_submols}\n')
        exit_file.write(f'#{len(all_submols)}')

    with open(f'{mf}/submols_{timestamp}/{metal}_all_submols_SMILES.txt', 'w+') as exit_file:
        exit_file.write(f'{metal}_all_submols_OG = {all_submols_og}\n')
        exit_file.write(f'#{len(all_submols_og)}')

    ## Sort the number of occurrences from the highest to lowest
    print(f'\n--------------- ** Grouped by number of molecules DESC ** ---------------\n')

    ds = dict(sorted(d.items(), key=lambda x: x[1], reverse=True))
    
    idx = 0
    for k in list(ds.keys()): ## Submol in dict
        idx += 1
        idx2 = str('{:03d}'.format(idx))
        temp = []
        lig_list = []
        print_2 = False
        for length, values in initial_dict.items():
            if k in values.keys():
                temp.append(int(length))
                lig_list.append(initial_dict[length][k])

                if k not in repeated.keys():
                    try:
                        submol_n = str('{:04d}'.format(done_submols.index(k)))
                        part1 = '\t' + str(idx2) + '.\t' + k + ' (submol-' + submol_n + ') '
                    except:
                        part1 = '\t' + str(idx2) + '.\t' + k + ' (submol-xxxx) '
                    part2a = ', ' + str(ds[k]) + ' molecule(s)'
                    print_2 = True

                else:
                    get = repeated[k] ## The one with ':'
                    if int(length) not in temp:
                        temp.append(int(length))
                    part1 = '\t' + str(idx2) + '.x2.\t' + k + ' '
                    part2a = ', ' + str(ds[k]) + ' molecule(s)'
                    try:
                        lig_list.extend(initial_dict[length][get])
                    except:
                        try:
                            for search_length in range(1, int(sorted(initial_dict.keys())[-1])):
                                lig_list.append(initial_dict[search_length][get])
                                break
                        except:
                            pass
                            
                    print_2 = True

        if print_2:
            part2 = str(temp) + part2a
            part3 = str(lig_list)
            space1 = 65 - len(k)
            space2 = 40 - len(part2)
            print(part1 + ' '*space1 + part2 + ' '*space2 + '::\t' + part3)
        else:
            print('ERROR:', k)
            
    return all_submols_og

##################################################################################################

### Files
cwd = os.getcwd()
mf = f'{cwd}/near' # Mother folder
mols = glob(f'{mf}/lig_files/*.sdf')
metal = cwd.rsplit('/', 1)[1]

PATH_TO_CIF_FILES = f'{mf}/cifs'
lig_list = os.listdir(PATH_TO_CIF_FILES)

PATH_TO_COORD = f'{mf}/lig_coord_2.6A'
distance_to_metal1 = PATH_TO_COORD.rsplit('/', 1)[1]
distance_to_metal = distance_to_metal1.rsplit('_')[2][:-1]
lig_w_coord1 = os.listdir(PATH_TO_COORD)
lig_w_coord = []

for file in lig_w_coord1:
    if file[-4:] == '.txt':
        lig = file.split('_', 2)[1]
        lig_w_coord.append(lig)
        
lig_w_coord = list(set(lig_w_coord)) # makes items unique

## Generate pictures
make_pic = True
if len(sys.argv) == 2: ### Add any word after the .py & no picture will be generated
    make_pic = False

print('There are currently', len(lig_list), 'unique ligands for', metal)
print('And', len(lig_w_coord), 'unique ligands have near atoms\n')
print(f'Submols including the atoms near {distance_to_metal} A of the metal will be generated')
print(f'Based on the shortest path between those atoms (2 or more atoms)')
print(f'Or on the presence of a ring (1 atom)')
print(f'No stereochemistry in the submols!')

### Make dictionary with ligands' SMILES
smi_dict = {}
match_smiles1 = ['SMILES_CANONICAL', 'CACTVS', '"']
match_smiles2 = ['SMILES', 'OpenEye OEToolkits', '"']
match_smiles3 = ['SMILES']
match_smiles4 = ['SMILES', 'CACTVS']
for _ in lig_list:
    _ = _[:-4]
    with open(f'{PATH_TO_CIF_FILES}/{_}.cif', 'r') as ref:
        for line in ref:
            if all(x in line for x in match_smiles3):
                smi2 = line.strip() #Remove spaces
                smi1 = list(filter(None, smi2.split(sep=' ')))
                smi = smi1[-1]
                if '"' in smi:
                    smi = smi[1:-1]
                key = _.split(sep="_")[0]
                if key in smi_dict and smi in list(smi_dict[key]):
                    continue
                else:
                    if not '@' in smi:
                        smi_dict.setdefault(key, []).append(smi)

### General
no_smi, no_rings, no_surroundings, empty_path = [], [], [], []
has_near, frag_done, frag_done_lig1 = [], [], []
only_1, more_than_1 = [], []
no_newMol, submol_dict, lig_x_near = {}, {}, {}
smiles_w_hs = {}
lig_pic = []
drew_near, drew_near1, drew_near2 = {}, {}, {}
near_no_h = {}
divergent_n_heavyatms = []
only = {}
last_smi, to_group = {}, {}

### Start
for mol in sorted(mols): ### All SDF files, different PDBs might = different modes of binding
    name = mol.rsplit('/', 1)[1][:-4]
    lig_mol = name.rsplit('_')[1]
    structure_id = name.rsplit('_')[0]
    metal_id = name.rsplit('_')[-1]
    new_name = f'{lig_mol}_{structure_id}_{metal_id}'

    ### Check if it has SMILES
    if lig_mol not in smi_dict:
        if lig_mol not in no_smi:
            no_smi.append(lig_mol)
        print('\n§§§ Current ligand: ' + lig_mol + ' (PDB ID: ' + structure_id + ', sdf: ' + name + ') NOT IN SMI_DICT')

    if lig_mol in smi_dict:
        has_newMol = False
        print('\n¤ Current ligand: ' + lig_mol + ' (PDB ID: ' + structure_id + ', sdf: ' + name + ')')

        ### Check if it has coordinates file with near atoms
        if os.path.isfile(f'{PATH_TO_COORD}/{name}_2.6A_near_atoms.txt'):
            coord_file = f'{PATH_TO_COORD}/{name}_2.6A_near_atoms.txt'

        if coord_file:
            if name not in has_near:
                has_near.append(name)
            mol1 = Chem.MolFromMolFile(mol, removeHs=False, sanitize=False) ### MAYBE CHANGE IT
            smiles_string = smi_dict[lig_mol]
            for idx in enumerate(smiles_string):
                print('SMILES (' + str(idx[0]+1) + '/' + str(len(smiles_string)) + '): ' + idx[1])
                template = Chem.MolFromSmiles(smiles_string[idx[0]])
                try:
                    template_heavyatms = template.GetNumHeavyAtoms()
                    newMol = AllChem.AssignBondOrdersFromTemplate(template, mol1)
                    has_newMol = True
                    mol2d1 = Chem.Mol(newMol)
                    AllChem.Compute2DCoords(mol2d1) ## FROM https://github.com/rdkit/rdkit/discussions/5316
                    break ### ???

                except Exception as e:
                    print('\tERROR:', lig_mol, '- Not possible to create newMol!', e)
                    no_newMol.setdefault(lig_mol, []).append(structure_id)

        ### From the new molecule object
        if has_newMol:
            print_later = []
            atom_list = []
            n_atoms = newMol.GetNumAtoms()
            conf = newMol.GetConformer()

            for i in range(n_atoms):
                atom_list.append([conf.GetAtomPosition(i)[0], conf.GetAtomPosition(i)[1], conf.GetAtomPosition(i)[2]])

            atoms_to_highlight = []

            with open (coord_file, 'r') as f:
                for line in f:
                    atoms_to_highlight.append([float(line.split(sep='  ')[0]), float(line.split(sep='  ')[1]), float(line.split(sep='  ')[2])])
                        
            ## make list with near atoms' indexes
            indices = []
            for i, v in enumerate(atom_list):
                if v in atoms_to_highlight:
                    check_near = Symbol(i)
                    atm_degree = newMol.GetAtomWithIdx(i).GetDegree()
                    if 'H' in check_near: ### Might have slipped away before
                        print(f'* {Symbol(i)} {i} is a hydrogen, excluding it from near atoms!')
                    elif 'P' in check_near: ## P doesn't appear chelating the metals in the PDBs, S does.
                        if atm_degree == 4:
                            print(f'* {Symbol(i)} {i} is a phosphorus (D4), excluding it from near atoms!')
                        else:
                            indices.append(i)
                    elif 'B' in check_near: ## B also doesn't appear chelating the metals in the PDBs.
                        if atm_degree == 4:
                            print(f'* {Symbol(i)} {i} is a boron (D4), excluding it from near atoms!')
                        else:
                            indices.append(i)
                    else:
                        indices.append(i)
     
            if make_pic:
                try:
                    if lig_mol not in lig_pic:
                        draw_lig(lig_mol, mol2d1)
                        lig_pic.append(lig_mol)

                    if not drew_near:
                        draw_near(lig_mol, mol2d1, new_name, indices)
                        drew_near.setdefault(lig_mol, []).append(indices)
                    if lig_mol not in drew_near.keys():
                        draw_near(lig_mol, mol2d1, new_name, indices)
                        drew_near.setdefault(lig_mol, []).append(indices)
                    if indices not in drew_near[lig_mol]:
                        draw_near(lig_mol, mol2d1, new_name, indices)
                        drew_near[lig_mol].append(indices)

                except Exception as e:
                    error1 = f'\tERROR: {lig_mol} - Not possible to draw it! {e}'
                    if error1 not in print_later:
                        print_later.append(error1)

            mol2d = Chem.Mol(newMol)
            AllChem.Compute2DCoords(mol2d) ## FROM https://github.com/rdkit/rdkit/discussions/5316

            try:
                near_atoms = Chem.MolFragmentToSmiles(newMol, indices)
                print(f'Near atoms ({len(indices)}): {near_atoms} {indices}')
            except Exception as e:
                print('\tERROR:', lig_mol, '- Not possible to create indices!', e)
                continue

            length = len(indices)
            has_rings = GetRingSystems(mol2d)
            all_ring_atoms = []
            submols_atm, submols_smiles, submols_smarts =[], [], []
            frag_fixed, frag_fixed2 = [], []
            add_term, add_term1, add_term2, excluded_term = [], [], [], []
            already_subgraph, has_terminal = False, False
            dont_addh = ['C', 'c', 'P']

            ### If it only has 1 near atom
            if length == 1:
                joined = False
                path_to_ring = []
                if lig_mol not in only_1:
                    only_1.append(lig_mol)

                if has_rings:
                    print(lig_mol, 'has rings:', has_rings)
                    all_ring_atoms += ListAtomsInRings(has_rings) ### Make list with all ring atoms
                    atm_in_ring = AtomInRing(has_rings, indices[0]) ### If the near atom is in a ring, gets list with ring atoms

                    ### Check if the near atom is in a ring or not
                    if atm_in_ring:
                        print('\t1.Atom', Symbol(indices[0]), indices[0], 'in a ring:', atm_in_ring, Chem.MolFragmentToSmiles(newMol, list(atm_in_ring), isomericSmiles=False))
                        submols_atm.append(atm_in_ring)

                    else:
                        for ring in list(has_rings):
                            for atom in ring:
                                path_to_ring.append(list(Chem.GetShortestPath(newMol, indices[0], atom)))

                        path_to_ring = [ele for ele in path_to_ring if ele != []] ### remove empty paths
                        path_to_ring.sort(key=len) ### To only get the shortest

                        ### Check distance to the ring
                        if len(path_to_ring) > 0:
                            if len(path_to_ring[0]) <= 3: ### Same distance as MakeSubgraph (2 bonds, 3 atoms), only gets the shortest
                                path_to_ring[0] = sorted(path_to_ring[0])
                                path_n_ring = list(path_to_ring[0]).copy()
                                n_check = len(path_n_ring)
                                for atom in list(path_n_ring):
                                    for ring in list(has_rings):
                                        if atom in list(ring) and not joined:
                                            path_n_ring += list(ring)
                                            joined = True
                        else:
                            print('@ Was it a ring? @')
                            
                        if joined:
                            path_n_ring = list(set(path_n_ring)) ### make it unique
                            path_n_ring = sorted(path_n_ring)
                            submols_atm.append(path_n_ring)
                            print('\t1.Atom', Symbol(indices[0]), indices[0], 'NOT in a ring, path to ring and ring added:', path_n_ring, Chem.MolFragmentToSmiles(newMol, path_n_ring, isomericSmiles=False))

                        ### If the ring is > 2 bonds away
                        if not joined:
                            try:
                                print('\t1.Atom', Symbol(indices[0]), indices[0], 'NOT in a ring, path to ring > 2 bonds (', str(len(path_to_ring[0]) - 1), '). Making subgraphs by fixed length.')
                            except:
                                print('\t1.Atom', Symbol(indices[0]), indices[0], 'NOT in a ring, path to ring EMPTY. It may not have a ring! Making subgraphs by fixed length.')
                            temp_sets = MakeSubgraph(indices, 3)
                            for key in temp_sets.keys():
                                print('\t1.Subgraph up to', key, 'bonds:', temp_sets[key], Chem.MolFragmentToSmiles(newMol, temp_sets[key], isomericSmiles=False))
                                if temp_sets[key] not in submols_atm:
                                    submols_atm.append(temp_sets[key])
                                    already_subgraph = True
                            
                if not has_rings:
                    print(lig_mol, 'has no rings. Making subgraphs by fixed length.')
                    if lig_mol not in no_rings:
                        no_rings.append(lig_mol)
                    temp_sets = MakeSubgraph(indices, 3)
                    for key in temp_sets.keys():
                        print('\t1.Subgraph up to', key, 'bonds:', temp_sets[key], Chem.MolFragmentToSmiles(newMol, temp_sets[key], isomericSmiles=False))
                        if temp_sets[key] not in submols_atm:
                            submols_atm.append(temp_sets[key])
                            already_subgraph = True
    #1
                ### "terminal"
                if not already_subgraph:
                    terminal = MakeSubgraph(indices, 2) ## or submols_atm??
                    if len(submols_atm) > 1:
                        print('HOW? RECHECK IT!')
                        print('$$submols_atm', submols_atm) ### remove after 1 successful run
                    if submols_atm[0]:
                        n_check = len(submols_atm[0])
                        add_term = submols_atm[0].copy()
                        for key in terminal.keys(): ## Will have only one key
                            for atm in terminal[key]:
                                if atm not in add_term:
                                    if atm not in all_ring_atoms:
                                        add_term.append(atm)
                                    else:
                                        excluded_term.append([Symbol(atm), atm])

                        add_term = list(set(add_term)) ## make unique
                        add_term = sorted(add_term)
                        if len(add_term) > n_check:
                            if add_term not in submols_atm:
                                submols_atm.append(add_term)
                                has_terminal = True
                                print("\t1.Atom's submol + terminal:", add_term, Chem.MolFragmentToSmiles(newMol, add_term, isomericSmiles=False))

            if length > 1:
                short_path, short_atm_in_ring = [], []
                do_extra, t_extra_bond, stop1 = False, False, False
                if lig_mol not in more_than_1:
                    more_than_1.append(lig_mol)

                remark1 = '\t' + str(length) + '.Shortest Path:'
                remark2 = '\t' + str(length) + ".At least one of Shortest Path's atoms in a ring, ring added to Shortest Path:"
                remark3 = remark1[:-1] + ' + terminal, 2 bonds:'
                remark4 = '\t' + str(length) + '.Shortest Path + ring + terminal, 2 bonds:'

                ### Get shortest path between all near atoms
                short_path_problem = []
                for x in range(0, length - 1):
                    for y in range(x+1, length):
                        short_path += Chem.GetShortestPath(newMol, indices[x], indices[y])
                        if not short_path:
                            short_path_problem.append([indices[x], indices[y]])
                            stop1 = True

                ### Empty short paths, should not continue! ¤ Current ligand: R1D (PDB ID: 7MPF, MN01)
                if stop1:
                    print('EMPTY SHORT PATH(S):', short_path_problem)
                    print('Will skip this molecule!')
                    short_path.clear()
                    empty_path.append(name)
                    continue

                ### Check if any of shortest path' atoms is in a ring
                if short_path:
                    short_path = list(set(short_path)) ## make unique
                    short_path = sorted(short_path)
                    n_check1 = len(short_path)
                    print(remark1, short_path, Chem.MolFragmentToSmiles(newMol, short_path, isomericSmiles=False))

                    if short_path not in submols_atm:
                        submols_atm.append(short_path)
                    if has_rings:
                        print(lig_mol, 'has rings:', has_rings)
                        all_ring_atoms += ListAtomsInRings(has_rings) ### Make list with all ring atoms
                        for atom in short_path:
                            in_ring = AtomInRing(has_rings, atom) ### If the near atom is in a ring, gets list with ring atoms
                            if in_ring:
                                short_atm_in_ring += in_ring

                    add_term1 = short_path.copy() ### "terminal"

                    if short_atm_in_ring:
                        short_atm_in_ring += short_path
                        short_atm_in_ring = list(set(short_atm_in_ring)) ## make unique
                        short_atm_in_ring = sorted(short_atm_in_ring)
                        n_check2 = len(short_atm_in_ring)
                        if short_atm_in_ring not in submols_atm:
                            submols_atm.append(short_atm_in_ring)
                            print(remark2, short_atm_in_ring, Chem.MolFragmentToSmiles(newMol, short_atm_in_ring, isomericSmiles=False))

                        add_term2 = short_atm_in_ring.copy() ### "terminal"

                else:
                    print(remark1, 'could NOT be created!')

                ### "terminal"
                if length <= 3: ### To avoid being too large!
                    terminal = MakeSubgraph(indices, 2)
                    for key in terminal.keys(): ## It'll have only one key
                        for atm in terminal[key]:
                            if atm not in short_path:
                                if atm not in all_ring_atoms:
                                    add_term1.append(atm)
                                else:
                                    excluded_term.append([Symbol(atm), atm])
                            if atm not in short_atm_in_ring:
                                if atm not in all_ring_atoms:
                                    add_term2.append(atm)
                                else:
                                    excluded_term.append([Symbol(atm), atm])

                try:
                    add_term1 = list(set(add_term1)) ## make unique
                    add_term1 = sorted(add_term1)
                    if len(add_term1) > n_check1:
                        if add_term1 not in submols_atm:
                            submols_atm.append(add_term1)
                            has_terminal = True
                            print(remark3, add_term1, Chem.MolFragmentToSmiles(newMol, add_term1, isomericSmiles=False))
                except Exception as e:
                    print('>ERROR1:', e) ### TAKE IT OUT AFTER 1 SUCCESSFULL RUN
                    pass

                if short_atm_in_ring:
                    try:
                        add_term2 = list(set(add_term2)) ## make unique
                        add_term2 = sorted(add_term2)
                        if len(add_term2) > n_check2:
                            if add_term2 not in submols_atm:
                                submols_atm.append(add_term2)
                                has_terminal = True
                                print(remark4, add_term2, Chem.MolFragmentToSmiles(newMol, add_term2, isomericSmiles=False))
                    except Exception as e:
                        print('>ERROR2:', e) ### TAKE IT OUT AFTER 1 SUCCESSFULL RUN
                        pass

#####################################################################
            ##################################################################### For all molecules made
            submols_atm.sort(key = len) ## Make sure that the largest submol is last
            near_w_h = []
            if excluded_term:
                print('\tExcluded terminal atom(s) in a ring:', excluded_term)

            smiles_w_hs[name] = Chem.MolToSmiles(newMol, isomericSmiles=False)

            for idx, sublist in enumerate(submols_atm):
                ## make list with explicit Hs for submol
                addhs = []
                try:
                    newMol = Chem.AddHs(newMol, explicitOnly=True, onlyOnAtoms=sublist) # explicitOnly=True ? ### to draw it
                    newMol = Chem.AddHs(newMol, onlyOnAtoms=indices) # To guarantee Hs in near atoms! ## N-H, O-H! ## remove it to make pictures OH, but it changes unique submols made
                    for atm in newMol.GetAtoms():
                        atm = atm.GetIdx()
                        if Symbol(atm) == 'H':
                            for sub_atm in sublist:
                                if Symbol(sub_atm) not in dont_addh: ## Don't add Hs for C or P (3L5, PDB ID: 4RAD)
                                    dist_to_sub_atm = Chem.GetShortestPath(newMol, atm, sub_atm)
                                    if len(dist_to_sub_atm) == 2:
                                        addhs.append(atm)
                                        near_w_h.append(sub_atm)

                    mol2d = Chem.Mol(newMol)
                    AllChem.Compute2DCoords(mol2d) ## FROM https://github.com/rdkit/rdkit/discussions/5316

                except:
                    print('Could not add Hs')
                    continue

                all_submol = sublist + addhs
                all_submol = list(set(all_submol))

                if make_pic: ### see if it actually only makes it once!!
                    try:
                        ## Make picture with submol only if unique set of atoms!!!
                        if not drew_near1:
                            draw_submols(lig_mol, mol2d, new_name, all_submol, idx)
                            drew_near1.setdefault(lig_mol, []).append(all_submol)
                        if lig_mol not in drew_near1.keys():
                            draw_submols(lig_mol, mol2d, new_name, all_submol, idx)
                            drew_near1.setdefault(lig_mol, []).append(all_submol)
                        if all_submol not in drew_near1[lig_mol]:
                            draw_submols(lig_mol, mol2d, new_name, all_submol, idx)
                            drew_near1.setdefault(lig_mol, []).append(all_submol)

                        track = [indices, all_submol]
                        if not drew_near2:
                            draw_near_submols(lig_mol, mol2d, new_name, all_submol, indices, idx)
                            drew_near2.setdefault(lig_mol, []).append(track)
                        if lig_mol not in drew_near2.keys():
                            draw_near_submols(lig_mol, mol2d, new_name, all_submol, indices, idx)
                            drew_near2.setdefault(lig_mol, []).append(track)
                        if track not in drew_near2[lig_mol]:
                            draw_near_submols(lig_mol, mol2d, new_name, all_submol, indices, idx)
                            drew_near2.setdefault(lig_mol, []).append(track)
                        
                    except Exception as e:
                        error2 = f'\tERROR: {name} - Not possible to draw it! {e}'
                        if error2 not in print_later:
                            print_later.append(error2)

                try:
                    frag = Chem.MolFragmentToSmiles(newMol, all_submol, isomericSmiles=False) ### no stereochemistry

                    ## Check if no extra atoms were added (6VB9_OXD_503_MG01)
                    try:
                        submol_heavyatms = Chem.MolFromSmiles(frag).GetNumHeavyAtoms()
                        if submol_heavyatms > template_heavyatms:
                            print('Atention! Submol heavy atoms (' + str(submol_heavyatms) + ') > ligand SMILES (' + str(template_heavyatms) + '). Will pass this ligand!')
                            divergent_n_heavyatms.append(name)
                            frag = None
                    except:
                        pass

                except Exception as e:
                    print('SUBMOL ERROR:', sublist, e)

                if frag:
                    ## fix C:C to c:c
                    if ':' in frag:
                        fixed = FixColon(frag)
                        frag = fixed[0]
                        frag_fixed.append(fixed[1])

                    if '[H]' in frag:
                        fixed2 = FixHs(frag)
                        frag_fixed2.append(fixed2[1])

                    if frag not in submols_smiles:
                        submols_smiles.append(frag)
                        if name not in frag_done:
                            frag_done.append(name)
                        if lig_mol not in frag_done_lig1:
                            frag_done_lig1.append(lig_mol)

            if frag_fixed:
                print('\t' + name + ' FIXED: ', frag_fixed)

            if frag_fixed2:
                print('\t' + name + ' WILL BE FIXED IN ALL_SUBMOLS.TXT: ', frag_fixed2)

            if print_later:
                for error in print_later:
                    print(error)

            near_w_h = list(set(near_w_h))
            for xx in indices:
                if xx not in near_w_h:
                    degree = newMol.GetAtomWithIdx(xx).GetDegree()
                    near_no_h.setdefault(name, []).append((Symbol(xx), xx, degree))

            ### All submolecules made with SMILES
            print('Submols (' + str(len(submols_smiles)) + '): ' + str(submols_smiles))
            if submols_smiles:
                submol_dict.setdefault(length, {}).setdefault(lig_mol, {}).setdefault(name, submols_smiles) ### nested dictionary
                lig_x_near.setdefault(lig_mol, []).append(length)

                ### "Broken" molecules
                BrokenSmiles(name, submols_smiles)

                ### Make dictionary with only the largest SMILES
                last_smi.setdefault(length, {}).setdefault(lig_mol, [])
                for smile01 in submols_smiles:
                    if smile01 == submols_smiles[-1]:
                        if '.' in smile01:
                            save_smi = submols_smiles[-2] ## Get second to last
                        else:
                            save_smi = smile01

                if save_smi not in last_smi[length][lig_mol]:
                    last_smi[length][lig_mol].append(save_smi) ### Will be printed

                to_group.setdefault(length, {}).setdefault(save_smi, []).append(lig_mol)

            print(f"\t     >>> Will be selected: '{save_smi}'")

            ### Get largest submolecule
            submol_sorted = []
            for submol in submols_smiles:
                n_char = 0
                for char in submol:
                    if char.isalpha():
                        n_char += 1
                submol_sorted.append((n_char, submol))
              
            submol_sorted = sorted(submol_sorted)

            if len(submol_sorted) > 1:
                try:
                    if submol_sorted[-1][1] != submols_smiles[-1]:
                        print('%%%NOT THE SAME%%') ## REMOVE AFTER 1 SUCCESSFUL RUN
                except Exception as e:
                    print('\tERROR! Not possible to get path at all!', e)
                    if name not in no_surroundings:
                        no_surroundings.append(name)

            print('\t\t\t\t\t\t\t\t\t\t#', len(frag_done_lig1), 'ligands with substructures made so far')

            ### Make new overview of number of near atoms
            only_key = f'only_{len(indices)}'
            only.setdefault(only_key, [])
            if lig_mol not in only[only_key]:
                only[only_key].append(lig_mol)


print('\n                    ┌───────────────┐')
print('                    │   ## END ##   │')
print('                    └───────────────┘\n')


############ ANALYSIS ############
has_near_lig, frag_done_lig = [], []

for name in has_near:
    lig_mol = name.rsplit('_')[1]
    if lig_mol not in has_near_lig:
        has_near_lig.append(lig_mol)

for name in frag_done:
    lig_mol = name.rsplit('_')[1]
    if lig_mol not in frag_done_lig:
        frag_done_lig.append(lig_mol)

### No newMol at all!
frag_missing = []
for k in sorted(no_newMol.keys()):
    if k not in list(frag_done_lig):
        if k not in list(frag_missing):
            frag_missing.append(k) 
    if len(no_newMol[k]) > 1:
        no_newMol[k] = sorted(list(set(no_newMol[k])))

for key in sorted(lig_x_near.keys()):
    lig_x_near[key] = list(set(lig_x_near[key]))
    lig_x_near[key] = sorted(lig_x_near[key])

### WHOLE DICTIONARY
whole_file = f'{metal}_all_submols.txt'

# Make file with substructures SMILES
os.makedirs(f'{mf}/submols_{timestamp}', exist_ok=True)
with open(f'{mf}/submols_{timestamp}/{whole_file}', 'w+') as exit_file:
    exit_file.write(f'{dt1_string}\n')
    exit_file.write(f'{mf}\n')
    exit_file.write(f'All substructures SMILES made for this set of ligands\n\n')
    for near_atoms, n_values in sorted(submol_dict.items()):
        if near_atoms == 1:
            exit_file.write(f'§§§ {near_atoms} near atom\n')
        if near_atoms > 1:
            exit_file.write(f'\n§§§ {near_atoms} near atoms\n')
        for lig, sdf_smiles in sorted(n_values.items()):
            if len(lig_x_near[lig]) == 1:
                if lig_x_near[lig] == 1:
                    appear_in = f'(Has only {str(lig_x_near[lig])[1:-1]} near atom)'
                if lig_x_near[lig] != 1:
                    appear_in = f'(Has only {str(lig_x_near[lig])[1:-1]} near atoms)'
            if len(lig_x_near[lig]) > 1:
                appear_in = f'(Can have {str(lig_x_near[lig][:-1])[1:-1]} or {str(lig_x_near[lig][-1])} near atoms)'
            exit_file.write(f'\t# {lig} {appear_in}\n')
            for sdf, smiles in sorted(sdf_smiles.items()):
                smiles_h = smiles_w_hs[sdf]
                space1 = 50 - len(str(smiles))
                if space1 >= 4:
                    space = ' '*space1
                if space1 < 4:
                    space = ' '*4
                exit_file.write(f'\t    {sdf}\t::\t{smiles}\t{space}:MOLECULE:\t{smiles_h}\n')
            exit_file.write('\n')

print(f'## Whole dictionary with SMILES saved as {mf}/submols_{timestamp}/{whole_file}')

with open(f'{mf}/submols_{timestamp}/submols_id.txt', 'w+') as exit_file:
    exit_file.write(f"## Dictionary with the pictures' submol id == SMILES\n")
    for idx, submol in enumerate(done_submols):
        submol_n = str('{:04d}'.format(idx))
        exit_file.write(f"\t\tsubmol-{submol_n}\t==\t'{submol}'\n")
	

print(f"## Dictionary with the pictures' submol id == SMILES saved as {mf}/submols_{timestamp}/submols_id.txt")

print('\n=======================================================================================================================================\n')

print(f'## Dictionary with unique submols by order and the ones that will be grouped per ligand')

unique_smiles = {}
only_1_submol = []

for near_atoms, n_values in sorted(submol_dict.items()):
    if near_atoms == 1:
        print('\n§§§', near_atoms, 'near atom')
    if near_atoms > 1:
        print('\n§§§', near_atoms, 'near atoms')
    for lig, sdf_smiles in sorted(n_values.items()):
        print_once = False
        length = str(near_atoms)
        unique_smiles.setdefault(length, {}).setdefault(lig, {})
        if len(lig_x_near[lig]) == 1:
            if lig_x_near[lig] == 1:
                appear_in = '(Has only ' + str(lig_x_near[lig])[1:-1] + ' near atom)'
            if lig_x_near[lig] != 1:
                appear_in = '(Has only ' + str(lig_x_near[lig])[1:-1] + ' near atoms)'
        if len(lig_x_near[lig]) > 1:
            appear_in = '(Can have ' + str(lig_x_near[lig][:-1])[1:-1] + ' or ' + str(lig_x_near[lig][-1]) + ' near atoms)'
        print('\t#', lig, appear_in)
        for sdf, smiles in sorted(sdf_smiles.items()):
            for idx, smile in enumerate(smiles):
                unique_smiles[length][lig].setdefault(idx, [])
                if smile not in unique_smiles[length][lig][idx]:
                    unique_smiles[length][lig][idx].append(smile)

                ### Make dictionary with only the largest SMILES
                if smile == smiles[-1]:
                    if idx == 0:
                        only_1_submol.append(lig)

                if sdf in near_no_h.keys():
                    if near_no_h[sdf][0][0] == 'O' and near_no_h[sdf][0][2] == 1:
                        if '=' not in smile and not print_once:
                            print('$$ Near atom has no Hs:', near_no_h[sdf])
                            print_once = True
                    if near_no_h[sdf][0][0] != 'O' and not print_once:
                        print('$$ Near atom has no Hs:', near_no_h[sdf])
                        print_once = True

        print('\t    ', str(unique_smiles[length][lig])[1:-1])
        print('')

only_1_submol = list(set(only_1_submol)) ### make it unique, will print in STATS
only_1_submol = sorted(only_1_submol)

print('\n=======================================================================================================================================\n')

print('------------------------------ Grouped substructures ------------------------------')

all_submols_grouped = Group(to_group)

print('\n=======================================================================================================================================\n')

print('## Ligands per number of near atoms after removal')
for key in sorted(list(only.keys())):
    number = key.split('_')[1]
    if number == '1':
        print('\n' + str(len(only[key])) + ' ligand(s) with ' + number + ' near atom:\t' + str(sorted(only[key])))
    else:
        print('\n' + str(len(only[key])) + ' ligand(s) with ' + number + ' near atoms:\t' + str(sorted(only[key])))

print('\n=======================================================================================================================================\n')

print('\n================================== STATS ==================================\n')

print('### Number of SDFs:', len(mols))
print('\n\t# SDFs with near atoms:', len(has_near), Percentage(has_near, mols))
print('\n\t# SDFs with generated substructures:', len(frag_done), Percentage(frag_done, mols))
print('\n\t# SDFs with empty paths and therefore no substructures:', len(empty_path), Percentage(empty_path, mols))
if len(empty_path) >= 1:
    print(empty_path)
print('\n\t# SDFs without surroundings:', len(no_surroundings), Percentage(no_surroundings, mols))
if len(no_surroundings) >= 1:
    print(no_surroundings)
divergent_n_heavyatms = list(set(divergent_n_heavyatms))
print("\n\t# SDFs with number of heavy atoms != ligand's SMILES:", len(divergent_n_heavyatms), Percentage(divergent_n_heavyatms, mols))
if len(divergent_n_heavyatms) >= 1:
    print(divergent_n_heavyatms)
for key, values in broken.items():
     if '.' in key:
        print(f"\n\t# Substructures per SDF with '{key}':", len(values), Percentage(values, mols))
        if len(values) >= 1:
            print(values)
     if ':' in key:
        print(f"\n\t# Substructures per SDF with '{key}':", len(values), Percentage(values, mols))

print('\n### Number of molecules:', len(lig_list))
print('\n\t# Molecules with near atoms:', len(has_near_lig), Percentage(has_near_lig, lig_list))
print('\n\t# Molecules with generated substructures:', len(frag_done_lig), Percentage(frag_done_lig, lig_list))
print('\t\t# Molecules with only 1 unique substructure:', len(only_1_submol), Percentage(only_1_submol, frag_done_lig), only_1_submol)
print('\n\t# Molecules without SMILES:', len(no_smi), Percentage(no_smi, lig_list))
print('\n\t# Molecules without newMol:', len(no_newMol.keys()), Percentage(list(no_newMol.keys()), lig_list))
print('\t\t# Molecules without NewMol in all instances:', len(frag_missing), Percentage(frag_missing, lig_list), sorted(frag_missing))
print('\n\t# Molecules with only 1 near atom:', len(only_1), Percentage(only_1, lig_list))
print('\t\t# Aliphatic molecules:', len(no_rings), Percentage(no_rings, only_1))
print('\n\t# Molecules with >= 2 near atoms:', len(more_than_1), Percentage(more_than_1, lig_list))

print('\n### ALL submols grouped:', len(all_submols_grouped))

print('\n===========================================================================\n')

end = datetime.now()
dt2_string = end.strftime('%d/%m/%Y %H:%M:%S')
print('\n-------------------------------\nTimestamp:', dt2_string)