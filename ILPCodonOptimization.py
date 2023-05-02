import itertools

import gurobipy as gb
import pandas as pd
import math
import argparse
import RNA

#
# Inverse table for standard genetic code
# Amino acid -> DNA codons list
# Stop codons are denoted by '#'
#

InvTableDNA = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'K': ['AAA', 'AAG'],
    'N': ['AAT', 'AAC'],
    'M': ['ATG'],
    'D': ['GAT', 'GAC'],
    'F': ['TTT', 'TTC'],
    'C': ['TGC', 'TGT'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'E': ['GAA', 'GAG'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'W': ['TGG'],
    'H': ['CAT', 'CAC'],
    'Y': ['TAT', 'TAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    '#': ['TAA', 'TGA', 'TAG']
}

#
# Reverse InvTableDNA
#

TableDNA = {i: k for k, v in InvTableDNA.items() for i in v}

#
# Codon normalized frequencies
#

#
# Organism:
# Escherichia coli O157:H7 EDL933 [gbbct]: 5347 CDS's (1611503 codons)
# Source: http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=155864
#

FrequencyNorm = {'TTT': 1.00, 'TTC': 0.71, 'TTA': 0.27, 'TTG': 0.25, 'CTT': 0.22,
                 'CTC': 0.20, 'CTA': 0.08, 'CTG': 1.00, 'ATT': 1.00, 'ATC': 0.80,
                 'ATA': 0.18, 'ATG': 1.00, 'GTT': 0.69, 'GTC': 0.56, 'GTA': 0.42,
                 'GTG': 1.00, 'TCT': 0.97, 'TCC': 1.00, 'TCA': 0.91, 'TCG': 0.99,
                 'CCT': 0.32, 'CCC': 0.25, 'CCA': 0.37, 'CCG': 1.00, 'ACT': 0.40,
                 'ACC': 1.00, 'ACA': 0.35, 'ACG': 0.66, 'GCT': 0.48, 'GCC': 0.78,
                 'GCA': 0.64, 'GCG': 1.00, 'TAT': 1.00, 'TAC': 0.74, 'TAA': 1.00,
                 'TAG': 0.14, 'CAT': 1.00, 'CAC': 0.73, 'CAA': 0.50, 'CAG': 1.00,
                 'AAT': 0.88, 'AAC': 1.00, 'AAA': 1.00, 'AAG': 0.32, 'GAT': 1.00,
                 'GAC': 0.58, 'GAA': 1.00, 'GAG': 0.48, 'TGT': 0.82, 'TGC': 1.00,
                 'TGA': 0.55, 'TGG': 1.00, 'CGT': 0.97, 'CGC': 1.00, 'CGA': 0.18,
                 'CGG': 0.30, 'AGT': 0.59, 'AGC': 1.00, 'AGA': 0.14, 'AGG': 0.09,
                 'GGT': 0.86, 'GGC': 1.00, 'GGA': 0.32, 'GGG': 0.42}
#
# Organism:
# Pichia pastoris [gbpln]: 137 CDS's (81301 codons)
# Source: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4922&aa=1&style=N
#
FrequencyNormPP = {'TTT': 1.00, 'TTC': 0.85, 'TTA': 0.48, 'TTG': 1.00, 'CTT': 0.48,
                 'CTC': 0.24, 'CTA': 0.33, 'CTG': 0.48, 'ATT': 1.00, 'ATC': 0.62,
                 'ATA': 0.36, 'ATG': 1.00, 'GTT': 1.00, 'GTC': 0.54, 'GTA': 0.35,
                 'GTG': 0.45, 'TCT': 1.00, 'TCC': 0.67, 'TCA': 0.62, 'TCG': 0.31,
                 'CCT': 0.83, 'CCC': 0.35, 'CCA': 1.00, 'CCG': 0.21, 'ACT': 1.00,
                 'ACC': 0.65, 'ACA': 0.60, 'ACG': 0.27, 'GCT': 1.00, 'GCC': 0.57,
                 'GCA': 0.51, 'GCG': 0.13, 'TAT': 0.88, 'TAC': 1.00, 'TAA': 1.00,
                 'TAG': 0.56, 'CAT': 1.00, 'CAC': 0.75, 'CAA': 1.00, 'CAG': 0.63,
                 'AAT': 1.00, 'AAC': 0.92, 'AAA': 0.88, 'AAG': 1.00, 'GAT': 1.00,
                 'GAC': 0.72, 'GAA': 1.00, 'GAG': 0.78, 'TGT': 1.00, 'TGC': 0.56,
                 'TGA': 0.39, 'TGG': 1.00, 'CGT': 0.35, 'CGC': 0.10, 'CGA': 0.20,
                 'CGG': 0.10, 'AGT': 0.51, 'AGC': 0.31, 'AGA': 1.00, 'AGG': 0.33,
                 'GGT': 1.00, 'GGC': 0.31, 'GGA': 0.75, 'GGG': 0.22}
#
# Organism:
# Arabidopsis thaliana [gbpln]: 80395 CDS's (31098475 codons)
# Source: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=3702
#
FrequencyNormAT = {'TTT': 1.00, 'TTC': 0.96, 'TTA': 0.53, 'TTG': 0.48, 'CTT': 1.00,
                 'CTC': 0.65, 'CTA': 0.42, 'CTG': 0.42, 'ATT': 1.00, 'ATC': 0.85,
                 'ATA': 0.58, 'ATG': 1.00, 'GTT': 1.00, 'GTC': 0.47, 'GTA': 0.37,
                 'GTG': 0.65, 'TCT': 1.00, 'TCC': 0.46, 'TCA': 0.71, 'TCG': 0.35,
                 'CCT': 1.00, 'CCC': 0.28, 'CCA': 0.86, 'CCG': 0.47, 'ACT': 1.00,
                 'ACC': 0.58, 'ACA': 0.91, 'ACG': 0.44, 'GCT': 1.00, 'GCC': 0.37,
                 'GCA': 0.62, 'GCG': 0.32, 'TAT': 1.00, 'TAC': 0.92, 'TAA': 0.81,
                 'TAG': 0.45, 'CAT': 1.00, 'CAC': 0.63, 'CAA': 1.00, 'CAG': 0.78,
                 'AAT': 1.00, 'AAC': 0.92, 'AAA': 0.96, 'AAG': 1.00, 'GAT': 1.00,
                 'GAC': 0.47, 'GAA': 1.00, 'GAG': 0.92, 'TGT': 1.00, 'TGC': 0.66,
                 'TGA': 1.00, 'TGG': 1.00, 'CGT': 0.48, 'CGC': 0.20, 'CGA': 0.34,
                 'CGG': 0.25, 'AGT': 0.57, 'AGC': 0.46, 'AGA': 1.00, 'AGG': 0.57,
                 'GGT': 0.91, 'GGC': 0.37, 'GGA': 1.00, 'GGG': 0.43}

tRNAAT  =   {'GCT': 21.67, 'GCC': 11.0, 'GCA': 17.67, 'GCG': 15.67, 'TGT': 4.67,
             'TGC': 14.0,'GAT': 8.33,'GAC': 25.0,'GAA': 14.33,'GAG': 16.33,
             'TTT': 5.0,'TTC': 15.0,'GGT': 12.67,'GGC': 26.67,'GGA': 20.67,
             'GGG': 16.0,'CAT': 3.0,'CAC': 9.0,'ATT': 18.67,'ATC': 7.33,
             'ATA': 10.67,'AAA': 18.67,'AAG': 21.33,'TTA': 9.33,'TTG': 12.0,
             'CTT': 15.0,'CTC': 7.67,'CTA': 13.67,'CTG': 9.67,'ATG': 9.0,
             'AAT': 5.0,'AAC': 15.0,'CCT': 31.0,'CCC': 21.0,'CCA': 49.67,
             'CCG': 24.33,'CAA': 10.0,'CAG': 11.33,'CGT': 12.33,'CGC': 6.33,
             'CGA': 10.33,'CGG': 9.0,'AGA': 11.33,'AGG': 10.0,'TCT': 38.67,
             'TCC': 15.33,'TCA': 20.0,'TCG': 18.0,'AGT': 3.67,'AGC': 11.0,
             'ACT': 14.33,'ACC': 7.67,'ACA': 13.0,'ACG': 11.0,'GTT': 19.0,
             'GTC': 9.67,'GTA': 14.33,'GTG': 15.0,'TGG': 13.0,'TAT': 20.33,'TAC': 61.0}
#
# Codon Adaptation Index costs to optimize CAI
#

Costs = {k: math.log(FrequencyNorm[k]) for k in FrequencyNorm.keys()}


#
# Functions
#


def read_fasta(fasta_file: str) -> pd.DataFrame:
    """
    Read fasta file using pandas
    :param: filename
    :return: pandas DataFrame with target proteins
    """
    df = pd.read_csv(fasta_file, header=None, comment='>', engine='python',
                     names=['Protein'])

    df['Protein'] = df['Protein'].apply(lambda x: x[:len(x) - len(x) % 3])

    df['Length'] = df['Protein'].str.len()

    print(fasta_file, ('has %d protein(s)' % df.shape[0]))

    return df


def read_motifs(forbidden_file: str,
                desired_file: str) -> (list, list):
    """
    Read forbidden and desired motifs
    :param: forbidden and desired file names
    :return: lists of forbidden and desired motifs
    """

    forbidden = []
    desired = []

    if forbidden_file is not None:
        with open(forbidden_file) as file:
            forbidden = file.readlines()

        forbidden = [x.strip() for x in forbidden if len(x.strip()) > 3]

        print(f'Number of forbidden motifs: {len(forbidden)}')

    if desired_file is not None:
        with open(desired_file) as file:
            desired = file.readlines()

        desired = [x.strip() for x in desired if len(x.strip()) > 3]

        print(f'Number of desired motifs: {len(desired)}')

    return forbidden, desired


def cds_from_solution_y(x: gb.tupledict) -> str:
    """
    Build the protein from a solver's solution (x variables only)
    :param: solution x
    :return: protein
    """
    cds = []
    for i in x:
        if x[i].X > 0.5:
            cds.append(i[1])

    return ''.join(cds)


def find_amino_from_codon(codon: str,
                          table: dict) -> dict:
    """
    Find all possible aminos that can be codified by a (partly empty) codon
    :param: a 3-chars string representing a codon. May contain '*' as don't care
    either at the beginning or at the end of the sequence
    :return: a dictionary with all possible amino acids that the sequence codifies
    """

    aminos = {}

    if codon.isalpha():
        if table[codon] in aminos:
            aminos[table[codon]].append(codon)
        else:
            aminos[table[codon]] = [codon]

    if codon[0] == '*':
        partial_codon = codon.lstrip('*')
        offset = len(codon) - len(partial_codon)
        for i in table:
            if partial_codon == i[offset:]:
                if table[i] in aminos:
                    aminos[table[i]].append(i)
                else:
                    aminos[table[i]] = [i]

    if codon[2] == '*':
        partial_codon = codon.rstrip('*')
        offset = len(codon) - len(partial_codon)

        for i in table:
            if partial_codon == i[:3 - offset]:
                if table[i] in aminos:
                    aminos[table[i]].append(i)
                else:
                    aminos[table[i]] = [i]

    return aminos


def find_amino_sequence_in_protein(sequence: str,
                                   aminos: str) -> list:
    """
    Check if a non-empty sequence of amino acids is
    contained in the aminos string
    :param: two sequences (strings) of amino acids
    :return: an (eventually empty) list of positions
    """
    feasiblepos = []

    pos = aminos.find(sequence)

    while pos != -1:
        restart = pos
        start = pos
        end = (pos + len(sequence) - 1)
        feasiblepos.append((start, end))
        pos = aminos.find(sequence, restart + 1)

    return feasiblepos


def bases_per_nucleotide(amino: str) -> dict:
    """
        Provides the bases compatible with each nucleotide in the given
        amino acid sequence
        :param: a sequence (string) of amino acids
        :return: a dictionary containing, for each base position, the set of admissible nucleotide
        """
    protein_length = len(amino) * 3
    bases = {i: set() for i in range(protein_length)}
    for idx, a in enumerate(amino):
        for codon in InvTableDNA[a]:
            for i in range(3):
                bases[idx * 3 + i].add(codon[i])
    return bases


def protein_from_solution(x: gb.tupledict,
                          protein: list) -> str:
    """
    Build the protein from a solver's solution (x variables only)
    :param: solution x
    :return: protein
    """
    for i in x:
        if x[i] > 0.5:
            position = i[0]
            codon = i[2]
            protein[position] = codon

    return ''.join(protein)

def cds_from_solution(x: gb.tupledict) -> str:
    """
    Build the protein from a solver's solution (x variables only)
    :param: solution x
    :return: protein
    """
    cds = ''
    for i in x:
        if x[i] > 0.5:
            base = i[1]
            cds += base

    return cds


def find_motif_position(motif: str,
                        aminos: str) -> (dict, dict):
    """
    Given a sequence of amino acids find positions in which
    a motif may be located
    :param:  motif, amino acids sequense
    :return: a dictionaries with feasible positions and
    corresponding codon encoding
    Positions starts from 0 and refer to basis position
    Codon encoding at first and last level may contain
    more than one alternative
    """

    feasible_positions = {}
    feasible_positions_basis = {}

    #
    # Test three possible positions
    #

    for offset in range(3):

        sequence = ''.join('*' * offset) + motif
        start_offset = offset

        if len(sequence) % 3 != 0:
            end_offset = (3 - (len(sequence) % 3))
            sequence += ('*' * (3 - (len(sequence) % 3)))
        else:
            end_offset = 0

        first_level = find_amino_from_codon(sequence[0:3], TableDNA)

        if len(sequence) > 3:
            last_level = find_amino_from_codon(sequence[-3:], TableDNA)
        else:
            raise ValueError('Motif must have length > 3')

        middle_levels = []

        for idx in range(3, len(sequence) - 3, 3):
            codon = sequence[idx: idx + 3]
            middle_levels.append(find_amino_from_codon(codon, TableDNA))

        #
        # Build the sequences starting from possible
        # encodings
        #

        for head in first_level:
            sequence = [{head: first_level[head]}]

            sequence += middle_levels

            amino_sequence = ''.join([[*i.keys()][0] for i in middle_levels])
            amino_sequence = head + amino_sequence

            for tail in last_level:
                sequence += [{tail: last_level[tail]}]
                amino_sequence += tail

                feasible_pos_list = find_amino_sequence_in_protein(amino_sequence, aminos)

                for pos in feasible_pos_list:
                    feasible_positions[pos] = sequence.copy()
                    start_basis_pos = pos[0] * 3 + start_offset
                    end_basis_pos = (pos[1] + 1) * 3 - end_offset - 1
                    feasible_positions_basis[start_basis_pos, end_basis_pos] = sequence.copy()

                del sequence[-1]
                amino_sequence = amino_sequence[:-1]

    return feasible_positions_basis


def build_yf_variables_index(forbidden: list,
                             aminos: str) -> (list, dict, dict):
    """
    Build the index set of yf variables (position variables referred to forbidden motifs)
    :param: list of forbidden motifs, amino acids sequence
    :return list of yf indexes, feasible positions for the basis, index of forbidden motifs
    """

    index_list = []
    feasible_positions_forbidden_basis = {}
    index_forbidden = {}

    count = 1

    for motif in forbidden:
        index_forbidden[motif] = count
        count += 1

        aux = find_motif_position(motif, aminos)
        if aux:

            feasible_positions_forbidden_basis[motif] = aux

            for pos in feasible_positions_forbidden_basis[motif]:
                index_list.append((index_forbidden[motif], str(pos).replace(" ", "")))

    return index_list, feasible_positions_forbidden_basis, index_forbidden


def f_gc(sequence):
    w = 30
    total = 0
    gc_ideal = 60
    gc_max = 80
    gc_min = 40
    gc_pan = 3000
    gc_diff = 200
    for i in range(len(sequence) - w):
        gc_local = sequence[i:i+w].count('G') + sequence[i:i+w].count('C')
        total += (abs(gc_ideal - gc_local) * gc_diff) if gc_min <= gc_local <= gc_max else (abs(gc_ideal - gc_local) * gc_diff + gc_pan)
    return total


def v_gc(sequence):
    w = 30
    total = 0
    gc_ideal = 51
    for i in range(len(sequence) - w):
        gc_local = (sequence[i:i + w].count('G') + sequence[i:i + w].count('C')) / w * 100
        total += math.pow(gc_local - gc_ideal, 2)
    return total / (len(sequence) - w)


def count_forbidden(protein: str,
                    forbidden: list) -> int:
    """
    Count the number of forbidden motifs in protein
    """
    total_forb = 0
    for motif in forbidden:
        found = protein.find(motif)
        count = 0
        while found != -1:
            count += 1
            start = found + 1
            found = protein.find(motif, start)
        total_forb += count

    return total_forb


def gen_protein_hierarchical_objectives(target_protein, lpercentage, upercentage, forbidden):

    # Setup the GUROBI Model

    model = gb.Model()
    model.Params.OutputFlag = 0

    index_list = []
    costs = {}

    protein_length = len(target_protein)

    aminos = ''

    # Build index and costs of x variables

    for h in range(protein_length // 3):
        amino = TableDNA[target_protein[h * 3: h * 3 + 3]]
        aminos += amino
        for i in enumerate(InvTableDNA[amino]):
            index_list.append((h, TableDNA[target_protein[h * 3: h * 3 + 3]], i[1]))
            costs[index_list[-1]] = Costs[i[1]] # CAI costs

    bases = bases_per_nucleotide(aminos)
    bases_list = []
    for key, elem in bases.items():
        for e in elem:
            bases_list.append((key, e))
            # costs[bases_list[-1]] = 1 if e == 'G' or e == 'C' else 0

    x = model.addVars(index_list, vtype='B', name='x')
    y = model.addVars(bases_list, vtype='B', name='y')
    model.setObjectiveN(x.prod(costs), 1, 0, name='CAI')
    model.ModelSense = -1

    index_list, feasible_positions_forbidden, index_forbidden = \
        build_yf_variables_index(forbidden, aminos)

    if index_list:
        yf = model.addVars(index_list, vtype='B', name='yf')

        for motif in feasible_positions_forbidden:

            for pos in feasible_positions_forbidden[motif]:
                lhs = gb.LinExpr()
                rhs = 0.0
                aux = [codon[h] for codon in feasible_positions_forbidden[motif][pos] for h in codon]
                for idx, codlist in enumerate(aux):
                    lhs += gb.quicksum([x[pos[0] // 3 + idx, TableDNA[cod], cod] for cod in codlist])
                    rhs += 1.0
                model.addConstr(lhs - yf[index_forbidden[motif], str(pos).replace(" ", "")] <= rhs - 1,
                                name='ForbMotif[' + str(index_forbidden[motif]) + ']' + str(pos).replace(" ", ""))
        model.setObjectiveN(yf.prod({i: -1.0 for i in index_list}), 0, 1, name='Forbidden')

    # Assignment constraints

    model.addConstrs((x.sum(i, '*', '*') == 1 for i in range(protein_length // 3)), name='XPos')
    model.addConstrs((y.sum(i, '*') == 1 for i in range(protein_length)), name='YPos')
    for x_var in x:
        model.addConstr(y[x_var[0]*3, x_var[2][0]] + y[x_var[0]*3 + 1, x_var[2][1]] + y[x_var[0]*3 + 2, x_var[2][2]]
                        >= 3*x[x_var], name='Cover')
    cg = {'C', 'G'}
    cg_list = []
    model.update()
    for key, elem in bases.items():
        for base in cg.intersection(elem):
            cg_list.append(y[key, base])

    #model.addConstr(100 * gb.quicksum(cg_list) <= protein_length * percentage, name='GC_Content')
    short_repeat_length = 10
    gc_nucleotide_window = 30
    ygc = model.addVars(range(protein_length), vtype='I', name='ygc')

    for i in range(protein_length):
        local_gc = gb.LinExpr()
        if i + gc_nucleotide_window >= protein_length:
            break
        for j in range(gc_nucleotide_window):
            for nucleotide in cg:
                if nucleotide not in bases[i+j]:
                    continue
                local_gc += y[i+j, nucleotide]
        model.addConstr(100 * local_gc <= gc_nucleotide_window * upercentage + ygc[i], name='GC_Content_Local_UB_'+str(i))
        model.addConstr(ygc[i] + 100 * local_gc >= gc_nucleotide_window * (lpercentage), name='GC_Content_Local_LB_'+str(i))
    model.setObjectiveN(ygc.prod({i: -1.0 for i in range(protein_length)}), 3, 3, name='GCOverLimit')

    ys = model.addVars(range(protein_length), vtype='B', name='ys')

    for nucleotide in ['A', 'C', 'G', 'T']:
        for i in range(protein_length):
            if i + short_repeat_length >= protein_length:
                break
            flag = False
            short_repeat = gb.LinExpr()
            for j in range(short_repeat_length):
                if nucleotide not in bases[i+j]:
                    flag = True
                    break
                short_repeat += y[i+j, nucleotide]
            if flag:
                continue
            model.addConstr(short_repeat - ys[i] <= short_repeat_length - 1, name="Short_repeat_"+str(nucleotide)+"_"+str(i))
    model.setObjectiveN(ys.prod({i: -1.0 for i in range(protein_length)}), 2, 2, name='SameConsecutiveBases')

    model._vars = x
    model._tentative_protein = list('***' for i in range(protein_length // 3))

    model._current_protein = ''
    # model.Params.PoolSolutions = 1000
    # model.Params.PoolSearchMode = 2
    model.optimize()
    model.write('model.lp')
    if model.Status != gb.GRB.Status.OPTIMAL:
        return '', model.Runtime, model.NodeCount, model.Status, 0

    x_sol = []

    for s in range(model.SolCount):
        model.params.SolutionNumber = s
        x_sol.append(model.getAttr('Xn', x))

    #  Get the best solution for CAI evaluation

    model.params.SolutionNumber = 0
    model.params.ObjNumber = 1
    cai_exp = model.ObjNVal

    final_protein = protein_from_solution(x_sol[0], model._tentative_protein)
    cai = math.pow(math.e, cai_exp / (len(final_protein) // 3))
    #  Get the best solution for CAI evaluation
    #cai_exp = model.ObjVal
    prev_char = ''
    count = 0
    final_protein = cds_from_solution_y(y)

    # for i in range(protein_length):
    #     if prev_char != final_protein[i]:
    #         count = 1
    #         prev_char = final_protein[i]
    #     else:
    #         count += 1
    #     if count >= short_repeat_length:
    #         assert False
    # for forb in forbidden:
    #     assert(forb not in final_protein)
    assert(len(forbidden) > 0)
    cai = math.pow(math.e, cai_exp / (len(final_protein) // 3))

    return final_protein, model.Runtime, model.NodeCount, model.Status, cai


def count_gc(seq):
    c = 0
    for s in seq:
        if s == 'C' or s == 'G':
            c += 1
    return c


def count_repeat(final_protein, short_repeat_length=10):
    prev_char = ''
    count = 0
    total_counter = 0
    flag = True
    for i in range(len(final_protein)):
        if prev_char != final_protein[i]:
            flag = True
            count = 1
            prev_char = final_protein[i]
        else:
            count += 1
        if count >= short_repeat_length and flag:
            flag = False
            total_counter += 1
    return total_counter


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('target_file', help='Target file (FASTA format)')
    parser.add_argument('--forbidden', '-f', help='File of forbidden motifs')
    parser.add_argument('--desired', '-d', help='File of desired motifs')
    parser.add_argument('--lowerpercentage', '-lp', help='Lower bound percentage of CG content')
    parser.add_argument('--upperpercentage', '-up', help='Upper bound percentage of CG content')
    parser.add_argument('--output', '-o', help='Output file', default='optimized')

    args = parser.parse_args()

    dataset = read_fasta(args.target_file)
    lpercentage = int(args.lowerpercentage)
    upercentage = int(args.upperpercentage)
    forbidden, desired = read_motifs(args.forbidden, args.desired) # desired not used atm

    count = 0
    print('Num. prot.  | Length  | Forb. bef. | Des. bef. | CAI bef. | Forb. af.| Des. af.| CAI af. | Time')
    for index, row in dataset.iterrows():
        count += 1

        protein = row['Protein']

        #
        # Calculate starting CAI and number of forbidden and desired motifs
        #
        cai_before = math.pow(math.e,
                              sum(Costs[protein[i * 3:i * 3 + 3]] for i in range(len(protein) // 3)) / (
                                      len(protein) // 3))

        final_protein, time, nodes, status, cai = \
            gen_protein_hierarchical_objectives(protein, lpercentage, upercentage, forbidden)
        if final_protein == '':
            print('PROBLEM INFEASIBLE')
            dataset.loc[index, 'GC-LowerLimit'] = lpercentage
            dataset.loc[index, 'GC-UpperLimit'] = upercentage
            dataset.loc[index, 'Time'] = time
            dataset.loc[index, 'FinalProtein'] = 'INFEASIBLE'
            dataset.loc[index, 'InitialCAI'] = cai_before
            dataset.loc[index, 'FinalCAI'] = cai
            dataset.loc[index, 'Percentage GC content'] = 0
            continue

        print('final prot', final_protein)
        dataset.loc[index, 'GC-LowerLimit'] = lpercentage
        dataset.loc[index, 'GC-UpperLimit'] = upercentage
        dataset.loc[index, 'Time'] = time
        dataset.loc[index, 'FinalProtein'] = final_protein
        dataset.loc[index, 'InitialCAI'] = cai_before
        dataset.loc[index, 'FinalCAI'] = cai
        dataset.loc[index, 'GCVariance'] = v_gc(final_protein)
        dataset.loc[index, 'UnwantedRepeatWF'] = count_repeat(final_protein)
        dataset.loc[index, 'FobriddenMotifWF'] = count_forbidden(final_protein, forbidden)
        perc_gc = 100*count_gc(final_protein)/len(final_protein)
        dataset.loc[index, 'Percentage GC content'] = perc_gc

        print('%11d' % count, '|%8d' % len(protein), ('|   %1.4f' % cai_before),
              ('| %1.4f' % dataset.loc[index, 'FinalCAI']), ('| %4.4f' % time), ('| %4.4f' % perc_gc))
        print('vgc', v_gc(final_protein))
        print('ur', count_repeat(final_protein))
        print('um', count_forbidden(final_protein, forbidden))
    dataset.to_csv(args.output + '_' + str(lpercentage) + '_' + str(upercentage) + '.csv', index=False)


if __name__ == '__main__':
    main()
