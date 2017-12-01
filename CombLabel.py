#!/usr/bin/python3

import tkinter as tk
import csv
import time
import copy
import re
import argparse
import socket # socket.gethostname()
from cl_errors import errors as err
#import residues as rs

# Globals

CITATION = """
#################################################
##
##   When using this software to calculate CSL schemes, please, cite the following article:
## 
##   'CombLabel: rational design of optimized combinatorial labeling schemes.
##   Application to backbone assignment of helical membrane proteins with
##   limited stability'
##   M.Yu. Myshkin, M.A. Dubinnyi, D.S. Kulbatskii, E.N. Lyukmanova, Z.O. Shenkarev 
##   Journal of Biomolecular NMR, 2018 (currnetly unpublished)
##
#################################################
#
"""
CONFIG_FILE = "comblabel.task"

parser = argparse.ArgumentParser()
parser.add_argument("config", default=CONFIG_FILE, help='Specify the config file')
parser.add_argument('--check', '-c', default='-', help='Check solution in the specified file')
parser.add_argument('--number', '-n', default=1, help='The number of solution to be checked in the specified file')
parser.add_argument('--samples', '-s', default=1, help='The number of samples to start with')
args = parser.parse_args()


# Classes

class Sequence:
    """ Sequence class

    it contains the protein sequence and calculates statistics like
    occurrence of each unique pair of residues, ranking of residues by number
    of unique pairs they have and creating the list of residues to be labeled
    with 'calculate stats' method
    based on amino acid stock, written in Stock class
    """

    def __init__(self, task, name="", sequence=""):
        self.name = name
        self.sequence = sequence # sequence itself
        self.got = True
        self.task = task  # task object
        self.bad_residues = []              # residues that don't have any CN labelings and
                                            # have both diagonal and 'Other' element in its column
        self.nitro_types = ("N", "D", "S", "T")         # labeling types with 15N
        self.carbon_types = ("C", "D", "A", "S", "T", "F")      # labeling types with 13C
        self.res_has_diagonal = {}    # dict of residues. If residue has diagonal element, value = True
        self.residue_pairs = [[]]

    def clean_sequence(self, lines):
        # makes all upper case letters in sequence
        # and clears up all letters, that are not one-letter AA codes
        joined = "".join(lines)
        sequence = ""
        for char in joined.upper():
            if char in self.task.res_types:
                sequence += char
        return sequence

    def read_from_file(self, sequence_lines):
        new_lines = []
        # dealing with FASTA files
        if sequence_lines[0][0] == ">":
            for line in sequence_lines[1:]:
                if line[0] == ">":
                    break
                elif line[0] == "#":
                    continue
                new_lines.append(line)
        else:
            new_lines = sequence_lines
        sequence = self.clean_sequence(new_lines)
        self.sequence = sequence
        if len(sequence) < 2:
            raise err.ReadConfigError(msg="ERROR! Sequence too short or absent")
        self.got = True

    def write_to_file(self, filename, mode='w'):
        ## ADD check if file exists

        f = open(filename, mode)
        f.write(self.sequence+"\n")
        f.flush()
        f.close()

    def calculate_stats(self, stock):
        # Stock class object required for calculating stats
        # Almost all methods should be rewritten as some of class variables
        # are defined and created in the methods other than '__init__'

        self.stock = stock
        self._rank_residues()       # rank all residues in sequence by the number of
                                    # distinct pairs (two or more equal pairs are counted as one)
        self._residues_to_label()   # creating a list of residues that are to be labeled based
                                    # on the residue ranks and presence of labels in stock
        self._count_pairs()         # counting all the pairs present in sequence
        self._check_carbon_pairs()  # recount pairs in order to skip empty rows
        self._count_nitrogens()     # counting non-epty cells in columns
                                    # for residues, that have 15N-labelin in stock
        self._find_bad_residues()   # find residues, that have diagonal pair and pair with
                                    # "Other" residues and don't have C/N labels in stock
                                    # so that you can't tell diagonal cell from pair with "Other"
        self._calculate_subtable_coordinates()      # precalculate coordinates of
                                                    # subtable cells for each residue
        self._calculate_under_cells()   # find cells that are under the particular one for crossing-out
        self._calculate_right_cells()   # NOT USED find cells that are
                                        # under the particular one for crossing-out

    def count_residues(self):
        residue_count = {}
        for residue in self.task.res_types:
            residue_count.update({residue : 0})
        for i in range(len(self.sequence)):
            residue_count[self.sequence[i]] += 1
        return residue_count

    def _rank_residues(self):
        # Bad code for initializing outside __init__

        residue_rank = [0 for _ in self.task.res_types]
        self.unique_pairs = []
        self.res_no_first = list(self.task.res_types)
        self.res_no_second = list(self.task.res_types)
        for i in range(len(self.sequence) - 1):
            pair = self.sequence[i:i+2]
            if pair[0] in self.task.res_types and pair[1] in self.task.res_types:
                if pair not in self.unique_pairs:
                    self.unique_pairs.append(pair)
                    if pair[1] != pair[0]:
                        residue_rank[self.task.res_types.index(pair[1])] += 1
                    residue_rank[self.task.res_types.index(pair[0])] += 1
                if pair[0] in self.res_no_first:
                    self.res_no_first.pop(self.res_no_first.index(pair[0]))
                if pair[1] in self.res_no_second:
                    self.res_no_second.pop(self.res_no_second.index(pair[1]))
        self.ranked_residues, self.rank_of_residue = self._sort_residues(self.task.res_types, residue_rank)
        self.res_no_second.append("P")

        # residues that are first in some pair in sequence
        self.residues_first = [res for res in self.ranked_residues if res not in self.res_no_first]

        # residues that are second in some pair in sequence
        self.residues_second = [res for res in self.ranked_residues if res not in self.res_no_second]

    def _residues_to_label(self):
        self.residues_nitro = []        # list of residues with 15N labels. Last one can be "Other"
        self.residues_not_nitro = []    # list of residues without 15N labels
        self.residues_carbon = []       # list of residues with 13C labels. Last one can be "Other"
        self.residues_not_carbon = []   # list of residues without 13C labels
        self.residues_to_label = []     # ordered list of residues with 13C or 15N labels
        self.non_labeled_residues = []  # residues without any labels

        for i in range(len(self.ranked_residues)):
            residue = self.ranked_residues[i]
            if self.rank_of_residue[i]:
                nitro = False
                for type in self.nitro_types:
                    if type in self.stock.label_options[residue]:
                        nitro = True
                        break
                if nitro and residue in self.residues_second:
                    self.residues_nitro.append(residue)
                else:
                    self.residues_not_nitro.append(residue)
                carbon = False
                for type in self.carbon_types:
                    if type in self.stock.label_options[residue]:
                        carbon = True
                        break
                if carbon and residue in self.residues_first:
                    self.residues_carbon.append(residue)
                else:
                    self.residues_not_carbon.append(residue)
                if residue in self.residues_carbon or residue in self.residues_nitro:
                    self.residues_to_label.append(residue)
                    self.res_has_diagonal[residue] = False
                else:
                    self.non_labeled_residues.append(residue)
        if len(self.residues_not_nitro):        # add "Other" column if there
            self.residues_nitro.append("Other") # are residues not labeled by 15N
        if len(self.residues_not_carbon):       # add "Other" row if there
            self.residues_carbon.append("Other")# are residues not labeled by 13C

    def _sort_residues(self, residue_list, residue_rank):
        # bubble sort for residues by residue rank.
        # used it because standard sorted() method gives randomized results
        residues = list(residue_list)
        for i in range(len(residues)-1):
            for j in range(len(residues)-1-i):
                if residue_rank[i] < residue_rank[i+j+1]:
                    temp_res = residues[i]
                    temp_rank = residue_rank[i]
                    residue_rank[i] = residue_rank[i + j + 1]
                    residues[i] = residues[i + j + 1]
                    residue_rank[i + j + 1] = temp_rank
                    residues[i + j + 1] = temp_res
        return residues, residue_rank

    def _check_carbon_pairs(self):
        carb_other = 0
        changed = False
        if self.residues_carbon and self.residues_carbon[-1] == "Other":
            carb_other = 1
        if self.residues_nitro and self.residues_nitro[-1] == "Other":
            carbons = len(self.residues_carbon) - carb_other
            for i in range(carbons):
                k = carbons - 1 - i
                got_pair = False
                for j in range(len(self.residues_nitro) - 1):
                    if self.residue_pairs[k][j]:
                        got_pair = True
                        break
                if not got_pair:
                    self.residues_not_carbon.append(self.residues_carbon[k])
                    residue = self.residues_carbon[k]
                    del self.residues_carbon[k]
                    if residue not in self.residues_nitro:
                        del self.residues_to_label[self.residues_to_label.index(residue)]
                    changed = True
        self.ranks = {}
        for res in self.residues_to_label:
            self.ranks.update({res: 0})
        if changed:
            self._count_pairs()
        for i in range(len(self.residues_nitro)):
            if self.residues_nitro[i] == "Other":
                break
            for j in range(len(self.residues_carbon)):
                if self.residue_pairs[j][i]:
                    if self.residues_carbon[j] != "Other":
                        self.ranks[self.residues_carbon[j]] += 1
                    if self.residues_carbon[j] != self.residues_nitro[i]:
                        self.ranks[self.residues_nitro[i]] += 1
        ranks = []
        residues = []
        for residue in self.residues_to_label:
            if residue in self.ranks:
                residues.append(residue)
                ranks.append(self.ranks[residue])
        self.residues_to_label, self.residue_ranks = self._sort_residues(residues, ranks)
        self.residues_nitro_new = []
        self.residues_carbon_new = []
        for res in self.residues_to_label:
            if res in self.residues_carbon:
                self.residues_carbon_new.append(res)
            if res in self.residues_nitro:
                self.residues_nitro_new.append(res)
        if "Other" in self.residues_carbon:
            self.residues_carbon_new.append("Other")
        if "Other" in self.residues_nitro:
            self.residues_nitro_new.append("Other")
        self.residues_carbon = self.residues_carbon_new
        self.residues_nitro = self.residues_nitro_new
        self._count_pairs()

    def _add_rank(self, residue):
        self.ranks[self.residues_to_label.index(residue)] += 1

    def _count_nitrogens(self):
        # counting non-empty cells in columns
        # for residues, that have 15N-labeling in stock
        self.min_nitrogens = []
        for i in range(len(self.residues_nitro)):
            pairs_count = 0
            for j in range(len(self.residues_carbon)):
                if self.residue_pairs[j][i]:
                    pairs_count += 1
            self.min_nitrogens.append(pairs_count)

    def _find_bad_residues(self):
        # find residues, that have diagonal pair and pair with
        # "Other" residues and don't have C/N labels in stock
        # so that you can't tell diagonal cell code from code of pair with "Other"

        for i in range(len(self.residues_nitro)):
            res = self.residues_nitro[i]
            if res == "Other":
                continue
            if ("D" not in self.stock.label_options[res]
                  and "S" not in self.stock.label_options[res]
                  and "T" not in self.stock.label_options[res]
                  and self.residues_carbon[-1] == 'Other'
                  and self.residues_nitro[i] in self.residues_carbon
                  and self.residue_pairs[self.residues_carbon.index(res)][i]
                  and self.residue_pairs[-1][i]):
                self.bad_residues.append(res)

    def _count_pairs(self):
        # table of residue pairs taken in account labeling stock
        self.residue_pairs = [[0 for col in range(len(self.residues_nitro))]
                              for row in range(len(self.residues_carbon))]
        # table of all residue pairs present in protein sequence
        self.all_residue_pairs = [[0 for col in range(len(self.residues_second))]
                                  for row in range(len(self.residues_first))]
        for i in range(len(self.sequence) - 1):  # Count all pairs in sequence,
            pair = self.sequence[i:(i + 2)]
            if pair[0] == '*' or pair[1] == '*':
                continue
            if pair[0] not in self.residues_carbon:
                first_res = len(self.residues_carbon) - 1
            else:
                first_res = self.residues_carbon.index(pair[0])
            if pair[1] not in self.residues_nitro:
                second_res = len(self.residues_nitro) - 1
            else:
                second_res = self.residues_nitro.index(pair[1])
            self.residue_pairs[first_res][second_res] += 1
            if pair[0] == pair[1]:
                self.res_has_diagonal[pair[0]] = True
            index_first = self.residues_first.index(pair[0])
            if pair[1] != "P":
                index_second = self.residues_second.index(pair[1])
                self.all_residue_pairs[index_first][index_second] += 1

    def _calculate_subtable_coordinates(self):
        # precalculate coordinates of
        # subtable cells for each residue

        self.subtable_coordinates = []
        check_other = (self.residues_carbon and self.residues_carbon[-1] == 'Other')
        carbon = 0
        nitro = 0
        total_cells = 0

        for i in range(len(self.residues_to_label)):
            self.new_coordinates = []
            residue = self.residues_to_label[i]
            if i == 0:
                if residue in self.residues_nitro:
                    if residue in self.residues_carbon:
                        if self.residue_pairs[0][0]:
                            self.new_coordinates.append((0, 0))
                        carbon += 1
                    if check_other and self.residue_pairs[-1][0] and residue not in self.bad_residues:
                        self.new_coordinates.append((len(self.residues_carbon) - 1, 0))
                    nitro += 1
            else:
                if residue in self.residues_carbon:
                    for j in range(nitro):
                        if self.residue_pairs[carbon][j]:
                            self.new_coordinates.append((carbon, j))
                    carbon += 1
                if residue in self.residues_nitro:
                    for j in range(carbon):
                        if self.residue_pairs[j][nitro]:
                            self.new_coordinates.append((j, nitro))
                    if check_other and self.residue_pairs[-1][nitro] and residue not in self.bad_residues:
                        self.new_coordinates.append((len(self.residues_carbon) - 1, nitro))
                    nitro += 1
            total_cells += len(self.new_coordinates)
            self.subtable_coordinates.append(self.new_coordinates)
        curr_cells = 0
        for i in range(len(self.subtable_coordinates)):
            add = len(self.subtable_coordinates[i])
            curr_cells += add

    def _calculate_under_cells(self):
        # find cells that are under the particular one for crossing-out

        self.under_cells = {}
        got_other = 0
        if self.residues_carbon and self.residues_carbon[-1] == "Other":
            got_other = 1
        rows = len(self.residues_carbon) - got_other
        for cell_list in self.subtable_coordinates:
            for cell in cell_list:
                under_cells = []
                row = cell[0] + 1
                if row > rows:
                    row = 0     # if the cell is in "Other" row,
                                #  then all other cells in the column are under it
                column = cell[1]
                for i in range(rows - row):
                    if self.residue_pairs[row+i][column]:
                        under_cells.append((row+i, column))
                self.under_cells.update({cell: under_cells})

    def _calculate_right_cells(self):
        # NOT USED, unless a new smart method for crossing out is designed
        # almost duplicate of "_calculate_under_cells"

        self.right_cells = {}
        got_other = 0
        if self.residues_nitro and self.residues_nitro[-1] == "Other":
            got_other = 1
        columns = len(self.residues_nitro) - got_other
        for cell_list in self.subtable_coordinates:
            for cell in cell_list:
                right_cells = []
                row = cell[0]
                column = cell[1] + 1
                for i in range(columns - column):
                    if self.residue_pairs[row][column+i]:
                        right_cells.append((row, column+i))
                self.right_cells.update({cell: right_cells})

    def _need_carbon(self, residue):
        carbon_number = self.residues_carbon.index(residue)
        if self.residues_carbon[-1] == "Other":
            scan_nitro = len(self.residues_nitro)
            if self.residues_nitro[-1] == "Other":
                scan_nitro -= 1
            if residue in self.residues_nitro:
                nitro_index = self.residues_nitro.index(residue)
                if self.residue_pairs[carbon_number][nitro_index] and self.residue_pairs[-1][nitro_index] and residue not in self.bad_residues:
                    return self._find_cheapest_option(residue, self.task.cn_label_types)
            for i in range(scan_nitro):
                if self.residue_pairs[carbon_number][i] and self.residue_pairs[-1][i]:
                    if residue == self.residues_nitro[i]:
                        if residue in self.bad_residues:
                            continue
                        else:
                            return self._find_cheapest_option(residue, self.task.cn_label_types)
                    else:
                        return self._find_cheapest_option(residue, self.carbon_types)
        return self._find_cheapest_option(residue, self.task.res_types)

    def _find_cheapest_option(self, residue, labeling_options):
        options = self.stock.label_options[residue]
        first_option = True
        lower_price = 0
        cheapest_option = ""
        for i in range(len(options)):
            option = options[i]
            if option in labeling_options:
                curr_price = self.stock.price_dict[residue][option]
                if first_option:
                    lower_price = curr_price
                    cheapest_option = option
                    first_option = False
                elif curr_price < lower_price:
                    cheapest_option = option
                    lower_price = curr_price
                else:
                    continue
        return cheapest_option

    def _need_nitrogen(self, residue, res_options):
        min_nitrogens = self.min_nitrogens[self.residues_nitro.index(residue)]
        all_options = self.stock.label_options[residue]
        label_powers = self.task.coding_table.label_power
        nitro_options = []
        self.variants = []
        current_power = 1
        for option in all_options:
            if option in self.nitro_types:
                nitro_options.append(option)
        for option in res_options:
            if option in self.nitro_types:
                current_power *= label_powers[option]
        if (min_nitrogens == 1 and current_power > 1) or (min_nitrogens > 1 and current_power >= min_nitrogens):
            return []
        else:
            self._add_nitrogen([], nitro_options, residue, min_nitrogens, current_power)
            best_result = []
            first_iteration = True
            for variant in self.variants:
                if first_iteration:
                    first_iteration = False
                    best_result = variant
                    best_price = self._calc_price_for_variant(residue, variant)
                else:
                    if self._calc_price_for_variant(residue, variant) < best_price:
                        best_result = variant
                        best_price = self._calc_price_for_variant(residue, variant)
            return best_result

    def _calc_price_for_variant(self, residue, variant):
        price = 0
        for type in variant:
            price += self.stock.price_dict[residue][variant]
        return price

    def _add_nitrogen(self, list1, nitro_options, residue, min_nitrogens, current_power):
        for option in nitro_options:
            list1.append(option)
            if self._check_power(list1, min_nitrogens, current_power):
                self.variants.append(list1)
            else:
                self._add_nitrogen(list1, nitro_options, residue, min_nitrogens, current_power)
            list1.pop()

    def _check_power(self, list1, min_nitrogens, current_power):
        power = current_power
        for option in list1:
            power *= self.task.coding_table.label_power[option]
        return power > min_nitrogens

    def _calculate_base_price(self, residue, options):
        price = 0
        for option in options:
            price += self.stock.price_dict[residue][option]
        return price

    def predict_prices(self):
        # predict prices for old algorithm
        # this function finds cheapest possible labeling for each residue

        residues = self.residues_to_label
        self.required_options = []
        self.base_prices = []
        self.add_prices = []
        self.max_requirement = 0
        for res in residues:
            res_options = []
            if res in self.residues_nitro:
                if res in self.residues_carbon:
                    res_options.append(self._need_carbon(res))
                nitro_required = self._need_nitrogen(res, res_options)
                for nitro in nitro_required:
                    res_options.append(nitro)
            else:
                res_options.append(self._need_carbon(res))
            self.required_options.append(res_options)
            if len(res_options) > self.max_requirement:
                self.max_requirement = len(res_options)
        for i in range(len(residues)):
            res = residues[i]
            cheapest_option = self._find_cheapest_option(res, self.task.res_types)
            add_price = self.stock.price_dict[res][cheapest_option]
            base_price = self._calculate_base_price(res, self.required_options[i])
            for j in range(self.max_requirement - len(self.required_options[i])):
                base_price += add_price
            self.base_prices.append(base_price)
            self.add_prices.append(add_price)
        return 0

    def __str__(self):
        return self.sequence


class Solution:
    '''
    Solution class

    Requires Sequence and Stock objects for initialization

    Has 2 main methods:
        1) Increment state: for each state this method modifies
            self.solution variable depending on whether all the checks
            are met
        2) Check_and_update: each solution is checked in each state

    All other methods are supportive for those two
    '''

    def __init__(self, task):
        self.task = task
        self.calc_code = 0
        self.name = self.task.name
        self.sequence = self.task.sequence
        self.coding_table = self.task.coding_table
        self.optimize_price = self.task.optimize_price
        if self.optimize_price:
            self.base_prices = self.sequence.base_prices
            self.add_prices = self.sequence.add_prices
            self.max_requirement = self.sequence.max_requirement
            self.predicted_prices = [0 for res in self.sequence.residues_to_label]
        self.solution = []
        self.samples_num = task.start_samples
        self.found = False
        self.good = False
        self.stock = self.task.stock
        self.symmetry = []
        self.depth = 0
        self.max_depth = 0
        self.parallel_run=True
        self.parallel_depth=3
        self.price = 0
        self.nitro = 0
        self.carbon = 0
        self.check_pattern_and_pairs = 0
        self.compare_to_previous_code_patterns = 0
        self.check_symmetry = 0
        self.check_subtable = 0
        self.time_increment = 0
        self.time_subtable = 0
        self.time_check = 0
        self.time_nitrogens = 0
        self.time_symmetry = 0
        self.time_patterns = 0
        self.current_patterns = []
        self.subtable = []
        self.got_other = 0
        self.check_other = (self.sequence.residues_carbon and self.sequence.residues_carbon[-1] == 'Other')
        if self.check_other:
            self.got_other = 1
        self.labeling_dictionary = {}
        self.best_price = 0
        self.last_price = []
        self.usage = self.stock.usage
#        self.rt = rs.make_rt()
        self.codes_dict = self.coding_table.codes_dict
        self.label_power = self.coding_table.label_power
        self.solution_number = 0
        self.stats = {}

    def __str__(self):
        output = ""
        if self.optimize_price:
            price = self.calculate_total_price()
            output += "Current state price: " + str(price) + "\n"
        output += "Number of samples:" + str(self.samples_num)+"\n"
        output += "Depth:" + str(self.depth) + "\n"
        for i in range(len(self.solution)):
            res = self.task.to_three_letter_code[self.sequence.residues_to_label[i]]
            output += res + "," + self.solution[i] + "\n"
        for res in self.sequence.non_labeled_residues:
            output += self.task.to_three_letter_code[res] + "," + self._generate_other_code() + "\n"
        output += "Codes calculated:" + str(self.calc_code)
        return output

    def copy_solution_to(self, solution):
        solution.sequence = self.sequence
        solution.labeling_dictionary = copy.deepcopy(self.labeling_dictionary)
        solution.samples_num = copy.copy(self.samples_num)
        solution.solution = list(self.solution)
        solution.price = copy.copy(self.price)
        solution.depth = copy.copy(self.depth)
        solution.found = copy.copy(self.found)
        solution.solution_number = copy.copy(self.solution_number)

    def calculate_total_price(self):
        price = 0
        for i in range(self.depth):
            residue = self.sequence.residues_to_label[i]
            labeling = self.solution[i]
            price += self.calculate_price_for_residue(residue, labeling)
        return price

    def read_existing_solution(self, filename):
        # rewrite, so that it can read self-made solutions
        self.labeling_dictionary = {}
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            d = list(reader)

            # ADD check if file format is correct

        for i in range(len(d)):
            self.labeling_dictionary.update({d[i][0]: d[i][1]})
        self.samples_num = len(d[0][1])
        self.labeling_dictionary.update({"Other": self._generate_other_code()})
        self.solution = []
        for residue in self.sequence.residues_to_label:
            self.solution.append(self.labeling_dictionary[residue])

    def increment_state(self):
        if self.solution == []:
            # samples ++; depth = 0
            self.samples_num += 1
            self.labeling_sets_dict = {}
            self.labeling_sets = self.generate_labeling_sets()
            for labeling_set in self.labeling_sets:
                self.labeling_sets_dict[labeling_set.residue] = labeling_set
            if not self.found:
                output = ""
                if self.samples_num > self.task.start_samples + 1:
                    output += "Max depth reached: {}\n".format(self.max_depth)
                else:
                    output += "\n"
                output += "Search in {} samples\n".format(self.samples_num)

                print(output)
                self.task.logfile.write(output + "\n")
                self.task.logfile.flush()

            # still actual until new prediction algorithm
            if self.optimize_price:
                if self.samples_num >= self.max_requirement:
                    self.predicted_prices = [(self.samples_num - self.max_requirement) * x + y for x, y in
                                             zip(self.add_prices, self.base_prices)]
                else:
                    self.predicted_prices = self.base_prices

            self.other_code = self._generate_other_code()
            self.labeling_dictionary.update({'Other': self.other_code})

            self.labeling_sets[0].cross_out_symmetry([self.samples_num])

            self.depth += 1
            self.labeling_sets[0].update_on_depth_change()
            if self.labeling_sets[0].is_empty():
                self.depth -= 1
                self.increment_state()
            else:
                self.solution.append(self.labeling_sets[0].get_labeling())

                if self.depth > self.max_depth:
                    self.max_depth = self.depth
                if self.sequence.residues_to_label[self.depth - 1] in self.sequence.residues_nitro:
                    self.nitro += 1
                if self.sequence.residues_to_label[self.depth - 1] in self.sequence.residues_carbon:
                    self.carbon += 1
                self.symmetry.append([self.samples_num])

        elif self.good and self.depth != len(self.sequence.residues_to_label):
            # depth ++
            self.symmetry.append(self._calculate_symmetry())
            last_residue = self.sequence.residues_to_label[self.depth - 1]
            last_label = self.solution[-1]
            self.labeling_dictionary.update({last_residue: last_label})
            self.subtable.append(set(self.new_codes))
            if self.optimize_price:
                self.last_price.append(self.calculate_price_for_residue(last_residue,
                                                                        last_label))
                self.price += self.last_price[-1]
            self.labeling_sets[self.depth].cross_out_symmetry(self.symmetry[-1])
            self.depth += 1
            if self.depth > self.max_depth:
                self.max_depth = self.depth
                #print("max depth: ", self.max_depth)

            # if last_residue in self.sequence.residues_nitro:
            #     self.current_patterns.append(self._convert_code_into_nitrogen_pattern(self.solution[-1]))

            self.labeling_sets[self.depth - 1].update_on_depth_change()
            self.solution.append(self.labeling_sets[self.depth - 1].get_labeling())
            if self.sequence.residues_to_label[self.depth - 1] in self.sequence.residues_nitro:
                self.nitro += 1
            if self.sequence.residues_to_label[self.depth - 1] in self.sequence.residues_carbon:
                self.carbon += 1

        elif self.labeling_sets[self.depth-1].has_last_index():
            # depth --
            if self.sequence.residues_to_label[self.depth - 1] in self.sequence.residues_nitro:
                self.nitro -= 1
                if self.current_patterns != []:
                    self.current_patterns.pop()
            if self.sequence.residues_to_label[self.depth - 1] in self.sequence.residues_carbon:
                self.carbon -= 1

            # add restoration of labeling sets
            for labeling_set in self.labeling_sets[self.depth:]:
                labeling_set.restore_last_depth()

            self.depth -= 1
            self.labeling_sets[self.depth].restore_symmetry()

            self.solution.pop()
            self.symmetry.pop()
            if self.optimize_price and self.last_price != []:
                self.price -= self.last_price[-1]
                self.last_price.pop()
            if self.subtable != []:
                self.subtable.pop()
            self.increment_state()
        else:
            # next index in the same depth
            self.labeling_sets[self.depth-1].increment_index()
            self.solution[-1] = self.labeling_sets[self.depth-1].get_labeling()
            for labeling_set in self.labeling_sets[self.depth:]:
                labeling_set.restore_last_depth()

    def generate_labeling_sets(self):
        list_of_labeling_sets = []
        for i in range(len(self.sequence.residues_to_label)):
            labeling_set = LabelList(self, i)
            list_of_labeling_sets.append(labeling_set)
        return list_of_labeling_sets

    def _increment_code(self, current_code, options_to_label):
        if len(current_code) == 0:
            return current_code
        if current_code[-1] == options_to_label[-1]:
            return ''.join([self._increment_code(current_code[:-1], options_to_label), options_to_label[0]])
        else:
            return ''.join([current_code[:-1],
                           options_to_label[options_to_label.index(current_code[-1]) + 1]])

    def _base_case(self):
        residue = self.sequence.residues_to_label[len(self.solution)]
        return ''.join([self.stock.label_options[residue][0] for _ in range(self.samples_num)])

    def _calculate_symmetry(self):
        symmetry = self.symmetry[-1]
        code = self.solution[-1]
        if len(symmetry) == self.samples_num:
            return symmetry
        new_symmetry = []
        start_point = 0
        for i in range(len(symmetry)):
            if symmetry[i] > 1:
                number = 1
                for j in range(symmetry[i] - 1):
                    if code[start_point + j] == code[start_point + j + 1]:
                        number += 1
                    else:
                        new_symmetry.append(number)
                        number = 1
                    if j == symmetry[i] - 2:
                        new_symmetry.append(number)
            else:
                new_symmetry.append(1)
            start_point += symmetry[i]
        return new_symmetry

    def check_and_update(self):
        result = (#self._check_symmetry()
                  self._check_subtable()
                  and self._cross_out()
                  and self._check_price_2()
                  )
        self.good = result
        return result

    def _cross_out(self):
        coordinates_list = self.sequence.subtable_coordinates[self.depth - 1]

        for i in range(len(coordinates_list)):
            cell = coordinates_list[i]
            code = self.new_codes[i]

            under_cells = self.sequence.under_cells[cell]

            for current_cell in under_cells:
                if current_cell[1] + 1 == self.nitro and current_cell[0] < self.carbon:
                    continue
                second_residue = self.sequence.residues_nitro[current_cell[1]]
                first_residue = self.sequence.residues_carbon[current_cell[0]]
                first_labeling_set = self.labeling_sets_dict[first_residue]
                second_labeling = self.labeling_dictionary[second_residue]
                cross_out_set = self._make_cross_out_set(first_labeling_set, second_labeling, code)
                first_labeling_set.cross_out(cross_out_set)
                if not first_labeling_set.label_set:
                    return False

            # right_cells = self.sequence.right_cells[cell]
            #
            # for current_cell in right_cells:
            #     if current_cell[0] + 1 == self.carbon and current_cell[1] < self.nitro:
            #         continue
            #     second_residue = self.sequence.residues_nitro[current_cell[1]]
            #     first_residue = self.sequence.residues_carbon[current_cell[0]]
            #     second_labeling_set = self.labeling_sets_dict[second_residue]
            #     first_labeling = self.labeling_dictionary[first_residue]
            #     cross_out_set = self._make_cross_out_set(second_labeling_set, first_labeling, code, under=False)
            #     second_labeling_set.cross_out(cross_out_set)
            #     if not second_labeling_set.label_set:
            #         return False
        return True

    def _make_cross_out_set(self, labeling_set, labeling, code, under=True):
        cross_out_set = set()
        second_labeling = labeling
        first_labeling = labeling
        for curr_labeling in labeling_set.label_set:
            if under:
                first_labeling = curr_labeling
                if code == self._calculate_code_multiple(first_labeling, second_labeling):
                    cross_out_set.add(first_labeling)
            else:
                second_labeling = curr_labeling
                second_residue = labeling_set.residue
                if self.sequence.res_has_diagonal[second_residue] and self._calculate_code_multiple(first_labeling, second_labeling) == self._calculate_code_multiple(second_labeling, second_labeling):
                    cross_out_set.add(first_labeling)
        return cross_out_set

    def _check_subtable(self):
        self.new_codes = []
        self.labeling_dictionary.update({self.sequence.residues_to_label[self.depth - 1]: self.solution[-1]})
        coordinates_list = self.sequence.subtable_coordinates[self.depth - 1]

        for coordinates in coordinates_list:
            first_residue = self.sequence.residues_carbon[coordinates[0]]
            second_residue = self.sequence.residues_nitro[coordinates[1]]
            first_residue_labeling = self.labeling_dictionary[first_residue]
            second_residue_labeling = self.labeling_dictionary[second_residue]
            code = self._calculate_code_multiple(first_residue_labeling,
                                                 second_residue_labeling)
            if code == "0" * self.samples_num:
                return False
            if self._code_used(code):
                return False
            self.new_codes.append(code)
        return True

    def _code_used(self, new_code):
        if new_code in self.new_codes:
            return True
        for code_list in self.subtable:
            if new_code in code_list:
                return True
        return False

    def _check_symmetry(self):
        if len(self.symmetry[-1]) == self.samples_num:
            return True
        symmetry = self.symmetry[-1]
        code = self.solution[-1]
        residue = self.sequence.residues_to_label[self.depth - 1]
        label_num_dict = self.stock.label_num_dict[residue]
        start_point = 0
        for i in range(len(symmetry)):
            if symmetry[i] > 1:
                for j in range(symmetry[i] - 1):
                    if (label_num_dict[code[start_point + j]]
                            > label_num_dict[code[start_point + j + 1]]):
                        return False
            start_point += symmetry[i]
        return True

    def _check_price(self):
        last_residue = self.sequence.residues_to_label[self.depth-1]
        last_label = self.solution[-1]
        predicted_price = 0
        for i in range(len(self.solution) - self.depth):
            predicted_price += self.predicted_prices[self.depth+i]

        price = (self.price
                 + self.calculate_price_for_residue(last_residue, last_label)
                 + predicted_price
                 )
        if price >= self.best_price:
            return False
        return True

    def _check_price_2(self):
        if not self.optimize_price or not self.found:
            return True
        last_residue = self.sequence.residues_to_label[self.depth-1]
        last_label = self.solution[-1]
        # predicted_price = 0
        # for i in range(len(self.solution) - self.depth):
        #     predicted_price += self.predicted_prices[self.depth+i]

        price = (self.price
                 + self.calculate_price_for_residue(last_residue, last_label)
                 + self._predict_price_for_sets()
                 )
        if price >= self.best_price:
            return False
        return True

    def _predict_price_for_sets(self):
        price = 0
        for label_set in self.labeling_sets[self.depth:]:
            price += label_set.predict_price()
        return price

    def _find_cheapest_option(self, residue):
        first = True
        cheapest = 0
        for option in self.stock.label_options[residue]:
            if first:
                cheapest = int(self.stock.price_dict[residue][option])
                first = False
            elif int(self.stock.price_dict[residue][option]) < cheapest:
                cheapest = int(self.stock.price_dict[residue][option])
        return cheapest

    def _calculate_code_multiple_2(self, first_aa_code, second_aa_code):
        samples = self.samples_num
        # code_string = ''.join([str(self._calculate_code(first_aa_code[i],
        #                                                 second_aa_code[i])) for i in range(samples)])
        code_string = ''.join([self._calculate_code(first_aa_code[i],
                                                        second_aa_code[i]) for i in range(samples)])
        return code_string

    def _calculate_code_multiple(self, first_aa_code, second_aa_code):
        self.calc_code += 1
        return ''.join(self.codes_dict[first_aa_code[i]][second_aa_code[i]] for i in range(self.samples_num))

    def _calculate_code(self, first_aa_code, second_aa_code):
        return self.codes_dict[first_aa_code][second_aa_code]

    def _compare_nitrogen_patterns(self, code1, code2):
        return (self._convert_code_into_nitrogen_pattern(code1)
                == self._convert_code_into_nitrogen_pattern(code2))

    def _convert_code_into_nitrogen_pattern(self, code):
        new_code = ''
        for i in range(len(code)):
            if code[i] not in self.task.n_label_types:
                new_code += '0'
            else:
                new_code += '1'
        return new_code

    def calculate_price_for_residue(self, residue, labeling):
        price = 0
        for i in range(len(labeling)):
            price += self.stock.price_dict[residue][labeling[i]]
        return price

    def _generate_last_solution(self):
        self.last_solution = []
        for residue in self.sequence.residues_to_label:
            current_code = ''.join([self.stock.label_options[residue][-1] for _ in range(self.samples_num)])
            self.last_solution.append(current_code)

    def _generate_other_code(self):
        code = ''
        for i in range(self.samples_num):
            code += 'X'
        return code

    def print_table(self, codes=True):
        aa_types_1 = self.sequence.residues_carbon
        aa_types_2 = self.sequence.residues_nitro
        aa_pairs = self.sequence.residue_pairs
        code_dict = self.labeling_dictionary
        num_of_samples = self.samples_num
        other_nitro = 0
        if aa_types_2[-1] == "Other":
            other_nitro = 1

        cell_height = 20
        cell_width = 50

        row_number = len(aa_types_1)
        column_number = len(aa_types_2) - other_nitro

        top = tk.Tk()
        top.title("CombLabel solution for {}".format(self.name))
        C = tk.Canvas(top, bg="black", height=cell_height * (row_number + 2),
                           width=cell_width * (column_number + 2))

        ### Draw the table and fill the headers
        for i in range(row_number):
            res = aa_types_1[i]
            if aa_types_1[i] != "Other":
                res = self.task.to_three_letter_code[aa_types_1[i]]
            C.create_line(0, cell_height * (i + 1), cell_width * (column_number + 2), cell_height * (i + 1),
                          fill="white")
            C.create_text(cell_width / 2, cell_height * (i + 2.5), text=res, fill="white")
            C.create_text(cell_width * 1.5, cell_height * (i + 2.5), text=code_dict[aa_types_1[i]], fill="white")
        for i in range(column_number):
            C.create_line(cell_width * (i + 1), 0, cell_width * (i + 1), cell_height * (row_number + 2),
                          fill="white")
            res = aa_types_2[i]
            if aa_types_2[i] != "Other":
                res = self.task.to_three_letter_code[aa_types_2[i]]
            C.create_text(cell_width * (i + 2.5), cell_height / 2, text=res, fill="white")

            column_codes = self.get_codes(aa_types_2, code_dict)
            codes_match = 0
            for k in range(column_number):
                if (k != i and self._compare_nitrogen_patterns(column_codes[i], column_codes[k])):
                    codes_match = 1
                if codes_match:
                    cell_color = 'red'
                else:
                    cell_color = 'green'
            C.create_rectangle(cell_width * (i + 2) + 1, cell_height + 1, cell_width * (i + 3) - 1,
                               cell_height * 2 - 1, fill=cell_color)
            C.create_text(cell_width * (i + 2.5), cell_height * 1.5, text=code_dict[aa_types_2[i]], fill="white")
        C.create_line(0, cell_height * (row_number + 1), cell_width * (column_number + 2),
                      cell_height * (row_number + 1), fill="white")
        C.create_line(cell_width * (column_number + 1), 0, cell_width * (column_number + 1),
                      cell_height * (row_number + 2), fill="white")

        all_aa_codes = [['' for col in range(column_number)] for row in range(row_number)]
        for i in range(row_number):
            for j in range(column_number):
                if aa_pairs[i][j]:
                    all_aa_codes[i][j] = self._calculate_code_multiple(code_dict[aa_types_1[i]], code_dict[aa_types_2[j]])

        ### fill the cells
        good_cells = 0
        bad_cells = 0
        for i in range(row_number):
            for j in range(column_number):
                if aa_pairs[i][j] > 0:
                    if codes:
                        codes_match = 0
                        for k in range(row_number):
                            if (k != i and all_aa_codes[i][j] == all_aa_codes[k][j]):
                                codes_match = 1
                        if codes_match:
                            cell_color = 'red'
                            bad_cells += 1
                        else:
                            cell_color = 'green'
                            good_cells += 1
                    else:
                        cell_color = 'green'
                        good_cells += 1
                    C.create_rectangle(cell_width * (j + 2) + 1, cell_height * (i + 2) + 1,
                                       cell_width * (j + 3) - 1, cell_height * (i + 3) - 1, fill=cell_color)
                    if codes:
                        C.create_text(cell_width * (j + 2.5), cell_height * (i + 2.5), text=all_aa_codes[i][j],
                                      fill="white")
                    else:
                        C.create_text(cell_width * (j + 2.5), cell_height * (i + 2.5), text = aa_pairs[i][j], fill="white")

        C.create_rectangle(0, 0, cell_width * 2 - 1, cell_height * 2 - 1, fill="black")
        # print(str(good_cells / (good_cells + bad_cells) * 100) + '%', (good_cells + bad_cells))

        C.pack()
        top.mainloop()

    def get_codes(self, aa_string, code_dict):
        aa_codes = []
        for i in range(len(aa_string)):
            aa_codes.append(code_dict[aa_string[i]])
        return aa_codes

    def write_codes_table(self, filename, mode='w'):
        self.codes_table = []
        self.meanings_table = []
        self.code_meanings = {}
        for res in self.sequence.non_labeled_residues:
            self.labeling_dictionary[res] = self._generate_other_code()
        for i in range(len(self.sequence.sequence) - 1):
            pair = self.sequence.sequence[i:i+2]
            residue1 = pair[0]
            residue2 = pair[1]
            scheme1 = self.labeling_dictionary[residue1]
            scheme2 = self.labeling_dictionary[residue2]
            code = self._calculate_code_multiple(scheme1, scheme2)
            meaning = "".join([residue1, str(i+1), " - ", residue2, str(i+2)])
            self.codes_table.append(code)
            self.meanings_table.append(meaning)
            if code not in self.code_meanings:
                self.code_meanings[code] = [pair]
            else:
                self.code_meanings[code].append(pair)
        LI_a = 0
        LI_p = 0
        LN1_a = 0
        LN1_p = 0
        LN2_a = 0
        LN2_p = 0
        LU_a = 0
        LA_p = 0
        LA_a = 0

        for code in self.code_meanings:
            meanings = self.code_meanings[code]
            if code == "0" * self.samples_num:
                LI_p += len(set(meanings))
                LI_a += len(meanings)
            elif len(meanings) == 1:
                LU_a += 1
            elif len(set(meanings)) == 1:
                LN2_p += 1
                LN2_a += len(meanings)
            elif len(set([pair[1] for pair in meanings])) == 1: # all second residues in pair are equal
                LN1_p += len(set(meanings))
                LN1_a += len(meanings)
            else:
                LA_a += len(meanings)
                LA_p += len(set(meanings))
        LU_p = LU_a

        self.stats.update({
            "LIa": LI_a,
            "LIp": LI_p,
            "LUa": LU_a,
            "LUp": LU_p,
            "LN1a": LN1_a,
            "LN1p": LN1_p,
            "LN2a": LN2_a,
            "LN2p": LN2_p,
            "LAa": LA_a,
            "LAp": LA_p
        })

        self.ranked_meanings = [mean for code, mean in sorted(zip(self.codes_table,
                                                                self.meanings_table))]
        self.ranked_codes_table = [code for code in sorted(self.codes_table)]
        with open(filename, mode) as f:
            f.write("#"*50+"\n")
            f.write("# The Spectrum codes dictionary for \'"+self.name+"\':\n")
            f.write("# Spectrum code: First AA - Second AA\n\n")
            for i in range(len(self.ranked_codes_table)):
                f.write(": ".join([self.ranked_codes_table[i], self.ranked_meanings[i]]) + "\n")
            f.flush()
            f.close()

    def write_solution_to_file(self, filename, comment=True, mode='w'):

        output=""
        if comment:
            output += "#"*50+"\n"
            output += "# The list of solutions found        \n"
            output += "#\n\n"
        output += "[solution]\n"

        if self.optimize_price:
            output += "% Solution number = " + str(self.solution_number) + "\n"
            output += "% Solution price  = " + "%.2f" % self.calculate_total_price() + "\n"
        output += "Res"
        for i in range(self.samples_num):
            output += ", S_" + str(i+1)
        output += "\n"
        for i in range(len(self.solution)):
            output += self.task.to_three_letter_code[self.sequence.residues_to_label[i]] + ",   " + ",   ".join(list(self.solution[i])) + "\n"
        for res in self.sequence.non_labeled_residues:
            output += self.task.to_three_letter_code[res] + ",   " + ",   ".join(list(self._generate_other_code())) + "\n"
        output += "\n"

        with open(filename, mode) as f:
            f.write(output)
            f.flush()
            f.close()

    def write_full_pairs_table(self, filename, mode='w'):
        output = CITATION
        output += "#"*50+"\n"
        output += "# Table of all amino acid pairs \n"
        output += "# \n"
        output += "# Number in the table represents how many times \n"
        output += "# the pair occurs in the sequence \n\n"
        output += "[full_pairs_table]\n"
        output+= "   ," + ",".join([self.task.to_three_letter_code[res] for res in self.sequence.residues_second]) + "\n"
        for i in range(len(self.sequence.residues_first)):
            output += self.task.to_three_letter_code[self.sequence.residues_first[i]]
            for j in range(len(self.sequence.residues_second)):
                # MAYBE REPLACE 0 BY SOME SYMBOL????
                output += "," + "{:>3}".format(self.sequence.all_residue_pairs[i][j])
            if i+1 < len(self.sequence.residues_first):
                output += "\n"
        with open(filename, mode) as f:
            f.write(output)
            f.flush()
            f.close()

    def write_pairs_table(self, filename, mode='w'):
        output = "\n\n"+"#"*50+"\n"
        output += "# Table of amino acid pairs used for \n"
        output += "# combinatorial labeling \n"
        output += "# \n"
        output += "# Number in the table represents how many times \n"
        output += "# the pair occurs in the sequence \n\n"
        output += "[pairs_table]\n"
        additional_output = "\n"
        if self.sequence.residues_carbon[-1] == "Other":
            output += "   "
        output += "Res,"
        if self.sequence.residues_nitro[-1] == "Other":
            output += ",".join([self.task.to_three_letter_code[res] for res in self.sequence.residues_nitro[:-1]])
            output += ",OtherN"
            additional_output += "\nOtherN: " + ",".join([self.task.to_three_letter_code[res] for res in self.sequence.residues_not_nitro])
        else:
            output += ",".join([self.task.to_three_letter_code[res] for res in self.sequence.residues_nitro])
        output += "\n"
        for i in range(len(self.sequence.residues_carbon)):
            res1 = self.sequence.residues_carbon[i]
            if res1 == "Other":
                output += "OtherC"
                additional_output += "\nOtherC: " + ",".join([self.task.to_three_letter_code[res] for res in self.sequence.residues_not_carbon])
            else:
                if self.sequence.residues_nitro[-1] == "Other":
                    output += "   "
                output += self.task.to_three_letter_code[res1]
            for j in range(len(self.sequence.residues_nitro)):
                # MAYBE REPLACE 0 BY SOME SYMBOL????
                output += "," + "{:>3}".format(self.sequence.residue_pairs[i][j])
            if i + 1 < len(self.sequence.residues_carbon):
                output += "\n"
        output += additional_output + "\n"
        with open(filename, mode) as f:
            f.write(output)
            f.close()

    def write_pairs_codes(self, filename, mode='w'):
        output  = "\n\n"+"#"*50+"\n"
        output += "# Spectrum codes of the labeled amino acid pairs \n#\n"
        output += "# Amino acid and labeling code strings \n"
        output += "# according to the number of samples are in the headers\n"
        output += "# Spectrum codes are in the table\n\n"
        output += "[pairs_codes]\n"
        additional_output = "\n"
        separator = ","
        if self.samples_num > 3:
            separator += " " * (self.samples_num - 3)
        otherC_spaces = ""
        if self.sequence.residues_carbon[-1] == "Other":
            otherC_spaces += "   "
        output += otherC_spaces + "   ," + " " * self.samples_num + separator
        if self.sequence.residues_nitro[-1] == "Other":
            res_list = self.sequence.residues_nitro[:-1]
        else:
            res_list = self.sequence.residues_nitro
        output += separator.join([self.task.to_three_letter_code[res] for res in res_list])
        if self.sequence.residues_nitro[-1] == "Other":
            output += ",OtherN"
            additional_output += "\nOtherN: " + ",".join([self.task.to_three_letter_code[res] for res in self.sequence.residues_not_nitro])
        output += "\n"
        output += otherC_spaces + "   ," + " " * self.samples_num
        for res in res_list:
            output += "," + self.labeling_dictionary[res]
        if self.sequence.residues_nitro[-1] == "Other":
            output += "," + "X" * self.samples_num
        output += "\n"
        for i in range(len(self.sequence.residues_carbon)):
            res1 = self.sequence.residues_carbon[i]
            if res1 == "Other":
                additional_output += "\nOtherC: " + ",".join([self.task.to_three_letter_code[res] for res in self.sequence.residues_not_carbon])
                output += "OtherC"
            else:
                output += otherC_spaces + self.task.to_three_letter_code[res1]
            output += "," + self.labeling_dictionary[res1]
            for j in range(len(self.sequence.residues_nitro)):
                res2 = self.sequence.residues_nitro[j]
                if self.sequence.residue_pairs[i][j]:
                    code = self._calculate_code_multiple(self.labeling_dictionary[res1],
                                                         self.labeling_dictionary[res2])
                else:
                    code = " " * self.samples_num
                output += "," + code
            if i + 1 < len(self.sequence.residues_carbon):
                output += "\n"
        output += additional_output
        output += "\n"
        with open(filename, mode) as f:
            f.write(output)
            f.flush()
            f.close()

    def write_stats(self, filename, mode='w'):
        output = "\n\n"+"#"*50+"\n"
        output += "# Calculation statistics\n"
        output += "# \n\n"
        output += "[stats]\n"
        output += ""
        output += "\n\n"
        self.calculate_stats()
        output += "# Statistics for PAIRS in amino acid sequence\n"
        output += "# The STOCK (availability of isotopically labeled amino acid)\n"
        output += "# is NOT accounted for in this statistics\n"
        output += "# The labeling scheme is NOT accounted too\n"
        output += "\n"
        output += "[stats,pairs]\n"
        output += "Par, Description,  Residues, Pairs\n"
        output += "N,   Total number, {:>8}, {:>5}\n".format(self.stats["Na"], self.stats["Np"])
        output += "PI,  Invisible,    {:>8}, {:>5}\n".format(self.stats["PIa"], self.stats["PIp"])
        output += "PU,  Unique,       {:>8}, {:>5}\n".format(self.stats["PUa"], self.stats["PUp"])
        output += "PN,  Non-unique,   {:>8}, {:>5}\n".format(self.stats["PNa"], self.stats["PNp"])
        output += "\n"
        output += "# Statistics for STOCK-available pairs in amino acid sequence\n"
        output += "# The STOCK is used to check whether the particular pairs are distinguishable \n"
        output += "# in principle with some labeling scheme unlimited in size with some NMR spectra\n"
        output += "# The particular labeling scheme, found by the program, is NOT accounted here\n"
        output += "\n"
        output += "[stats,stock]\n"
        output += "Par, Description,     Residues, Pairs\n"
        output += "N,   Total number,    {:>8}, {:>5}\n".format(self.stats["Na"], self.stats["Np"])
        output += "SI,  Invisible,       {:>8}, {:>5}\n".format(self.stats["SIa"], self.stats["SIp"])
        output += "SU,  Unique code,     {:>8}, {:>5}\n".format(self.stats["SUa"], self.stats["SUp"])
        output += "SN2, AA type of both, {:>8}, {:>5}\n".format(self.stats["SN2a"], self.stats["SN2p"])
        output += "SN1, AA type of last, {:>8}, {:>5}\n".format(self.stats["SN1a"], self.stats["SN1p"])
        output += "\n"
        output += "# Statistics for LABELING CODES\n"
        output += "# The pairs are distinguishable, if their labeling codes are different\n"
        output += "# Both sequence, stock, NMR spectra and particular labeling scheme is accounted here\n"
        output += "\n"
        output += "[stats,labeling]\n"
        output += "Par, Description,     Residues, Pairs\n"
        output += "N,   Total number,    {:>8}, {:>5}\n".format(self.stats["Na"], self.stats["Np"])
        output += "LI,  Invisible,       {:>8}, {:>5}\n".format(self.stats["LIa"], self.stats["LIp"])
        output += "LU,  Unique code,     {:>8}, {:>5}\n".format(self.stats["LUa"], self.stats["LUp"])
        output += "LN2, AA type of both, {:>8}, {:>5}\n".format(self.stats["LN2a"], self.stats["LN2p"])
        output += "LN1, AA type of last, {:>8}, {:>5}\n".format(self.stats["LN1a"], self.stats["LN1p"])
        output += "LA,  Ambiguous code,  {:>8}, {:>5}\n".format(self.stats["LAa"], self.stats["LAp"])

        with open(filename, mode) as f:
            f.write(output)
            f.close()

    def write_results(self, filename, mode='w'):
        self.write_full_pairs_table(filename, mode)
        self.write_pairs_table(filename, mode='a')
        self.coding_table.write_codes_to_file(filename, mode='a')
        self.coding_table.write_pairs_to_file(filename, mode='a')
        if self.found:
            self.write_solution_to_file(filename, comment=True, mode='a')
            self.write_pairs_codes(filename, mode='a')
            self.write_stats(filename, mode='a')

    def calculate_stats(self):
        N_aa = len(self.sequence.sequence)
        N_pairs = len(self.sequence.unique_pairs)
        PI_p = N_pairs - len([pair for pair in self.sequence.unique_pairs if pair[1] != "P"])
        PI_a = N_aa - len([res for res in self.sequence.sequence if res != 'P'])
        PU_p = 0
        PN_a = 0
        PN_p = 0
        for i in range(len(self.sequence.residues_first)):
            for j in range(len(self.sequence.residues_second)):
                second_res = self.task.res_types
                number_of_pairs = self.sequence.all_residue_pairs[i][j]
                if second_res != "P" and number_of_pairs:
                    if number_of_pairs == 1:
                        PU_p += 1
                    else:
                        PN_p += 1
                        PN_a += number_of_pairs
        PU_a = PU_p
        SI_a = 0
        SI_p = 0
        SU_p = 0
        SN2_a = 0
        SN2_p = 0
        SN1_a = 0
        SN1_p = 0
        #print(len(self.sequence.residues_carbon))
        #print(self.sequence.residues_carbon)
        #print(len(self.sequence.residues_nitro))
        #print(self.sequence.residues_nitro)
        #print(self.sequence.residue_pairs)
        for i in range(len(self.sequence.residues_carbon)):
            for j in range(len(self.sequence.residues_nitro)):
                #print("i=%2d j=%2d" % (i,j))
                cell_value = self.sequence.residue_pairs[i][j]
                res1 = self.sequence.residues_carbon[i]
                res2 = self.sequence.residues_nitro[j]
                if cell_value:
                    seq = self.sequence.sequence
                    if res2 == 'Other':
                        diff_pairs = len(set([seq[s+1] for s in range(len(seq) - 1)
                                              if seq[s] == res1
                                              and seq[s+1] in self.sequence.residues_not_nitro]))
                        if res1 == 'Other':
                            diff_pairs = len(set([seq[s + 1] for s in range(len(seq) - 1)
                                                  if seq[s] in self.sequence.residues_not_carbon
                                                  and seq[s + 1] in self.sequence.residues_not_nitro]))
                        SI_a += cell_value
                        SI_p += diff_pairs
                    elif res1 == 'Other':
                        diff_pairs = len(set([seq[s] for s in range(len(seq) - 1)
                                              if seq[s + 1] == res2
                                              and seq[s] in self.sequence.residues_not_carbon]))
                        if res2 in self.sequence.bad_residues:
                            SN1_p += diff_pairs
                            SN1_a += cell_value
                        elif cell_value == 1:
                            SU_p += 1
                        elif diff_pairs > 1:
                            SN1_p += diff_pairs
                            SN1_a += cell_value
                        else:
                            SN2_p += 1
                            SN2_a += cell_value
                    elif (res1 == res2
                          and res2 in self.sequence.bad_residues):
                        SN1_p += 1
                        SN1_a += cell_value
                    elif cell_value == 1:
                        SU_p += 1
                    else:
                        SN2_p += 1
                        SN2_a += cell_value
        SU_a = SU_p
        self.stats.update({
            "SIa": SI_a,
            "SIp": SI_p,
            "SUa": SU_a,
            "SUp": SU_p,
            "SN1a": SN1_a,
            "SN1p": SN1_p,
            "SN2a": SN2_a,
            "SN2p": SN2_p,
            "Na": N_aa,
            "Np": N_pairs,
            "PIa": PI_a,
            "PIp": PI_p,
            "PUa": PU_a,
            "PUp": PU_p,
            "PNa": PN_a,
            "PNp": PN_p
        })


class Stock:
    '''
Stock object

Contains the amino acid stock
and also prices for each type of labeling
'''

    def __init__(self, task):
        self.task = task
        self.label_dict = {}
        self.label_options = {}
        self.label_num_dict = {}
        self.price_dict = {}
        self.usage = {}
        self.label_types = []
        self.price_label_types = []
        for res in self.task.res_types:
            self.usage.update({res: 1})

    def read_stock(self, lines):
        d = []
        first_entry = True
        row_len = 0
        for line in lines:
            if line[0] != "#" and line != "":
                s = [x.strip() for x in line.split(",")]
                if first_entry:
                    row_len = len(s)
                    first_entry = False
                elif len(s) != row_len:
                    message = "ERROR in stock file!\nThe lenths of rows are not equal"
                    raise err.ReadConfigError(msg=message)
                d.append(s)

        try:
            self.label_types = d[0][1:]
            for label_type in self.label_types:
                if label_type not in self.task.label_types:
                    raise err.ReadConfigError(msg="ERROR! Wrong labeling type in stock file.")
            for i in range(len(d) - 1):
                res_name = d[i + 1][0]
                if res_name in self.task.res_types:
                    res = res_name
                elif res_name.upper() in [r.upper() for r in self.task.res_types_3]:
                    res = self.task.to_one_letter_code[res_name]
                else:
                    raise err.ReadConfigError(msg="ERROR! Wrong residue name in stock file.")
                self.label_dict.update({res: ''.join(d[i + 1][1:])})
            self._generate_label_options()
        except IndexError:
            raise err.ReadConfigError(msg="ERROR in stock file.\nPlease check the length of each row.\nUse comma as separator")
        except err.ReadConfigError as e:
            raise err.ReadConfigError(msg=e.msg)

    def read_from_table(self, table):
        #self.label_types = ['X', 'N', 'C', 'D']
        for i in range(len(self.task.res_types)):
            res = self.task.res_types[i]
            table_part = [table[j][i] for j in range(len(table))]
            self.label_dict.update({res: ''.join(["1" if cell else "0" for cell in table_part])})
        self._generate_label_options()

    def read_prices_from_table(self, table):
        pass

    def _generate_label_options(self):
        for residue in self.task.res_types:
            if residue in self.label_dict:
                stock_to_change = self.label_dict[residue]
                option = []
                ## Labeling types are in preferred order,
                ##but due to success with this software the order doesn't matter.
                for label in self.task.label_types:
                    if label in self.label_types and label in self.task.coding_table.label_types_list:
                        label_index = self.label_types.index(label)
                        if stock_to_change[label_index] == '1':
                            option.append(label)
                label_num_dict = {}
                for i in range(len(option)):
                    label = option[i]
                    label_num_dict[label] = i
                self.label_options.update({residue: option})
                self.label_num_dict.update({residue: label_num_dict})

    def read_prices(self, lines):
        d = []
        for line in lines:
            a = line[0]
            if line[0] != "#" and line != "":
                s = [x.strip() for x in line.split(",")]
                d.append(s)
        self.price_label_types = d[0][1:]
        for label in self.label_types:
            if label not in self.price_label_types:
                raise err.ReadConfigError(msg="ERROR! Prices are not specified for '" + label + "' label type")

        # ADD check if file format is correct
        try:
            for i in range(len(d) - 1):
                curr_dict = {}
                for j in range(len(self.price_label_types)):
                    try:
                        price = float(d[i + 1][j + 1])
                    except ValueError:
                        raise err.ReadConfigError(msg="ERROR! Price must be set in digits\nPlease check price file (row {}; col {})".format(i+2, j+2))
                    curr_dict.update({self.price_label_types[j]: price})
                res_name = d[i + 1][0]
                if res_name in self.task.res_types:
                    residue_type = res_name
                elif res_name.upper() in [r.upper() for r in self.task.res_types_3]:
                    residue_type = self.task.to_one_letter_code[res_name]
                else:
                    raise err.ReadConfigError(msg="ERROR! Wrong residue name in price file.")
                for label_type in self.label_options[residue_type]:
                    if curr_dict[label_type] < 0:
                        raise err.ReadConfigError(msg="ERROR! Price is not specified or negative for '"
                              + residue_type + "'-residue's '"
                              + label_type + "' label type")
                self.price_dict.update({residue_type: curr_dict})
        except IndexError:
            raise err.ReadConfigError(msg="ERROR in price file.\nPlease check the length of each row.\nUse comma as separator")
        for res in self.label_options:
            for option in self.label_options[res]:
                try:
                    self.price_dict[res][option]
                except KeyError:
                    raise err.ReadConfigError(msg="ERROR in price file.\n"
                                                  "Price is not specified for {} label "
                                                  "option of {}".format(option, self.task.to_three_letter_code[res]))

    def write_to_file(self, filename, mode='w'):

        #ADD check if file exists

        f = open(filename, mode)
        f.write('Res\t' + '\t'.join(self.task.label_types) + '\n')
        for residue in self.task.res_types:
            line = residue
            for label in self.task.label_types:
                line += '\t'
                if residue in self.price_dict:
                    line += str(self.price_dict[residue][label])
                else:
                    line += '0'
            line += '\n'
            f.write(line)
        f.flush()
        f.close()

    def log_print(self):
        i = 1
        for key in self.label_dict:
            self.task.logfile.write (str(i) + "." + key + ":" + self.label_dict[key])
            i += 1
        self.task.logfile.write("\n")
        self.task.logfile.flush()


class Task:

    """
    Parent object that reads the config file, creates Solution, Stock, CodingTable and Solution
    objects and also performs main search cycle in "solve" method.
    """

    def __init__(self, args):
        self.start_time=time.time();
        self.res_types = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
             "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
        self.res_types_3 = ("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
                          "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
        self.to_one_letter_code = {
            "Ala":"A", "Cys":"C",
            "Asp":"D", "Glu":"E",
            "Phe":"F", "Gly":"G",
            "His":"H", "Ile":"I",
            "Lys":"K", "Leu":"L",
            "Met":"M", "Asn":"N",
            "Pro":"P", "Gln":"Q",
            "Arg":"R", "Ser":"S",
            "Thr":"T", "Val":"V",
            "Trp":"W", "Tyr":"Y"
        }
        self.to_three_letter_code = {}
        for code3 in self.to_one_letter_code:
            code1 = self.to_one_letter_code[code3]
            self.to_three_letter_code[code1] = code3
        self.label_types = ("X", "N", "C", "D", "A", "S", "T", "F")
        self.cn_label_types = ("D", "S", "T")
        self.n_label_types = ("N", "D", "S", "T")
        self.c_lable_types = ("C", "D", "A", "S", "T", "F")
        self.sequence = Sequence(self)
        self.stock = Stock(self)
        self.optimize_price = False
        self.spectra_vector = []
        self.labeling_types_vector = []
        self.start_samples = 0
        self.set = False
        self.name = ""
        self.sequence_ok = False
        self.stock_ok = False
        self.price_ok = False
        self.spectra_ok = False
        self.label_types_ok = False
        self.check_file = args.check
        self.check = True
        if self.check_file == '-':
            self.check = False
        else:
            self.check_number = int(args.number)
        try:
            self.read_config(args.config)
        except err.ReadConfigError as e:
            raise err.ReadConfigError(msg=e.msg)
        if int(args.samples) > 1:
            self.start_samples = int(args.samples) - 1
        self.create_solution()

    def read_config(self, config_file):
        try:
            lines = self.get_lines_from_file(config_file)
        except err.ReadLinesError as e:
            raise err.ReadConfigError(msg="ERROR! Config file '" + e.filename + "' not found")
        parameters = {}
        for line in lines:
            split_comment = line.split("#")[0]
            split_line = split_comment.split()
            try:
                parameters[split_line[0]] = split_line[1:]
            except IndexError:
                pass

        try:
            job_name = parameters["job_name"]
        except KeyError:
            raise err.ReadConfigError(msg="ERROR! Job name is not specified in the config file")
        if job_name == []:
            raise err.ReadConfigError(msg="ERROR! Job name is not specified in the config file")
        self.name = job_name[0]

        try:
            spec_vec = parameters["spectra_types"]
        except KeyError:
            raise err.ReadConfigError(msg="ERROR! Spectra types are not specified in the config file")
        if spec_vec == []:
            raise err.ReadConfigError(msg="ERROR! Spectra types are not specified in the config file")
        try:
            self.read_spec_vec(spec_vec)
        except err.ReadConfigError as e:
            raise err.ReadConfigError(msg=e.msg)

        if not self.check:
            try:
                lab_vec = parameters["labeling_types"]
            except KeyError:
                raise err.ReadConfigError(msg="ERROR! Labeling types are not specified in the config file")
            if lab_vec == []:
                raise err.ReadConfigError(msg="ERROR! Labeling types are not specified in the config file")
            try:
                self.read_label_vec(lab_vec)
            except err.ReadConfigError as e:
                raise err.ReadConfigError(msg=e.msg)
        else:
            try:
                self.read_solution(self.check_file)
            except err.ReadConfigError as e:
                raise err.ReadConfigError(msg=e.msg)

        try:
            self.make_coding_table()
        except err.LabelPowerError as e:
            raise err.ReadConfigError(msg=e.msg)

        if not self.check:
            try:
                stock_file = parameters["stock_file"]
            except KeyError:
                raise err.ReadConfigError(msg="ERROR! Stock file is not specified in the config file")
            if stock_file == []:
                raise err.ReadConfigError(msg="ERROR! Stock file is not specified in the config file")
            try:
                self.read_stock(stock_file[0])
            except err.ReadConfigError as e:
                raise err.ReadConfigError(msg=e.msg)

        else:
            self.stock.label_options = self.label_options
            self.stock.label_num_dict = self.label_num_dict

        self.coding_table.recalculate_coding_table()

        try:
            optimize_price = parameters["optimize_price"]
        except KeyError:
            raise err.ReadConfigError(msg="ERROR! Price optimization is not specified")
        if optimize_price[0] == "Yes":
            self.set_price_optimization(True)
            try:
                price_file = parameters["prices_file"]
            except KeyError:
                raise err.ReadConfigError(msg="ERROR! Price file is not specified in the config file")
            if price_file == []:
                raise err.ReadConfigError(msg="ERROR! Price file is not specified in the config file")
            try:
                self.read_prices(price_file[0])
            except err.ReadConfigError as e:
                raise err.ReadConfigError(msg=e.msg)
        elif optimize_price[0] == "No":
            self.set_price_optimization(False)
        else:
            raise err.ReadConfigError(msg="ERROR! 'optimize_price' field must have value 'Yes' or 'No'")

        try:
            sequence_file = parameters["sequence_file"]
        except KeyError:
            raise err.ReadConfigError(msg="ERROR! Sequence file is not specified in the config file")
        if sequence_file == []:
            raise err.ReadConfigError(msg="ERROR! Sequence file is not specified in the config file")
        try:
            self.read_sequence(sequence_file[0])
        except err.ReadConfigError as e:
            raise err.ReadConfigError(msg=e.msg)

        try:
            self.start_samples = int(parameters["start_samples"][0])
            self.start_samples -= 1
            if self.start_samples < 0:
                self.start_samples = 0
        except KeyError:
            self.start_samples = 0

    def get_lines_from_file(self, filename):
        try:
            f = open(filename, 'r', encoding='UTF-8')
            lines = f.readlines()
            f.close()
        except IOError:
            raise err.ReadLinesError(filename)
        new_lines = []
        for line in lines:
            curr_line = line.rstrip()
            if curr_line:
                new_lines.append(curr_line)
        return new_lines

    def import_sequence(self, sequence):
        self.sequence.sequence = sequence
        self.sequence_ok = True

    def read_stock(self, stock_file):
        try:
            lines = self.get_lines_from_file(stock_file)
        except err.ReadLinesError as e:
            raise err.ReadConfigError(msg="ERROR! Stock file '" + e.filename + "' not found")
        try:
            self.stock.read_stock(lines)
        except err.ReadConfigError as e:
            raise err.ReadConfigError(msg=e.msg)
        self.stock_ok = True

    def read_spec_vec(self, spec_vec):
        if len(spec_vec) != 7:
            raise err.ReadConfigError(msg="ERROR! Not for all (7) spectra the usage is specified\nor specified for more than 7")
        for item in spec_vec:
            if item != "1" and item != "0":
                raise err.ReadConfigError(msg="ERROR! Only '1' or '0' must be used to specify spectra")
        self.spectra_vector = [int(item) for item in spec_vec]
        self.spectra_ok = True

    def read_label_vec(self, lab_vec):
        if len(lab_vec) != 8:
            raise err.ReadConfigError(msg="ERROR! Not for all (8) labeling types specified\nor specified for more than 8")
        for item in lab_vec:
            if item != "1" and item != "0":
                raise err.ReadConfigError(msg="ERROR! Only '1' or '0' must be used to specify labeling types")
        self.labeling_types_vector = [int(item) for item in lab_vec]
        self.label_types_ok = True

    def read_solution(self, solution_file):
        try:
            lines = self.get_lines_from_file(solution_file)
        except err.ReadLinesError as e:
            raise err.ReadConfigError(msg="ERROR! Solution file '" + e.filename + "' not found")
        start_found = False
        solution_found = False
        self.got_solution = {}
        self.labeling_types_vector = [0 for _ in range(len(self.label_types))]
        self.label_options = {}
        number = 0
        for line in lines:
            if not start_found and "[solution]" in line:
                number += 1
                if number == self.check_number:
                    start_found = True
            if start_found:
                entries = re.split(", *", line)
                if entries[0].upper() in self.res_types \
                        or entries[0].upper() in [res.upper() for res in self.res_types_3]:
                    solution_found = True
                    if entries[0].upper() in self.res_types:
                        self.got_solution[entries[0]] = "".join(entries[1:])
                    elif entries[0].upper() in [res.upper() for res in self.res_types_3]:
                        residue = self.to_one_letter_code[entries[0]]
                        self.got_solution[residue] = "".join(entries[1:])
                elif solution_found:
                    break
        if not start_found or not solution_found:
            if self.check_number > 1:
                raise err.ReadConfigError(msg="ERROR! In the file '{}' there are less then {} solutions".format(solution_file, self.check_number))
            else:
                raise err.ReadConfigError(msg="ERROR! No solution found in solution file '{}'".format(solution_file))
        self.check_samples = 0
        first = True
        for residue in self.got_solution:
            value = self.got_solution[residue]
            if first:
                self.check_samples = len(value)
                first = False
            else:
                if len(value) != self.check_samples:
                    raise err.ReadConfigError(msg="ERROR! Not equal number of samples for residues in solution file '{}'".format(solution_file))
            options = []
            for type in self.label_types:
                if type in value:
                    if type not in options:
                        options.append(type)
                    self.labeling_types_vector[self.label_types.index(type)] = 1
            self.label_options[residue] = options
        self.got_solution['Other'] = 'X' * self.check_samples
        for residue in self.res_types:
            if residue not in self.label_options:
                self.label_options[residue] = ['X']
        self.label_num_dict ={}
        for i in range(len(self.label_types)):
            self.label_num_dict[self.label_types[i]] = i

    def read_sequence(self, seq_file):
        try:
            lines = self.get_lines_from_file(seq_file)
        except err.ReadLinesError as e:
            raise err.ReadConfigError(msg="ERROR! Sequence file '" + e.filename + "' not found")
        try:
            self.sequence.read_from_file(lines)
        except err.ReadConfigError as e:
            raise err.ReadConfigError(msg=e.msg)
        self.sequence_ok = True

    def set_price_optimization(self, value):
        self.optimize_price = value
        self.price_ok = not value

    def read_prices(self, price_file):
        try:
            lines = self.get_lines_from_file(price_file)
        except err.ReadLinesError as e:
            raise err.ReadConfigError(msg="ERROR! Price file '" + e.filename + "' not found")
        try:
            self.stock.read_prices(lines)
        except err.ReadConfigError as e:
            raise err.ReadConfigError(msg=e.msg)
        self.price_ok = True

    def make_coding_table(self):
        try:
            self.coding_table = CodingTable(self)
        except err.LabelPowerError:
            raise err.LabelPowerError

    def create_solution(self):
        self.sequence.calculate_stats(self.stock)
        if self.optimize_price:
            self.sequence.predict_prices()
        self.solution = Solution(self)
        self.best_solution = Solution(self)
        self.solution.copy_solution_to(self.best_solution)

    def solve(self):
        if self.sequence_ok and self.stock_ok and self.price_ok and self.spectra_ok and self.label_types_ok:
            t0 = time.time()
            samples_number = 0
            with open(self.name + "_logfile.log", "w") as f:
                log = "\nSearch for {} solution started\n".format(self.name)
                log += "Host name is                      {}\n".format(socket.gethostname())
                log += "Start clock time:                 {}\n".format(time.strftime("%d-%m-%Y %H:%M:%S", time.gmtime(t0)))
                f.write(log)
                f.flush()
                f.close()
                print(log)
            self.iteration = 0

            solutions = 0
            final_depth = len(self.solution.sequence.residues_to_label)
            print("\nSearch for {} solution started".format(self.name))
            print("Host name is:                     {}".format(socket.gethostname()))
            print("Start clock time:                 {}".format(time.strftime("%d-%m-%Y %H:%M:%S", time.gmtime(t0))))

            self.logfile = open(self.name + "_logfile.log", "a")
            self.write_results()

            while True:

                #   MAIN CYCLE, increments Solution and checks it

                self.iteration += 1
                self.solution.increment_state()
                self.solution.check_and_update()

                # Check if solution is found
                if self.solution.depth == final_depth and self.solution.good:

                    # block that works if we optimize the solution by price
                    if self.optimize_price:
                        if not self.solution.found:
                            self.solution.solution_number += 1
                            samples_number = self.solution.samples_num  # remember samples number
                            self.solution.best_price = self.solution.calculate_total_price()
                            current_output = "\nSolutions found: {}\nBest price: {}\nIterations: {}".format(
                                self.solution.solution_number,
                                "%.2f" % self.solution.best_price,
                                self.iteration)
                            print(current_output)
                            self.logfile.write(self.log_time(t0) + current_output + "\n")
                            self.logfile.flush()
                            self.solution.copy_solution_to(self.best_solution)
                            self.solution.write_solution_to_file(self.name + "_all_solutions.txt")
                            self.solution.found = True
                            self.write_results()
                        if self.solution.best_price > self.solution.calculate_total_price():
                            self.solution.solution_number += 1
                            self.solution.copy_solution_to(self.best_solution)
                            self.solution.write_solution_to_file(self.name + "_all_solutions.txt", mode='a', comment=False)
                            self.solution.best_price = self.solution.calculate_total_price()
                            current_output = "\nSolutions found: {}\nBest price: {}\nIterations: {}".format(
                                self.solution.solution_number,
                                "%.2f" % self.solution.best_price,
                                self.iteration)
                            print(current_output)
                            self.logfile.write(self.log_time(t0) + current_output+"\n")
                            self.logfile.flush()
                            self.write_results()
                    else:
                        self.solution.solution_number += 1
                        self.solution.found = True
                        self.solution.copy_solution_to(self.best_solution)
                        print("\nSolution found!")
                        print(self.solution)
                        break
                    solutions += 1

                # break the cycle if number of samples increased and the solution is already found
                # bad idea to continue the search anyways
                if self.solution.found and samples_number != self.solution.samples_num:
                    break

                # TESTING BLOCK

                if self.iteration % 100000 == 0:
                    current_output = ""
                    current_output += self.log_time(t0)
                    current_output += "\niteration: " + str(self.iteration)+"; "

                    current_output += "iterations/min: "+str(round(self.iteration * 60 / (time.time() - t0), 1))+"\n"
                    if self.optimize_price:
                        current_output += "curr best price: {}\n".format(self.solution.best_price)
                    print("\nJob name: {}".format(self.name) + current_output)
                    current_output += str(self.solution)+"\n"

                    self.logfile.write(current_output +"\n")
                    self.logfile.flush()

            # End of solution search
            # Print final output to the log file
            self.logfile.write("\nEvaluation Finished!\n")
            self.logfile.write("Final clock time:                " + time.strftime("%d-%m-%Y %H:%M:%S", time.gmtime())+"\n")
            self.logfile.write(self.log_time(t0,
                               "Final evaluation time DDHHMMSS:  ") + "\n")
            self.logfile.write("Final evaluation time (seconds): " + str(time.time() - t0) + "\n")
            self.logfile.write("Final number of samples:         " + str(samples_number) + "\n")
            self.logfile.write("Final number of iterations:      " + str(self.iteration)+"\n")
            self.logfile.write("Final number of codes evaluated: " + str(self.solution.calc_code) + "\n")
            self.logfile.close()
            self.best_solution.print_table()
            return self.iteration

    def log_time(self, t0, logtext="Elapsed time:"):
        t1 = time.time()
        t_diff = time.gmtime(t1 - t0)
        #output = "\nStart time: " + time.strftime("%d-%m-%Y %H:%M:%S", time.gmtime(t1)) + "\n"
        output = ("\n" + logtext + " {:>2} days {}").format(t_diff.tm_yday-1, time.strftime("%H:%M:%S", t_diff))
        return(output)

    def write_results(self):
        solution_filename = self.name + "_solution.txt"
        code_dict_filename = self.name + "_code_dictionary.txt"
        if self.solution.found:
            self.best_solution.write_codes_table(code_dict_filename)
        else:
            pass
        self.best_solution.write_results(solution_filename, mode='w')

            # solution_filename = self.name + "_solution.txt"
            # code_dict_filename = self.name + "_code_dictionary.txt"
            # codes_filename = self.name + "_codes.txt"
            # labeling_types_pairs_filename = self.name + "_labeling_pairs.txt"
            # pairs_table_filename = self.name + "_pairs.txt"
            # full_table_filename = self.name + "_all_pairs.txt"
            # pairs_codes_table_filename = self.name + "_pairs_codes.txt"
            # self.best_solution.write_to_file(solution_filename)
            # self.best_solution.write_codes_table(code_dict_filename)
            # self.coding_table.write_to_file(codes_filename)
            # self.coding_table.write_pairs_to_file(labeling_types_pairs_filename)
            # self.best_solution.write_full_pairs_table(full_table_filename)
            # self.best_solution.write_pairs_table(pairs_table_filename)
            # self.best_solution.write_pairs_codes(pairs_codes_table_filename)

    def check_solution(self):
        self.solution.solution = []
        self.solution.samples_num = self.check_samples
        self.solution.labeling_dictionary = self.got_solution
        for residue in self.sequence.residues_to_label:
            self.solution.solution.append(self.got_solution[residue])
        self.solution.depth = len(self.sequence.residues_to_label)
        self.solution.found = True
        self.best_solution = self.solution
        self.write_results()
        self.solution.print_table(True)

    def std_output(self, time):
        iterations_per_min = "%.0f" % (self.iteration / (time) * 60)
        regular_output = "\nIteration: {}\nIterations/min: {}".format(self.iteration, iterations_per_min)
        if self.solution.found:
            regular_output += "\nSolutions found: {}\nBest price: {}".format(self.solution.solution_number,
                                                                             "%.2f" % self.solution.best_price)
        else:
            regular_output += "\nMax depth reached: {}".format(self.solution.max_depth)
        print(regular_output)
        self.logfile.write("\n" + regular_output + "\n")
        self.logfile.write(str(self.solution))
        self.logfile.flush()


class CodingTable:

    """
    CodingTable object contains the table of codes.
    Main variable is codes_dict - a dictionary that returns code for given
    first and second albeling type. This dictionary is calculated based on
    the set of spectra and label types in use. Also the codes_dict is recalculated
    if ticked label types in config file are not present in stock for all residues.
    Each labeling type is assigned "label_power" - the number of distinct codes
    that can be obtained if this label type is second in pair
    """


    def __init__(self, task):
        self.task = task
        self.label_vec = task.labeling_types_vector
        self.spec_vec = task.spectra_vector
        self.spec_types = ("HSQC", "HNCO", "HNCA", "HNCOCA", "CO-HNCA", "DQ-HNCA", "HNCACO")
        self.spec_list = [self._HSQC, self._HNCO, self._HNCA, self._HNCOCA, self._COHNCA, self._DQHNCA, self._HNCACO]
        self.label_meanings = {
            "X": "000",
            "N": "100",
            "C": "001",
            "D": "111",
            "A": "010",
            "S": "101",
            "T": "110",
            "F": "011",
        }
        self.label_types = self.task.label_types  # THIS order of labels sets their priority in search for solution
        self.letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']
        self.codes_dict = {}
        self.label_power = {}
        self.vectors = []
        self.spectra_numbers = []
        self.label_types_list = []
        try:
            self._make_coding_table()
        except err.LabelPowerError:
            raise err.LabelPowerError


    def _HSQC(self, atom_list):
        return int(atom_list[3])

    def _HNCO(self, atom_list):
        return int(atom_list[3] and atom_list[2])

    def _HNCA(self, atom_list):
        return int(atom_list[3] and (atom_list[1] or atom_list[4]))

    def _HNCOCA(self, atom_list):
        return int(atom_list[3] and atom_list[2] and atom_list[1])

    def _COHNCA(self, atom_list):
        return int(atom_list[3] and atom_list[1] and not atom_list[2])

    def _DQHNCA(self, atom_list):
        return int(atom_list[3] and atom_list[1] and atom_list[4])

    def _HNCACO(self, atom_list):
        return int(atom_list[3] and atom_list[4] and atom_list[5])

    def _make_coding_table(self):
        self.label_types_list = []
        for i in range(len(self.label_types)):
            if self.label_vec[i]:
                self.label_types_list.append(self.label_types[i])
        spec_types_list = []
        for i in range(len(self.spec_list)):
            if self.spec_vec[i]:
                spec_types_list.append(self.spec_list[i])
                self.spectra_numbers.append(i)
        codes_table = [[0 for i in range(len(self.label_types_list))] for j in range(len(self.label_types_list))]
        self.vectors = [[0 for _ in spec_types_list]]
        for i in range(len(self.label_types_list)):
            label_2 = self.label_types_list[i]
            for j in range(len(self.label_types_list)):
                label_1 = self.label_types_list[j]
                atom_list = self._make_atom_list(label_1, label_2)
                vector = []
                for spectrum in spec_types_list:
                    vector.append(spectrum(atom_list))
                if vector in self.vectors:
                    code = self.vectors.index(vector)
                else:
                    code = len(self.vectors)
                    self.vectors.append(vector)
                if code > 9:
                    codes_table[j][i] = self.letters[code-10]
                else:
                    codes_table[j][i] = str(code)
        for i in range(len(self.label_types_list)):
            result = []
            for row in codes_table:
                result.append(row[i])
            power = len(set(result))
            self.label_power[self.label_types_list[i]] = power
            if power == 1 and self.label_types_list[i] in self.task.n_label_types:
                raise err.LabelPowerError
        for i in range(len(self.label_types_list)):
            label_1 = self.label_types_list[i]
            subdict = {}
            for j in range(len(self.label_types_list)):
                label_2 = self.label_types_list[j]
                subdict[label_2] = codes_table[i][j]
            self.codes_dict[label_1] = subdict

    def _make_atom_list(self, first_type, second_type):
        atom_string = self.label_meanings[first_type] + self.label_meanings[second_type]
        return [int(symbol) for symbol in atom_string]

    def write_pairs_to_file(self, filename, mode='w'):
        output = "\n\n"+"#"*50+"\n"
        output+= "# Spectrum code for each labeling pair \n#\n"
        output+= "# One-letter codes in the headers of columns"
        output+= "# and rows are the labeling types \n"
        output+= "# Don't confuse with one-letter amino acid codes\n\n"
        output+= "[code_pairs]\n "
        for type in self.label_types_list:
            output += ", " + type
        output += "\n"
        for type_1 in self.label_types_list:
            output += type_1
            for type_2 in self.label_types_list:
                output += ", " + self.codes_dict[type_1][type_2]
            output += "\n"
        output += "\n"

        with open(filename, mode) as f:
            f.write(output)
            f.flush()
            f.close()

    def write_codes_to_file(self, filename, mode='w'):
        output = "\n\n"+"#"*50+"\n"
        output += "# Spectrum codes table\n#\n"
        output += "# The spectrum code is in the first column\n"
        output += "# Flag of the peak presence (0 or 1) for each spectrum\n\n"
        output += "[codes]\n"
        output += "Code"
        for i in range(len(self.vectors[0])):
            output += "," + "{:>7}".format(self.spec_types[self.spectra_numbers[i]])
        output += "\n"
        for i in range(len(self.vectors)):
            if i < 10:
                code = "{:>4}".format(i)
            else:
                code = "{:>4}".format(self.letters[i-10])
            output += code + ","
            output += ",".join(["{:>7}".format(item) for item in self.vectors[i]])
            if i+1 < len(self.vectors):
                output += "\n"

        output += "\n"
        with open(filename, mode) as f:
            f.write(output)
            f.flush()
            f.close()

    def recalculate_coding_table(self):
        # recalculate coding table to make sure, that in main algorithm
        # only labelings that are present in stock are used
        # (in case some labelings are ticked in config file and not present in stock)

        self.label_vec = []
        current_types = []
        for residue in self.task.res_types:
            for type in self.task.stock.label_options[residue]:
                if type not in current_types:
                    current_types.append(type)
        for type in self.label_types:
            if type in current_types:
                self.label_vec.append(1)
            else:
                self.label_vec.append(0)
        self._make_coding_table()
        # self.write_pairs_to_file(self.task.name+"_coding_pairs.txt")


class LabelList:

    """
    This object is generated for each residue in residues_to_label list
    for each number of samples in solution.
    Initial label_set is generated.
    Then at once the labelings with not enough 15N labels are crossed out from label_set 
    During the search for solution, labelings that knowingly will not lead to solution are
    crossed out from the label_set and stored in cross_out_sets.
    If the residue is being labeled, at first the symmetrical labelings are crossed-out
    and then the label_list is generated based on current label_set. the the main algorithm
    goes through all these labelings in the label_list using index variable.
    Each time current solution doesn't pass the checks all crossed out labelings are restored.  
    """

    def __init__(self, solution, residue_depth):
        self.solution = solution    # solution object
        self.residue_depth = residue_depth  # number of the residue in residues_to_label list
        # one-letter code for LabelList's residue
        self.residue = self.solution.sequence.residues_to_label[residue_depth]

        # labeling options for this particular residue
        self.label_options = self.solution.stock.label_options[self.residue]

        self.samples_number = self.solution.samples_num
        self.max_depth = len(self.solution.sequence.residues_to_label)

        # main labeling set, from which all cross-outs are made
        self.label_set = self._generate_initial_set(self.samples_number)

        # label list that is generated from current label_set when the residue is being labeled
        self.label_list = []

        # index of current labeling in the label_list when the residue is being labeled
        self.index = 0

        # label powers of each label type, e.a. how many
        # different codes it gives when it is second in pair
        self.label_power = self.solution.label_power

        # cross out all labelings that don't give enough info if residue has 15N labeling
        self._cross_out_N_power()

        # generate empty cross-out sets
        self.crossed_out_sets = [set() for _ in range(len(self.solution.sequence.residues_to_label))]

        # generate empty symmetry cross-out set
        self.symmetry_cross_out = set()

    def _generate_initial_set(self, samples):
        # recursive function, that generates all possible combinations
        # of labels given the number of samples for the given residue

        if samples == 0:
            new_set = set()
            new_set.add("")
            return new_set
        current_set = self._generate_initial_set(samples - 1)
        new_set = set()
        for item in current_set:
            for option in self.label_options:
                new_set.add(item + option)
        return new_set

    def _cross_out_N_power(self):
        # cross out all labelings that are knowingly will not give
        # enough info to distinguish all cells in the column

        if self.residue in self.solution.sequence.residues_nitro:
            self.N_cross_out = set()
            for labeling in self.label_set:
                if not self._check_N_power(labeling):
                    self.N_cross_out.add(labeling)
            self.label_set = self.label_set.difference(self.N_cross_out)

    def _check_N_power(self, labeling):
        # check if there is enough N labels in the labeling
        # coding theory is used

        max_pairs = 1
        got_nitro = False
        for label in labeling:
            max_pairs *= self.label_power[label]
            if label in self.solution.sequence.nitro_types:
                got_nitro = True
        nitro_index = self.solution.sequence.residues_nitro.index(self.residue)
        return got_nitro and max_pairs >= self.solution.sequence.min_nitrogens[nitro_index]

    def cross_out(self, label_set_to_cross):
        # main cross out method
        # each crossed-out label is stored in self.crossed_out_sets list
        # at index, corresponding to solution depth, so that each time
        # you can restore original labeling set at each depth

        if not self.is_frozen():
            crossed_out = self.label_set.intersection(label_set_to_cross)
            current_depth = self.solution.depth
            self.crossed_out_sets[current_depth-1] = self.crossed_out_sets[current_depth-1].union(crossed_out)
            self.label_set = self.label_set.difference(crossed_out)
            self.index = 0

    def cross_out_symmetry(self, symmetry):
        # cross out symmetrical labelings to skip similar solutions with swapped samples

        # crossed out labelings are stored in the following set, so you can restore
        self.symmetry_cross_out = set()

        if len(symmetry) == self.solution.samples_num:
            return
        residue = self.residue
        label_num_dict = self.solution.stock.label_num_dict[residue]
        for label in self.label_set:
            start_point = 0
            for i in range(len(symmetry)):
                if symmetry[i] > 1:
                    for j in range(symmetry[i] - 1):
                        if (label_num_dict[label[start_point + j]]
                                > label_num_dict[label[start_point + j + 1]]):
                            self.symmetry_cross_out.add(label)
                start_point += symmetry[i]
        self.label_set = self.label_set.difference(self.symmetry_cross_out)

    def restore_symmetry(self):
        # method to restore symmetry

        self.label_set = self.label_set.union(self.symmetry_cross_out)

    def restore_last_depth(self):
        # restore all crossed out labelings at the previous solution depth

        current_depth = self.solution.depth
        self.label_set = self.label_set.union(self.crossed_out_sets[current_depth-1])
        self.crossed_out_sets[current_depth-1] = set()
        self.index = 0

    def is_frozen(self):
        # "frozen" residues are those, that are already labeled and the one being labeled
        return self.solution.depth >= (self.residue_depth + 1)

    def increment_index(self):
        # increment index for residue being labeled to get the next labeling

        if self.is_frozen() and not self.has_last_index():
            self.index += 1

    def get_labeling(self):
        if self.is_frozen():
            return self.label_list[self.index]

    def update_on_depth_change(self):
        # method that creates list of labeling from the current labeling set
        # for residue that is the next to be labeled
        # After this main cycle goes through the list changing the index
        # and getting the next index

        if self.is_frozen():
            self.index = 0
            alphabet = {"X": 3, "N": 1, "C":2, "D":4, "A":5, "S":6, "T":7, "F":8}
            self.label_list = sorted(self.label_set, key=lambda word: [alphabet[c] for c in word])

    def has_last_index(self):
        return self.index+1 == len(self.label_list)

    def is_empty(self):
        if len(self.label_set):
            return False
        else:
            return True

    def predict_price(self):
        # find the price of cheapest labeling left in the labeling set and use it
        # as bottom estimate of price
        first = True
        best_price = 0
        for labeling in self.label_set:
            price = self.solution.calculate_price_for_residue(self.residue, labeling)
            if first:
                best_price = price
                first = False
            else:
                if price < best_price:
                    best_price = price
        return best_price


# FUNCTIONS

def main(args):
    try:
        task = Task(args)
    except err.ReadConfigError as e:
        print(e.msg)
        return
    if task.check:
        task.check_solution()
    else:
        task.solve()
        task.write_results()


if __name__ == '__main__':
    main(args)
