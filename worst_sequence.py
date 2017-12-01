#!/usr/bin/python3

import argparse

parser = argparse.ArgumentParser(description='Generate worst sequence from selected amino acids')

parser.add_argument('-s', '--seq', default="ACDEFGHIKLMNPQRSTVWY")
args=parser.parse_args()
vars=vars(args)
SEQ=vars['seq']
#print(SEQ)

#RES_TYPES = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
#             "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

RES_TYPES=SEQ.strip('//')

pair_list = []
for res1 in RES_TYPES:
    for res2 in RES_TYPES:
        pair = res1 + res2
        pair_list.append(pair)

print("".join(pair_list))
