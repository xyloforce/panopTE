from fileinput import lineno
import re
import sys

te_fam = open(sys.argv[1])
te_species = open(sys.argv[2])
output = open(sys.argv[3], "w")

fam_to_equivalents = dict()
for line in te_fam:
    line = line.strip().split("\t")
    fam_to_equivalents["root;" + line[0]] = line[1:]

final_file = list()
for i in range(0,30): # skip first header lines
    te_species.readline()

for line in te_species:
    line = line.strip().split()
    length = re.findall("len=(\d+)", line[-1])[0] # formerly 3
    family = line[2]
    try:
        element_name = re.findall("'([\.\w-]+)'", line[1])[0]
    except IndexError:
        print(line[1])
        raise IndexError
    try:
        rm_type = fam_to_equivalents[family][4]
    except KeyError:
        print(line)
        raise KeyError
    count = fam_to_equivalents[family][3]
    try:
        rm_subtype = fam_to_equivalents[family][5]
    except IndexError:
        rm_subtype = rm_type
    final_file.append("\t".join((element_name, rm_type, rm_subtype, length, count)) + "\n")

output.writelines(final_file)
