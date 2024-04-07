import sys

te_to_fam = dict()
for line in open(sys.argv[2]):
    line = line.strip().split("\t")
    te_to_fam[line[0]] = line[2]

aggreg_counts = dict()
count = 0
for line in open(sys.argv[1]):
    line = line.split("\t")
    if te_to_fam[line[1]] not in aggreg_counts:
        aggreg_counts[te_to_fam[line[1]]] = dict()
    if line[0] not in aggreg_counts[te_to_fam[line[1]]]:
        aggreg_counts[te_to_fam[line[1]]][line[0]] = 0
    aggreg_counts[te_to_fam[line[1]]][line[0]] += 1
    count += 1
    if count % 10 == 0:
        print(count, end = "\r")


output_file = open(sys.argv[3], "w")
output_list = list()
for key in aggreg_counts:
    for pos in aggreg_counts[key]:
        output_list.append("\t".join((key, pos, str(aggreg_counts[key][pos]))))

output_file.writelines(output_list)