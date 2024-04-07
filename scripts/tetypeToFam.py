import sys

full = False
if len(sys.argv) > 5:
    if sys.argv[5] == "full":
        full = True
    else:
        pass

te_to_fam = dict()
for line in open(sys.argv[3]):
    line = line.strip().split("\t")
    te_to_fam[line[0]] = (line[2], line[3])

output_content = list()
bed_file = open(sys.argv[1]).readlines()
hit_file = open(sys.argv[2]).readlines()
for line, line_fam in zip(bed_file, hit_file[1:]): # skip first line of hitfile as its only headers
    line = line.split("\t")
    if full:
        line_fam = line_fam.split("\t")
        model_end = int(line_fam[7])
        fam_len = int(te_to_fam[line[3]][1])
        min_f = fam_len - (0.2 * fam_len)
        max_f = fam_len + (0.1 * fam_len)
        if min_f > model_end:
            line[3] = te_to_fam[line[3]][0]
            output_content.append("\t".join(line))
    else:
        line[3] = te_to_fam[line[3]][0]
        output_content.append("\t".join(line))

output = open(sys.argv[4], "w")
output.writelines(output_content)