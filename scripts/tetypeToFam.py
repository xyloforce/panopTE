import sys

full = False
if len(sys.argv) > 5:
    if sys.argv[5] == "full":
        full = True
    else:
        pass

species = ""
if len(sys.argv) > 6:
    species = sys.argv[6].replace("_", " ")

te_to_fam = dict()
for line in open(sys.argv[2]):
    line = line.strip().split("\t")
    # print(line)
    te_to_fam[line[0]] = (line[1], line[3])
    # create dict of TE type : fam, average_size

output_content = list()
bed_file = open(sys.argv[1]).readlines()
if full:
    repbase_dict = dict()
    for line in open(sys.argv[3]):
        line = line.strip().split("\t")
        if len(line) > 2:
            if line[2] == species:
                repbase_dict[line[0]] = line[-1]

    for line in bed_file:
        line = line.strip().split("\t")
        interval_size = int(line[2]) - int(line[1])
        try:
            ref_length = int(repbase_dict[line[3]])
        except KeyError:
            ref_length = float(te_to_fam[line[3]][-1]) # check against average size (bad but hey)
        interval_length = (ref_length - (0.1 * ref_length),
                           ref_length + (0.1 * ref_length))
        if (interval_size > interval_length[0]) & \
           (interval_size < interval_length[1]):
            line[3] = te_to_fam[line[3]][0]
            output_content.append("\t".join(line) + "\n")

else:
    for line in bed_file:
        line = line.strip().split("\t")
        line[3] = te_to_fam[line[3]][0]
        output_content.append("\t".join(line) + "\n")

output = open(sys.argv[4], "w")
output.writelines(output_content)
