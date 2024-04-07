import re
import sys
import multiprocessing

header = ""
total_entry = ""
outputH = open(sys.argv[2], "w")

def annotatePolyA(fasta_entry, window_size = 10, A_richness = 0.8):
    seq = fasta_entry[1].upper()
    size = len(seq)
    subset = seq[size - window_size:]
    result = ""
    if (subset.count("A"))  / window_size > A_richness:
        end = size
        start = size - window_size
        while (subset.count("A")) / (size - start) > A_richness:
            start -= 1
            subset = seq[start:]
        if start > (0.5 * size):
            header = fasta_entry[0]
            infos = re.findall("([.\w]+):(\d+)-(\d+)\(([+-])\)", header)[0]
            if infos[3] == "-":
                end = int(infos[2]) - start
                start = infos[1]
            else:
                start = int(infos[1]) + start
                end = infos[2]
            result = "\t".join((infos[0], str(start), str(end), header, ".", infos[3] + "\n"))
    return(result)

entries = list()
total_entry = ""
for line in open(sys.argv[1]):
    if line.startswith(">"):
        if header != "":
            entries.append((header, total_entry))
        header = line.strip()[1:]
        print(header, end = "\r")
        total_entry = ""
    else:
        total_entry += line.strip()

entries.append((header, total_entry)) # add the last seq

with multiprocessing.Pool(processes = int(sys.argv[3])) as pool:
    results = pool.map(annotatePolyA, entries)
# results = [annotatePolyA(x) for x in entries]

outputH.writelines(results)
